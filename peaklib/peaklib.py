import itertools
import logging
import math
import os

from idr import idr
from scipy import signal
import numpy as np
import pandas as pd
import pybedtools

np.random.seed(0)

logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)
logger.setLevel(logging.DEBUG)
pybedtools.set_tempdir(os.path.abspath(os.path.dirname(__file__)))


class PeakCaller(object):
    def __init__(self, bedgraph, strand="+"):
        """
        Class to handle peak-calling using scipy.signal.find_peaks.
        """
        self.bedgraph = bedgraph
        self.df = pd.read_table(self.bedgraph, names=["chrom", "start", "stop", "val"])
        self.x = self.df["val"]
        self.strand = strand

    def call_peaks(self, prominence=(None, None), width=(1, None), rel_height=0.75):
        """
        Calls peaks per chromosome
        """

        self.peaks_by_chrom = {}
        self.meta_by_chrom = {}
        self.starts_by_chrom = {}
        self.chroms = []

        for chrom in sorted(self.df["chrom"].unique()):

            idx = self.df["chrom"] == chrom

            # TODO: should we expand this out to have all positions in the
            # chromosome in memory? Otherwise two, 1-bp peaks that are 10kb away
            # from each other will appear to the peak caller as adjacent. In
            # practice, at least using the existing data, this doesn't seem to
            # be a problem.
            x = self.df.loc[idx, "val"]
            starts = self.df.loc[idx, "start"]

            peaks, meta = signal.find_peaks(
                x, prominence=prominence, width=width, rel_height=rel_height
            )

            self.peaks_by_chrom[chrom] = peaks
            self.meta_by_chrom[chrom] = meta
            self.starts_by_chrom[chrom] = starts
            self.chroms.append(chrom)
        return self

    def peaks_to_bed(self):
        """
        Call peaks into the internal format, and then output as a narrowPeak
        file.

        Returns
        -------
        pybedtools.BedTool object sorted by score
        """
        logger.info(self.bedgraph)

        def gen():
            for chrom in self.chroms:
                logger.info(chrom)
                starts = self.starts_by_chrom[chrom]
                left_ips = self.meta_by_chrom[chrom]["left_ips"]
                right_ips = self.meta_by_chrom[chrom]["right_ips"]
                left_bases = self.meta_by_chrom[chrom]["left_bases"]
                right_bases = self.meta_by_chrom[chrom]["right_bases"]
                prominences = self.meta_by_chrom[chrom]["prominences"]
                widths = self.meta_by_chrom[chrom]["widths"]
                peaks = self.peaks_by_chrom[chrom]

                # Results from find_peaks are in the coordinate system of
                # integer indices, so we need to interpolate back out to
                # genomic coordinates.
                xp = np.arange(len(starts))
                ileft_ips = np.interp(left_ips, xp, starts).round().astype(int)
                iright_ips = np.interp(right_ips, xp, starts).round().astype(int)
                ipeaks = np.interp(peaks, xp, starts).round().astype(int)

                idx = ileft_ips <= iright_ips
                ileft_ips = ileft_ips[idx]
                iright_ips = iright_ips[idx]
                ipeaks = ipeaks[idx]
                widths = widths[idx]
                prominences = prominences[idx]

                n_removed = sum(~idx)
                if n_removed:
                    logger.info(
                        "Peaks removed due to start/stop problems: {0}".format(n_removed)
                    )

                for start, stop, peak, prominence, width in zip(
                    ileft_ips, iright_ips, ipeaks, prominences, widths
                ):
                    # This uses the promience as the score.
                    p = str(prominence)

                    # TODO: evaluate the usefulness of increasing the score for
                    # punctate peaks:
                    # p = str(prominence / math.sqrt(width))

                    yield pybedtools.create_interval_from_list(
                        [
                            chrom,
                            str(start),
                            str(stop),
                            ".",
                            p,
                            self.strand,
                            p,
                            "-1",
                            "-1",
                            str(peak - start),
                        ]
                    )

        # Ensure we're coord-sorted for the merging step
        x = pybedtools.BedTool(gen()).sort()
        x = merge_narrowbed(x, self.strand)

        # But the output needs to be sorted by score
        return sort_by_score(x)


class MergedSignalPeakCaller(PeakCaller):
    def __init__(self, bedgraphs, strand="+"):
        """
        Class to handle averaging of multiple bedgraphs

        Parameters
        ----------

        bedgraphs : list
            Filenames or BedTool objects of bedGraphs to be averaged

        """
        logger.info("Unioning bedgraphs...")
        df = pybedtools.BedTool().union_bedgraphs(i=bedgraphs).to_dataframe()
        logger.info("Averaging bedgraphs...")

        avg = df.iloc[:, 3:].mean(axis=1)
        df_merged = df.loc[:, ["chrom", "start", "end"]]
        df_merged["value"] = avg
        df_merged["value"].fillna(0)
        bedgraph = pybedtools.BedTool.from_dataframe(df_merged)
        super().__init__(bedgraph.fn, strand=strand)


class MultiPeakIDR(object):
    def __init__(self, peaks, oracle, strand="+"):
        """
        Class to handle running IDR with more than 2 replicates, which default
        IDR does not handle. Here we run all pairwise IDR, and then select the
        min number of peaks under the IDR threshold and then return that many
        from the provided oracle.

        Parameters
        ----------
        peaks : list
            List of narrowPeak files or pybedtools.BedTool objects pointing to
            narrowPeak files

        oracle : string or pybedtools.BedTool
            Peaks to pull from, generally from original peaks that have been
            merged in some way.

        strand : +, -, .
            Assumes the entire object represents a single strand; specify it
            here.
        """
        #: list of peaks
        self.peaks = peaks

        #: BedTool of merged peaks to uses as oracle
        self.oracle = pybedtools.BedTool(oracle)

        #: This object represents a single strand indicated here
        self.strand = strand

        # Simplified from idr.load_samples()
        self.signal_type = "signal.value"
        self.signal_index = 6
        self.peak_merge_fn = sum
        self.summit_index = 9

        #: Peaks loads as internal IDR data structures
        self.fps = [
            idr.load_bed(open(fn), self.signal_index, self.summit_index)
            for fn in self.peaks
        ]

        self.oracle_peaks = idr.load_bed(
            open(self.oracle.fn), self.signal_index, self.summit_index
        )

        # self._build_oracle()

        #: Holds information from running IDR.
        #: Keys are frozenset([i, j]) indicating the pairwise IDRs between
        #: peaks i and j.
        self.idrs = {}

    def _build_oracle(self):
        """
        Attempts as building an oracle. Deprecated, but retaining as fodder.
        """
        logger.info("Building oracle peaks...")

        # cat-and-merge strategy
        if 0:
            oracle = (
                pybedtools.BedTool.from_dataframe(
                    pybedtools.BedTool(self.peaks[0])
                    .cat(*self.peaks[1:], o="sum", c=5)
                    .to_dataframe()
                    .sort_values("name", ascending=False)
                )
                .each(to_narrowpeak)
                .saveas()
            )

        # multiintersect strategy
        if 0:
            h = pybedtools.BedTool().multi_intersect(i=self.peaks, cluster=True)

            lim = str(len(self.peaks))

            def filt(x):
                if x[3] != lim:
                    return
                return pybedtools.create_interval_from_list(
                    [x.chrom, str(x.start), str(x.stop)]
                )

            oracle = h.each(filt).saveas()

        # clustered strategy
        if 1:
            clustered = (
                pybedtools.BedTool(self.peaks[0])
                .cat(*self.peaks[1:], postmerge=False)
                .sort()
                .cluster()
                .to_dataframe()
            )

            def gen():
                for _, group in clustered.groupby("blockSizes"):
                    score = group["score"].sum()
                    start = group["start"].min()
                    stop = group["end"].max()
                    chrom = group["chrom"].unique()[0]
                    yield pybedtools.create_interval_from_list(
                        [
                            chrom,
                            str(start),
                            str(stop),
                            ".",
                            ".",
                            self.strand,
                            str(score),
                            "-1",
                            "-1",
                            "-1",
                        ]
                    )

            oracle = sort_by_score(pybedtools.BedTool(gen()).saveas())

        # IDR internal strategy
        if 0:
            oracle = self._multiway_merge()

        # By the time we get here, should have `oracle`
        self.oracle = oracle
        self.oracle_peaks = idr.load_bed(
            open(oracle.fn), self.signal_index, self.summit_index
        )

    def _multiway_merge(self):
        """
        Run IDR's internal routine for merging peaks.

        Uses self._multiway_merge_bed() to convert this to a BED file.
        """
        return idr.merge_peaks(
            self.fps,
            self.peak_merge_fn,
            self.oracle_peaks,
            use_nonoverlapping_peaks=False,
        )

    def _multiway_merge_bed(self):
        """
        Returns a BED6 of the multiway-merge object.
        """

        def gen0():
            for i, m_pk in enumerate(self._multiway_merge()):

                # from idr.build_idr_output_line_with_bed6
                yield pybedtools.create_interval_from_list(
                    [
                        m_pk.chrm,
                        str(m_pk.start),
                        str(m_pk.stop),
                        ".",
                        str(m_pk.merged_signal),
                        self.strand,
                    ]
                )

        return pybedtools.BedTool(gen0())

    def _build_merged(self, idx1, idx2):
        """
        Initial stage used by IDR.

        Uses IDR's internal routine for merging peaks. This is intended to be
        called by self.idr, which only works with 2 replicates at a time, hence
        the hard-coding of idx1 and idx2. See self._multiway_merge() for
        merging more than 2 replicates.

        Parameters
        ----------

        idx1, idx2 : int
            Indexes into self.peaks


        Returns
        -------
        idr
        """
        logger.info(f"Merging peaks for {self.peaks[idx1]} and {self.peaks[idx2]}")
        fn1 = self.peaks[idx1]
        fn2 = self.peaks[idx2]
        f1, f2 = [
            idr.load_bed(open(fp), self.signal_index, self.summit_index)
            for fp in [fn1, fn2]
        ]
        merged_peaks = idr.merge_peaks(
            [f1, f2],
            self.peak_merge_fn,
            self.oracle_peaks,
            use_nonoverlapping_peaks=False,
        )
        return merged_peaks

    def idr(self, idx1, idx2):
        """
        Run IDR between two sets of peaks

        Parameters
        ----------

        idx1, idx2 : int
            Indexes into self.peaks

        Returns
        -------
        None, but as a side effect this method populates the self.idrs
        dictionary for the key frozenset((idx1, idx2)). The value is another
        dictionary containing keys "IDRs", "localIDRs", and "merged_peaks". The
        values of these are the corresponding internal idr package data
        structures.
        """
        key = frozenset([idx1, idx2])
        if key in self.idrs:
            raise ValueError(f"key {key} exists")
        merged_peaks = self._build_merged(idx1, idx2)
        logger.info(f"Calcluating IDR for {self.peaks[idx1]} and {self.peaks[idx2]}")
        r1, r2 = idr.build_rank_vectors(merged_peaks)
        localIDRs = idr.fit_model_and_calc_local_idr(
            r1,
            r2,
            starting_point=(
                idr.idr.DEFAULT_MU,
                idr.idr.DEFAULT_SIGMA,
                idr.idr.DEFAULT_RHO,
                idr.idr.DEFAULT_MIX_PARAM,
            ),
            max_iter=idr.idr.MAX_ITER_DEFAULT,
            convergence_eps=idr.idr.CONVERGENCE_EPS_DEFAULT,
            fix_mu=False,
            fix_sigma=False,
        )
        IDRs = idr.calc_global_IDR(localIDRs)
        self.idrs[key] = dict(IDRs=IDRs, localIDRs=localIDRs, merged_peaks=merged_peaks)

    def _output(self, idx1, idx2):
        """
        Runs IDR's output routine

        Returns
        -------
        Generator of narrowPeak lines
        """
        key = frozenset([idx1, idx2])
        if key not in self.idrs:
            self.idr(idx1, idx2)
        d = self.idrs[key]
        IDRs = d["IDRs"]
        localIDRs = d["localIDRs"]
        merged_peaks = d["merged_peaks"]
        for localIDR, IDR, merged_peak in zip(localIDRs, IDRs, merged_peaks):
            line = idr.build_idr_output_line_with_bed6(
                merged_peak, IDR, localIDR, "narrowPeak", self.signal_type
            )
            yield line

    def npeaks_below_idr(self, thresh=0.05):
        """
        Dictionary of peak counts falling below IDR threshold.

        Returns
        -------
        Dictionary of the number of peaks falling below `thresh` in each
        pairwise IDR run.
        """
        counts = {}
        for i, j in itertools.combinations(range(len(self.peaks)), 2):
            c = 0
            for line in self._output(i, j):
                toks = line.split("\t")
                local_idr = float(toks[10])
                if local_idr >= -math.log10(thresh):
                    c += 1
            counts[(i, j)] = c
        return counts

    def final(self, thresh=0.05, use="oracle"):
        """
        Generate final peaks.

        Specificially, this extracts the top N peaks from the oracle peaks
        where N is determined by the minimum number of peaks below the IDR
        threshold across all pairwise IDR runs between replicates.

        Parameters
        ----------

        thresh : float
            IDR threshold

        use : oracle | idr-merged
            If "oracle", the final peaks will be selected from self.oracle.
            If "idr-merged", use IDR's internal merging routine, which allows
            multi-way merging if using their internal API.

        Returns
        -------
        BedTool of final peaks.
        """
        n = min(self.npeaks_below_idr(thresh).values())
        limit = n - 1
        if use == "oracle":

            def gen():
                for i, feature in enumerate(self.oracle):
                    if i >= limit:
                        break
                    yield feature

        elif use == "idr-merged":

            def gen():
                for i, m_pk in enumerate(self._multiway_merge()):
                    if i >= limit:
                        break

                    # from idr.build_idr_output_line_with_bed6
                    yield pybedtools.create_interval_from_list(
                        [m_pk.chrm, str(m_pk.start), str(m_pk.stop), self.strand]
                    )

        return pybedtools.BedTool(gen()).saveas()


def idr_peak_calls(
    bedgraphs,
    strand,
    thresh=0.05,
    oracle_fn="oracle.narrowPeak",
    final_fn="final.narrowPeak",
):
    """
    Returns oracle peaks and final peaks meeting IDR cutoff.

    Parameters
    ----------

    bedgraphs : list
        filenames of bedGraph files from replicates. Expected to be already
        normalized.

    strand : str
        One of "+", "-", or ".". Used to fill in the strand field of created
        narrowPeak files.

    thresh : float
        IDR threshold.

    oracle_fn : str
        Filename for "oracle" peaks. These are the peaks called after merging
        together all provided input bedGraphs.

    final_fn : str
        Filename for final thresholded peaks.
    """
    peak_callers = [PeakCaller(bedgraph, strand) for bedgraph in bedgraphs]
    mg = MergedSignalPeakCaller([i.bedgraph for i in peak_callers], strand=strand)
    oracle = mg.call_peaks().peaks_to_bed().moveto(oracle_fn)

    peaks = [pc.call_peaks().peaks_to_bed() for pc in peak_callers]
    m = MultiPeakIDR([p.fn for p in peaks], oracle=oracle, strand=strand)
    f = m.final(thresh=thresh).saveas(final_fn)


def to_narrowpeak(f):
    """
    Convert a feature into narrowPeak format, with signal and pval equivalent
    to the score.
    """
    return pybedtools.create_interval_from_list(
        [
            f.chrom,
            str(f.start),
            str(f.stop),
            ".",
            f.name,
            f.strand,
            f.name,
            f.name,
            "-1",
            str(int((f.stop - f.start) / 2)),
        ]
    )


def sort_by_score(x):
    """
    Sort a BedTool object by the score column.
    """
    df = pybedtools.BedTool(x).to_dataframe()
    df = df.sort_values("score", ascending=False)
    return pybedtools.BedTool.from_dataframe(df)


def merge_narrowbed(peaks, strand, additional_kwargs={"d": 1, "o": "max"}):
    """
        Method for merging narrowPeak files with bedtools.

        Using basic bedtools merge, merging narrowPeak files gets awkward if we
        want to output a valid narrowPeak. Here it's handled via pandas. Note
        that any narrowPeak summit positions are reset to -1 since it's not
        clear how to meaningfully aggregate them.

        Parameters
        ----------
        peaks : pybedtools.BedTool object
            Peaks to merge

        strand : str
            One of '+', '-', '.' to be set as the strand for each merged
            feature.

        additional_kwargs : dict
            Additional kwargs to send to pybedtools.BedTool.merge. By default,
            this merges features overlapping by 1 bp, and aggregates them by
            taking the maximum value. During testing, 'max' seemed to give
            better results than 'mean' because the latter tended to wash out
            strong peaks near smaller peaks.
        """
    x = (peaks.cut([0, 1, 2, 4]).merge(c=4, **additional_kwargs)).to_dataframe()
    x["score"] = "."
    x["strand"] = strand
    y = pybedtools.BedTool.from_dataframe(x)
    return y.each(to_narrowpeak).saveas()

