#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import pybedtools
from pybedtools import featurefuncs
import csv
import subprocess as sp
import gzip
import pysam




class CuratePeaks:
    def __init__(self, sample, narrowPeak, bw, cluster_length, output, region=None):
        """
        Return a 1bp-coordinate narrowPeak file corresponding to the highest score coordinate
        within cluster distance

        Parameters
        ----------

        sample : str
            Label for output

        narrowPeak : str
            Path to strand-specific narrowPeak file, where each interval represents the full
            size of the detected peak

        bw : list
            List of strand-specific bigWigs corresponding to the narrowPeak argument

        cluster_length : int
            If multiple peaks are found within this distance of each other,
            then select the highest peak among them and only report that one.

        output : str
            Path to output directory. Will be created if it does not already
            exist.

        region : str
            Subset region to use for this object. Useful for spot-checking and
            testing.
        """

        self.region = region
        self.output = output

        self.narrowpeak_df = pybedtools.BedTool(narrowPeak).to_dataframe()

        # incoming narrowPeak intervals may be large; select the single-bp
        # highest point within each interval
        self.onebp = self.precise_peak(sample, self.narrowpeak_df, bw)

        # there may be several 1-bp peaks in the above output that are within
        # `cluster_length`. Select only the highest within that range.
        self.curated = self.max_in_cluster(self.onebp, cluster_length)

    def precise_peak(self, sample, narrowpeak_df, bw):
        """
        Return a 1bp-coordinate dataframe corresponding to the highest score
        within narrowPeak intervals.

        If this object was provided a `region`, only that subsetted region will
        be used.
        """
        if not os.path.exists(self.output):
            os.makedirs(self.output)

        npzfn = os.path.join(self.output, sample + ".npz")
        rawfn = os.path.join(self.output, sample + ".raw.tsv")

        # Use multiBigwigSummary to average the bigwigs, and select the highest
        # single-bp position within the interval to be reported

        if self.region:
            region_arg = ['--region', self.region]
        else:
            region_arg = []
        sp.run(
            [
                'multiBigwigSummary', 'bins', '-bs', '1', '-b'
            ] + bw + region_arg + [
                '-o', npzfn, '--outRawCounts', rawfn
            ],
            check=True)

        raw = pd.read_csv(rawfn, sep="\t")
        # raw is a dataframe like:
        #
        #  chr start end  bw1 bw2 bw3
        #  chr 0     1    5   6   10
        #  chr 1     2    8   3   4

        raw["sum"] = raw.iloc[:, -len(bw) :].sum(axis=1)
        #  chr start end  bw1 bw2 bw3 sum
        #  chr 0     1    5   6   10  21
        #  chr 1     2    8   3   14  25

        raw["'end'"] = pd.to_numeric(raw["'end'"])

        def highest(interval, raw):
            subraw = raw[
                (raw["#'chr'"] == interval['chrom']) & \
                (raw["'end'"] >= interval['start']) & \
                (raw["'end'"] <= interval['end'])
            ]
            subraw = subraw.loc[subraw['sum'].idxmax(),:]
            return subraw[["'end'", 'sum']]

        # narrowpeak_df is a dataframe like:
        #
        # chrom  start  end  name  score  strand  ...
        # chr    1      3    .     1000   +       ...
        narrowpeak_df[['highest_coord', 'sum']] = narrowpeak_df.apply(
            lambda x: highest(x, raw), axis=1)
        # chrom  start  end  name  score  strand  ...  highest_coord  sum
        # chr    1      3    .     1000   +       ...  2              25

        narrowpeak_df['start'] = narrowpeak_df['highest_coord']
        narrowpeak_df['end'] = narrowpeak_df['highest_coord']
        # chrom  start  end  name  score  strand  ...  highest_coord  sum
        # chr    2      2    .     1000   +       ...  2              25

        # Replace the name '.' by a combination of chr_start_end
        narrowpeak_df.iloc[:, 3] = narrowpeak_df.iloc[:, 0:3].apply(
            lambda x: "_".join(x.map(str)), axis=1
        )
        # Sort by chromosome then by start coordinate
        narrowpeak_df = narrowpeak_df.sort_values(["chrom", "start"])
        return narrowpeak_df

    def max_in_cluster(self, onebp, cluster_length):
        """
        Return the peak with highest score among peaks within cluster range
        """
        df = (
            pybedtools.BedTool.from_dataframe(onebp, header=None)
            .cluster(d=cluster_length)
            .to_dataframe(header=None)
        )
        df.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                      'na1', 'na2', 'na3', 'na4', 'highest_coord', 'sum', 'cluster']
        idx = df.groupby('cluster')["score"].transform(max) == df["score"]
        return df[idx].iloc[:, 0:10]


class PrepGtfs:
    def __init__(self, gtf):
        """
        Splits the annotation file in 2 dataframes of ORFs minus their 3'end coordinate,
        and another of the 3'end coordinates only.

        Parameters
        ----------

        gtf : str
            Path to annotation gtf file, assumes it contains the mRNAs, sRNAs, tRNAs and rRNAs

        intergenic_min_size: integer
            Minimun intergenic region size under which the 2 same-strand adjacent ORFs
            will be considered too close to assign the peak to a promoter

        """
        self.full = pybedtools.BedTool(gtf)
        self.orf = self.split_3end(gtf)[0]
        self.orf3end = self.split_3end(gtf)[1]

    def split_3end(self, gtf):
        """
        Return 2 dataframes corresponding to the ORFs minus the 3' end coordinate (gtfintern)
        and to the 3'end coordinate only (gtf3prime).
        """
        gtf = pybedtools.BedTool(gtf)
        gtf3prime = gtf.each(
            featurefuncs.three_prime,
            downstream=0,
            upstream=1).saveas()
        gtfintern = gtf.subtract(gtf3prime, s=True)
        return [gtfintern.to_dataframe(), gtf3prime.to_dataframe()]



class ClassAssign:
    def __init__(
        self,
        peakdf,
        gtfs,
        fasta,
        param_antisense,
        param_down,
        trRNA,
        param_upstart,
        param_downstart,
        param_upstop,
        param_fasta,
        param_downfasta,
        output,
        sample,
        offset=None
    ):
        """
        Assign peaks to particular classes:
        - primary: within 3'end of any ORF (mRNA, tRNA, rRNA, sRNA) included and param_down-bp
          downstream on the same strand AND has the highest readcount of all such peaks
          associated within the same region.
        - secondary: fulfills the above criteria with respect to location BUT
          is NOT the peak with the highest readcount.
        - antisense: located within param_antisense-bp upstream, downstream or in an
          ORF of the opposite strand.
        - internal: within an any ORF (mRNA, tRNA, rRNA, sRNA) coordinates, excluding the
          3'end coordinate on the same strand.
        - orphan: not associated with any of the above categories.
        Peaks can have multiple classifications.
        Also returns lists of peaks within param_upstart-bp upstream of start codon to
        param_downstart-bp downstream of start codon, and within param_upstop-bp downstream
        of start codon to the stop codon.

        Parameters
        ----------

        peakdf: dataframe
            List of curated 1bp-peaks. Typically this is the output `curated` of class CuratePeaks

        gtfs: list of 2 dataframes
            A dataframe of ORFs minus their 3'end coordinate, and a dataframe of the 3'end
            coordinates only. Typically, this is the output of class PrepGtfs

        fasta: str
            Filename for genome fasta

        param_antisense: int
            Threshold distance to ORF to assign to class antisense

        param_down: int
            Threshold distance downstream of 3'end of ORF to assign to class
            Primary or Secondary
        trRNA: str
            Filename to file containing the list of tRNAs and rRNAs. Names can be exact or
            regex.
        param_upstart: int
            Threshold distance upstream of 5'end of an mRNA ORF to include when
            returning the list of peaks in promoters

        param_downstart: int
            Threshold distance downstream of 5'end of an mRNA ORF to include when
            returning the list of peaks in promoters

        param_upstop: int
            Threshold distance downstream of 5' end of an mRNA ORF to exclude when
            return the list of peaks upstream of 3' of mRNA ORF

        param_fasta: int
            Length upstream of the 3'end of ORF for which the fasta sequence will be returned

        param_downfasta: int
            Length downstream of the 3'end of ORF for which the fasta sequence will be returned

        output: str
            Path to output directory

        sample: str
            Name of sample
        """
        self.peakdf = peakdf
        self.fai = self.index_func(fasta)
        self.prepeaks = pybedtools.BedTool.from_dataframe(
            peakdf[["chrom", "start", "end", "name", "score", "strand"]]
        )
        self.prebedmRNA = pybedtools.BedTool.from_dataframe(gtfs.orf)
        self.prebedmRNA3 = pybedtools.BedTool.from_dataframe(gtfs.orf3end)

        # if the offset is set, shift all coordinates
        self.offset = offset
        self.bedpeaks = self.prepeaks.shift(g=self.fai[0], s=-self.offset) if \
            self.offset else self.prepeaks
        self.bedmRNA = self.prebedmRNA.shift(g=self.fai[0], s=-self.offset) if \
            self.offset else self.prebedmRNA
        self.bedmRNA3 = self.prebedmRNA3.shift(g=self.fai[0], s=-self.offset) if \
            self.offset else self.prebedmRNA3
        self.gtfsfull = gtfs.full.shift(g=self.fai[0], s=-self.offset) if \
            self.offset else gtfs.full

        self.antisense = self.antisense_func(
            self.peakdf, self.bedpeaks, self.bedmRNA, param_antisense, self.fai[0]
        )
        self.internal = self.internal_func(
            self.peakdf, self.bedpeaks, self.gtfsfull, self.fai[0]
        )
        self.primary = self.primary_func(
            self.peakdf, self.bedpeaks, self.bedmRNA3, param_down, trRNA, self.fai, output, sample
        )
        self.details_class = self.details_class_func(
            self.primary, self.antisense, self.internal
        )
        self.upstart_upstop = self.upstart_upstop_func(
            self.peakdf,
            self.bedpeaks,
            trRNA,
            gtfs.orf,
            self.bedmRNA3,
            param_upstart,
            param_downstart,
            param_upstop,
        )
        self.fasta = self.fasta_func(
            self.peakdf, self.fai, param_fasta, param_downfasta
        )

    def index_func(self, fasta):
        """ Return index file after creating it if did not exist"""
        if not os.path.exists(fasta + '.fai'):
            # checks if fasta is gzipped
            with gzip.open(fasta, 'r') as fh:
                try:
                    fh.read(1)
                    # if gzipped, unzip
                    os.system('gunzip -c ' + fasta + ' > ' + fasta.replace('.gz', ''))
                    fasta = fasta.replace('.gz', '')
                finally:
                    pysam.faidx(fasta)
        return [fasta + '.fai', fasta]

    intsct_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                        'intsct_chrom', 'insct_source', 'intsct_feature',
                        'intsct_start', 'intsct_end', 'intsct_score',
                        'intsct_strand', 'intsct_na', 'intsct_name']

    def antisense_func(self,
                       peakdf,
                       bedpeaks,
                       bedmRNA,
                       param_antisense,
                       fai,
                       intsct_cols=intsct_cols):
        """ Determine peaks antisense from ORFs"""
        temp = bedpeaks.intersect(
                [
                    bedmRNA.slop(g=fai, b=param_antisense),
                ],
                loj=True,
                S=True,
            ).to_dataframe(header=None)
        temp.columns = intsct_cols
        # temp is a dataframe like:
        # chrom    start   end   strand   name           ...  instc_name
        # chr      227     230   +        chr_227_230_+  ...  gene_id "ID1";
        # chr      227     230   +        chr_227_230_+  ...  gene_id "ID2";

        peakdf["antisense"] = temp \
            .groupby("name")["intsct_name"] \
            .apply(lambda x: ", ".join(x)) \
            .to_frame()["intsct_name"]
        # peakdf:
        # name          chrom    start   end   strand   ...  antisense
        # chr_227_230_+ chr      227     230   +        ...  gene_id "ID1"; gene_id "ID2";

        peakdf["antisense"] = np.where(
            peakdf["antisense"] == ".", np.nan, ("A:" + peakdf["antisense"])
        )
        # peakdf:
        # name          chrom    start   end   strand   ...  antisense
        # chr_227_230_+ chr      227     230   +        ...  A:gene_id "ID1"; gene_id "ID2";

        return peakdf["antisense"]

    def internal_func(self, peakdf, bedpeaks, gtffull, fai, intsct_cols=intsct_cols):
        """ Determine peaks internal to ORF, in same orientation"""
        # need to shift internal coord of minus strand ORF 3'end by -1 because of bedtools intersect nt-base
        # use bedtools shift
        shifted = gtffull.shift(g=fai, p=0,m=1)
        temp = bedpeaks.intersect(
            shifted,
            loj=True, s=True, f=1
        ).to_dataframe(header=None)
        temp.columns = intsct_cols
        # temp is a dataframe like:
        # chrom    start   end   strand   name           ...  instc_name
        # chr      227     230   +        chr_227_230_+  ...  gene_id "ID3";
        # chr      227     230   +        chr_227_230_+  ...  gene_id "ID4";

        peakdf["internal"] = (
            temp
            .groupby("name")["intsct_name"]
            .apply(lambda x: ", ".join(x))
            .to_frame()["intsct_name"]
        )
        # peakdf:
        # name          chrom    start   end   strand   ...  internal
        # chr_227_230_+ chr      227     230   +        ...  gene_id "ID3"; gene_id "ID4";

        peakdf["internal"] = np.where(
            peakdf["internal"] == ".", np.nan, ("I:" + peakdf["internal"])
        )
        # peakdf:
        # name          chrom    start   end   strand   ...  internal
        # chr_227_230_+ chr      227     230   +        ...  I:gene_id "ID3"; gene_id "ID4";

        return peakdf["internal"]

    def primary_func(self,
                     peakdf,
                     bedpeaks,
                     bedmRNA3,
                     param_down,
                     trRNA,
                     fai,
                     output,
                     sample,
                     intsct_cols=intsct_cols):
        """
        Determine primary and secondary peaks:
        - Primary: within 3'end of any ORF (mRNA, tRNA, rRNA, sRNA) included
        and param_down-bp downstream on the same strand AND has the highest readcount of all such peaks
        associated within the same region.
        - Secondary: identical to Primary but not the highest readcount.
        """

        shifted = bedmRNA3.shift(g=fai[0], p=0,m=1)

        peaksmrna = bedpeaks.intersect(
                shifted.each(
                    featurefuncs.three_prime,
                    downstream=param_down,
                    upstream=0),
                loj=True,
                s=True,
            ).to_dataframe(header=None)
        peaksmrna.columns = intsct_cols
        # peaksmrna is a dataframe like:
        # chrom    start   end   strand   name           ...  instc_name
        # chr      227     230   +        chr_227_230_+  ...  gene_id "ID5";
        # chr      227     230   +        chr_227_230_+  ...  gene_id "ID6";

        peaksmrna["origin"] = np.where(peaksmrna["intsct_name"] != ".", "mRNA", np.nan)
        # chrom    start   end   strand   name           ...  instc_name        origin
        # chr      227     230   +        chr_227_230_+  ...  gene_id "ID5";    mRNA
        # chr      227     230   +        chr_227_230_+  ...  gene_id "tRNA1";  mRNA

        # remove rRNA and tRNA from mRNA origin
        for i in pd.read_csv(trRNA, sep="\t", header=None)[0]:
            peaksmrna["origin"] = np.where(
                peaksmrna["intsct_name"].str.contains(i), np.nan, peaksmrna["origin"]
            )
        # chrom    start   end   strand   name           ...  instc_name        origin
        # chr      227     230   +        chr_227_230_+  ...  gene_id "ID5";    mRNA
        # chr      227     230   +        chr_227_230_+  ...  gene_id "tRNA1";  NaN

        peaksmrna["origin"] = peaksmrna["origin"].replace("nan", np.NaN)

        # the primary peak is defined as the highest score peak within primary range
        # of a mRNA, tRNA, rRNA or sRNA
        idx = (
            peaksmrna.groupby("intsct_name")["score"].transform(max) == peaksmrna["score"]
        )
        peaksmrna["primary"] = np.where(
            (idx) & (peaksmrna["intsct_name"] != "."), peaksmrna["intsct_name"], np.nan
        )
        # chrom    start   end   strand   name           ...  instc_name       origin  primary
        # chr      227     230   +        chr_227_230_+  ...  gene_id "ID5";   mRNA    gene_id "ID5";
        # chr      227     230   +        chr_227_230_+  ...  gene_id "tRNA1"; NaN

        # calculate distances between gene 3'end and termination peak
        peaksmrna["distance"] = np.where(
            peaksmrna["strand"] == "+",
            peaksmrna["start"] - peaksmrna["intsct_start"] +1,
            peaksmrna["intsct_end"] - peaksmrna["start"] +1,
        )
        peaksmrna["distance"] = np.where(
            peaksmrna["intsct_name"] != ".", peaksmrna["distance"], np.nan
        )
        # chrom    start   end   ...  instc_name       origin  primary          distance
        # chr      227     230   ...  gene_id "ID5";   mRNA    gene_id "ID5";   5
        # chr      227     230   ...  gene_id "tRNA1"; NaN     gene_id "tRNA1"; 12

        # distance to mRNAs only (exclude tRNA and rRNA)
        peaksmrna["distance_to_mRNA"] = np.where(
            peaksmrna["origin"] == "mRNA", peaksmrna["distance"], np.nan
        )

        # save the distance matrix
        fn = os.path.join(output, sample + "_primary_distances.tsv")
        peaksmrna[
            [
                "chrom",
                "start",
                "end",
                "name",
                "intsct_name",
                "distance",
                "distance_to_mRNA",
            ]
        ].to_csv(fn, sep="\t", index=None, header=None, quoting=csv.QUOTE_NONE)

        # copy peaksmrna results over to peakdf primary
        peakdf["primary"] = (
            peaksmrna.groupby("name")["primary"]
            .apply(lambda x: ",".join(set(x.dropna())))
            .to_frame()["primary"]
        )
        # peakdf is a dataframe like:
        #               chrom   start   end  ...    primary
        # chr_228_228   chr     228     228  ...    gene_id "ID5";

        peakdf["primary"] = np.where(
            peakdf["primary"] == "", np.nan, ("P:" + peakdf["primary"])
        )
        #               chrom   start   end  ...    primary
        # chr_228_228   chr     228     228  ...    P: gene_id "ID5";

        peakdf["primary_origin"] = (
            peaksmrna.groupby("name")["origin"]
            .apply(lambda x: ",".join(set(x.dropna())))
            .to_frame()["origin"]
        )
        #               chrom   start   end  ...    primary            primary_origin
        # chr_228_228   chr     228     228  ...    P: gene_id "ID5";  mRNA

        peakdf["primary_distance_to_mRNA"] = (
            peaksmrna.groupby("name")["distance_to_mRNA"]
            .apply(lambda x: ",".join(set(x.dropna().astype(int).astype(str))))
            .to_frame()["distance_to_mRNA"]
        )
        #               chrom   start   end  ...    primary_origin  primary_distance_to_mRNA
        # chr_228_228   chr     228     228  ...    mRNA            5

        # secondary peaks
        peaksmrna["secondary"] = np.where(
            (~idx) & (peaksmrna["intsct_name"] != "."), peaksmrna["intsct_name"], np.nan
        )
        # peaksmrna is a dataframe like:
        # chrom start   end  ... primary        distance  distance_to_mRNA  secondary
        # chr   227     230 ...  gene_id "ID5"; 5         5
        # chr   227     230  ...                12                          gene_id "tRNA1";

        peakdf["secondary"] = (
            peaksmrna.groupby("name")["secondary"]
            .apply(lambda x: ",".join(x.dropna()))
            .to_frame()["secondary"]
        )
        #               chrom   start   ...    primary_origin  primary_distance_to_mRNA    secondary
        # chr_228_228   chr     228     ...    mRNA            5                           gene_id "tRNA1";

        peakdf["secondary"] = np.where(
            peakdf["secondary"] == "", np.nan, ("S:" + peakdf["secondary"])
        )
        #               chrom   start   ...    primary_origin  primary_distance_to_mRNA    secondary
        # chr_228_228   chr     228     ...    mRNA            5                           S: gene_id "tRNA1";

        return peakdf[
            ["primary", "primary_origin", "primary_distance_to_mRNA", "secondary"]
        ]

    def details_class_func(self, primary, antisense, internal):
        """ Concatenate the results from primary, antisense and internal """
        df = pd.concat([primary[["primary", "secondary"]], internal, antisense], axis=1)
        df["details"] = df.apply(lambda x: ",".join(x.dropna().map(str)), axis=1)
        # df is a dataframe like:
        #               primary             secondary             antisense
        # chr_228_228   P: gene_id "ID5";   S: gene_id "tRNA1";   A:gene_id "ID1"; gene_id "ID2";
        #
        #               internal                            details
        #               I:gene_id "ID3"; gene_id "ID4";     P: gene_id "ID5";S: gene_id "tRNA1";A:gene_id "ID1"; gene_id "ID2";I:gene_id "ID3"; gene_id "ID4";

        # add 'primary', 'secondary', ... to the classification
        df["classification"] = np.where(~df["primary"].isna(), "primary, ", "")
        df["classification"] = np.where(
            ~df["secondary"].isna(),
            (df["classification"] + "secondary, "),
            df["classification"],
        )
        df["classification"] = np.where(
            ~df["internal"].isna(),
            (df["classification"] + "internal, "),
            df["classification"],
        )
        df["classification"] = np.where(
            ~df["antisense"].isna(),
            (df["classification"] + "antisense, "),
            df["classification"],
        )
        df["classification"] = np.where(
            df["classification"] == "", "orphan", df["classification"]
        )
        #               primary             ...     details                                 classification
        # chr_228_228   P: gene_id "ID5";   ...     P: gene_id "ID5"; ... gene_id "ID4";    primary, sedcondary, antisense, internal
        # chr_450_450                                                                       orphan

        return df[["details", "classification"]]

    def upstart_upstop_func(
        self,
        peakdf,
        bedpeaks,
        trRNA,
        orf,
        bedmRNA3,
        param_upstart,
        param_downstart,
        param_upstop,
    ):
        # make bed mRNA without tRNA, rRNA (nor sRNA)
        mRNAonly = orf.copy()
        for i in pd.read_csv(trRNA, sep="\t", header=None)[0]:
            mRNAonly = mRNAonly[~mRNAonly["attributes"].str.contains(i)]
        premRNAonly = pybedtools.BedTool.from_dataframe(mRNAonly)
        bedmRNAonly = premRNAonly.shift(g=self.fai[0], s=-self.offset) if \
            self.offset else premRNAonly

        # label peaks within 10bp downstream of any ORF
        peakstmp = bedpeaks.intersect(
            [
                bedmRNA3.each(featurefuncs.three_prime, downstream=10, upstream=0),
            ],
            loj=True,
            s=True,
        ).to_dataframe(header=None)
        peakstmp.columns = [
            "chrom", "start", "end",
            "name", "score", "strand",
            "chr-intersect", "source-intersect", "feature-intersect",
            "start-intersect", "end-intersect", "score-intersect",
            "strand-intersect", "na-intersect", "name-intersect"
        ]
        # peakstmp is a dataframe like:
        # chrom     start   end     ...     name-intersect
        # chr       228     228     ...     gene_id "ID10"
        # chr       450     450     ...     .

        # discard peaks within 10bp downstream of an ORF
        peakstmp = peakstmp[peakstmp["chr-intersect"] != "."]
        # chrom     start   end     ... name-intersect
        # chr       228     228     ... gene_id "ID10"

        peakdf["within_10bp-downstream"] = (
            peakstmp.groupby("name")["name-intersect"]
            .apply(lambda x: ", ".join(x))
            .to_frame()["name-intersect"]
        )
        # peakdf is a dataframe like:
        #                   chrom     start   end     ...     within_10bp_dowsntream
        #  chr_228_228      chr       228     228     ...     gene_id "ID10"
        #  chr_450_450      chr       450     450     ...     NaN

        # calculate distances peaks to start of mRNA only ORF
        peakstmp = pd.concat(
            [
                bedpeaks.sort()
                .closest(
                    bedmRNAonly.each(
                        featurefuncs.five_prime, downstream=1, upstream=0
                    ).sort(),
                    D="b",
                    k=10,
                    s=True,
                )
                .to_dataframe(header=None)
            ]
        )
        peakstmp.columns = [
            "chrom", "start", "end",
            "name", "score", "strand",
            "chr-intersect", "source-intersect", "feature-intersect",
            "start-intersect", "end-intersect", "score-intersect",
            "strand-intersect", "na-intersect", "name-intersect", "distance"
        ]
        # peakstmp is a dataframe like:
        # chrom     start   end     ...     name-intersect  distance
        # chr       228     228     ...     gene_id "ID10"  19
        # chr       630     630     ...     gene_id "ID11"  299

        # add 1 to distances because closest calculate the length of interval between features
        peakstmp["distance"] = np.where(
            peakstmp["distance"] >= 0,
            peakstmp["distance"] + 1,
            peakstmp["distance"] - 1,
        )
        # start 200bp: subset peakstmp to peaks within param_upstart upstream and param_downstart
        # downstream of start of ORF
        peakstmp = peakstmp[
            (
                (
                    (peakstmp["distance"].abs() <= param_upstart)
                    & (peakstmp["distance"] < 0)
                )
                | (
                    (peakstmp["distance"].abs() <= param_downstart)
                    & (peakstmp["distance"].abs() >= 0)
                )
            )
        ]
        # chrom     start   end     ...     name-intersect  distance
        # chr       228     228     ...     gene_id "ID10"  19

        # copy results over to peakdf
        peakdf["upstart200"] = (
            peakstmp.groupby("name")["name-intersect"]
            .apply(lambda x: ", ".join(x))
            .to_frame()["name-intersect"]
        )
        peakdf["dist_upstart200"] = (
            peakstmp.groupby("name")["distance"]
            .apply(lambda x: ", ".join(map(str, x)))
            .to_frame()["distance"]
        )
        # peakdf is a dataframe like:
        #                   chrom     start   end     ...     upstart200        dist_upstart200
        #  chr_228_228      chr       228     228     ...     gene_id "ID10"    19
        #  chr_450_450      chr       450     450     ...     NaN               NaN

        # stop: within mRNA only ORF but not in the first 50bp
        mRNAonly51 = mRNAonly.copy()
        # bed of param_upstop + 1bp downstream of ORF start to ORF end
        mRNAonly51["start"] = np.where(
            mRNAonly51["strand"] == "+",
            mRNAonly51["start"] + param_upstop,
            mRNAonly51["start"],
        )
        mRNAonly51["end"] = np.where(
            mRNAonly51["strand"] == "+",
            mRNAonly51["end"],
            mRNAonly51["end"] - param_upstop,
        )
        bedmRNAonly51 = pybedtools.BedTool.from_dataframe(mRNAonly)

        #  peaks intersecting with mRNAonly51
        bedtmp = bedpeaks.sort().intersect(bedmRNAonly51, s=True, loj=True)
        peakstmp = bedtmp.closest(
            bedmRNAonly.each(  # distance to start codon
                featurefuncs.five_prime, downstream=1, upstream=0
            ).sort(),
            D="b",
            k=10,
            s=True,
        ).to_dataframe(header=None)

        peakstmp.columns = [
            "chrom", "start", "end",
            "name", "score", "strand",
            "chr-intersect", "source-intersect", "feature-intersect",
            "start-intersect", "end-intersect", "score-intersect",
            "strand-intersect", "na-intersect", "name-intersect",
            "chr-closest", "source-closest", "feature-closest",
            "start-closest", "end-closest", "score-closest",
            "strand-closest", "na-closest", "name-closest",
            "distance"
        ]
        # chrom     start ... name-intersect    ...     name-closest    distance
        # chr       93        gene_id "ID30";   ...     gene_id "ID30"; 63
        # chr       842       -1                ...     gene_id "ID40"; 183

        # only keep distances for same gene than intersect
        peakstmp = peakstmp[
            (
                (peakstmp["name-intersect"] == peakstmp["name-closest"])
                & (peakstmp["distance"] >= param_upstop)
            )
        ]
        # chrom     start ... name-intersect    ...     name-closest    distance
        # chr       93        gene_id "ID30";   ...     gene_id "ID30"; 63

        # fix distances to ORF, because closest give to outside of ORF
        peakstmp["distance"] = np.where(
            peakstmp["distance"] >= 0,
            peakstmp["distance"] + 1,
            peakstmp["distance"] - 1,
        )
        # copy over results to peakdf
        peakdf["upstop"] = (
            peakstmp.groupby("name")["name-closest"]
            .apply(lambda x: ", ".join(x))
            .to_frame()["name-closest"]
        )
        peakdf["dist_upstop"] = (
            peakstmp.groupby("name")["distance"]
            .apply(lambda x: ", ".join(map(str, x)))
            .to_frame()["distance"]
        )
        return peakdf[
            [
                "within_10bp-downstream",
                "upstart200",
                "dist_upstart200",
                "upstop",
                "dist_upstop",
            ]
        ]

    def fasta_func(self, peakdf, fai, param_fasta, param_downfasta):
        # add fasta sequences
        dfplus = pybedtools.BedTool.from_dataframe(
            peakdf[peakdf["strand"] == "+"], header=None
        )
        dfminus = pybedtools.BedTool.from_dataframe(
            peakdf[peakdf["strand"] == "-"], header=None
        )
        if self.offset:
            dfplus = dfplus.shift(g=fai[0], s=-self.offset)
            dfminus = dfminus.shift(g=fai[0], s=-self.offset)
        genome = {}
        with open(fai[0]) as f:
            for line in f:
                chrom = line.split()
                genome[chrom[0]] = (0, int(chrom[1]))
        # want params.fasta bp upstream of the 3'end, 3'end then 10 bp downstream, so need to add 1 upstream
        dfplus2 = dfplus.each(
            featurefuncs.three_prime,
            downstream=param_downfasta,
            upstream=int(param_fasta) + 1,
            genome=genome,
        ).saveas()
        # same for minus strand, but was shifted by 1bp when on the minus strand
        dfminus2 = dfminus.each(
            featurefuncs.three_prime,
            downstream=int(param_downfasta) + 1,
            upstream=param_fasta,
            genome=genome,
        ).saveas()
        # concatenate both strands
        df2 = dfplus2.cat(dfminus2, postmerge=False).sort()
        # truncate to chromosome size in case an interval would go over
        df2 = df2.truncate_to_chrom(genome)
        # get sequence
        df3 = df2.sequence(fi=fai[1], s=True)
        fatmp = pd.read_csv(df3.seqfn, header=None).iloc[1::2, :]
        fatmp["idx"] = list(peakdf.index)
        fatmp = fatmp.set_index("idx")
        fatmp.columns = ["upstream_3end_seq"]
        return fatmp


def assign(
    sample,
    narrowPeak,
    bw,
    fasta,
    gtf,
    trRNA,
    output='data',
    cluster_length=100,
    param_antisense=50,
    param_down=50,
    param_upstart=200,
    param_downstart=50,
    param_upstop=50,
    param_fasta=50,
    param_downfasta=10,
    region=None,
    offset=None
):
    """
    Parameters
    ----------

    sample : str
        Label used to compose output files

    narrowPeak : str
        Path to narrowPeak output file created by termseq-peaks.

    bw : list
        List of paths to bigWig files representing the 3' ends of reads. See
        docs for details; often this is from unique reads taking the first bp
        of R1 using deeptools. These should be normalized already (see docs)

    fasta : str
        Path to reference genome fasta, used to extract sequence around peaks

    gtf : str
        Path to gtf annotations.

    trRNA: str
        Path to file containing the names of tRNAs and rRNAs. Can be exact gene
        names or regex.

    output: str
        Path to output directory. Default is `data`.

    cluster_length: int
        If multiple peaks are found within this distance of each other,
        then select the highest peak among them and only report that one. Default is 100.

    param_antisense: int
        Threshold distance to ORF to assign to class antisense. Default is 50.

    param_down: int
        Threshold distance downstream of 3'end of ORF to assign to class
        Primary or Secondary. Default is 50.

    param_upstart: int
        Threshold distance upstream of 5'end of an mRNA ORF to include when
        returning the list of peaks in promoters. Default is 200.

    param_downstart: int
        Threshold distance downstream of 5'end of an mRNA ORF to include when
        returning the list of peaks in promoters. Default is 50.

    param_upstop: int
        Threshold distance downstream of 5' end of an mRNA ORF to exclude when
        return the list of peaks upstream of 3' of mRNA ORF. Default is 50.

    param_fasta: int
        Length upstream of the 3'end of ORF for which the fasta sequence will be returned.
        Default is 50.

    param_downfasta: int
        Length downstream of the 3'end of ORF for which the fasta sequence will be returned.
        Default is 10.
    """
    peaks = CuratePeaks(sample, narrowPeak, bw, cluster_length, output, region)
    gtfs = PrepGtfs(gtf)
    peak = peaks.curated.set_index("name", drop=False)

    peakdf = peak.copy()
    # The class ClassAssign returns each of the assignment columns as class properties
    assigned = ClassAssign(
        peakdf=peakdf,
        gtfs=gtfs,
        fasta=fasta,
        param_antisense=param_antisense,
        param_down=param_down,
        trRNA=trRNA,
        param_upstart=param_upstart,
        param_downstart=param_downstart,
        param_upstop=param_upstop,
        param_fasta=param_fasta,
        param_downfasta=param_downfasta,
        output=output,
        sample=sample,
        offset=offset
    )

    # Concatenate all the AssignClass properties into one full dataframe
    full = pd.concat(
        [
            peak,
            assigned.antisense,
            assigned.internal,
            assigned.primary,
            assigned.details_class,
            assigned.upstart_upstop,
            assigned.fasta,
        ],
        axis=1,
    )
    # save the .narrowPeak files
    narrowPeakCols = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "na1",
        "na2",
        "na3",
        "na4",
    ]
    # narrowPeak of full dataframe
    full[narrowPeakCols].sort_values(by="score", ascending=False).to_csv(
        os.path.join(output, "all." + sample + ".narrowPeak"),
        sep="\t",
        index=None,
        header=None,
        quoting=csv.QUOTE_NONE,
    )
    # narrowPeak of peaks in the promoters but not within 10bp dowsntream of start
    full[full["upstart200"].notnull() & full["within_10bp-downstream"].isnull()][
        narrowPeakCols
    ].sort_values(by="score", ascending=False).to_csv(
        os.path.join(
            output, "all.200bp-start.minus10bpdownstream." + sample + ".narrowPeak"
        ),
        sep="\t",
        index=None,
        header=None,
        quoting=csv.QUOTE_NONE,
    )
    # narrowPeak of peaks in the promoters
    full[full["upstart200"].notnull()][narrowPeakCols].sort_values(
        by="score", ascending=False
    ).to_csv(
        os.path.join(output, "all.200bp-start." + sample + ".narrowPeak"),
        sep="\t",
        index=None,
        header=None,
        quoting=csv.QUOTE_NONE,
    )
    # narrowPeak of peaks within ORF but not within param_upstop of start
    full[full["upstop"].notnull()][narrowPeakCols].sort_values(
        by="score", ascending=False
    ).to_csv(
        os.path.join(output, "all.51bp-stop." + sample + ".narrowPeak"),
        sep="\t",
        index=None,
        header=None,
        quoting=csv.QUOTE_NONE,
    )
    # save .tsv of the full dataframe
    full[
        [
            "chrom",
            "start",
            "strand",
            "classification",
            "details",
            "primary_origin",
            "score",
            "upstream_3end_seq",
            "upstart200",
            "dist_upstart200",
            "within_10bp-downstream",
            "upstop",
            "dist_upstop",
        ]
    ].sort_values(by="score", ascending=False).to_csv(
        os.path.join(output, "all." + sample + ".tsv"),
        sep="\t",
        index=None,
        header=True,
        quoting=csv.QUOTE_NONE,
    )


def table_output(sample, assigned, assigned2, kinefold_scores, output='data'):
    """
    Parameters
    ----------

    assigned : str
        Path to file with curated peaks. Typically this is the output of the assign function.
        Can be for a single strand or both.

    assigned2 : str
        Optional. Path to file with curated peaks of opposite strand if parameter assigned
        was for a single strand. Typically this is the output of the assign function.

    kinefold_scores : str
        Optional. Path to Kinefold output file corresponding to the curated peaks.

    output: str
        Path to output directory. Default is `data`.
    """

    score_cols = [
        "3' end",
        "Strand",
        "3' end_seq",
        "Kinefold_structure",
        "A-tract",
        "Hairpin",
        "Hp_structure",
        "Loop_seq",
        "U-Tract",
        "dGU",
        "dGL",
        "dGH",
        "dGHA",
        "dGA",
        "dGB",
        "TS",
        "NA",
    ]
    if assigned2:
        df = pd.concat(
            [pd.read_csv(assigned, sep="\t"), pd.read_csv(assigned2, sep="\t")]
        )
    else:
        df = pd.read_csv(assigned, sep="\t")
    df = df.reset_index()
    if kinefold_scores:
        scores = pd.read_csv(kinefold_scores, sep="\t", names=score_cols)
        scores = scores.iloc[1:, :].reset_index()
        tableS1 = df[
            [
                "chrom",
                "start",
                "strand",
                "score",
                "classification",
                "details",
                "upstream_3end_seq",
            ]
        ].join(
            scores[
                [
                    "3' end",
                    "TS",
                    "Kinefold_structure",
                    "A-tract",
                    "Hairpin",
                    "Hp_structure",
                    "Loop_seq",
                    "U-Tract",
                ]
            ]
        )
    else:
        tableS1 = df[
            [
                "chrom",
                "start",
                "strand",
                "score",
                "classification",
                "details",
                "upstream_3end_seq",
            ]
        ]

    tableS1.sort_values(["chrom", "start", "strand"]).to_csv(
        os.path.join(output, sample + "_TableS1.tsv"),
        sep="\t",
        header=True,
        index=False,
        quoting=csv.QUOTE_NONE,
    )
