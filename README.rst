Term-seq peak-caller
====================

This tool was designed to call 3'-end peaks from term-seq data in bacteria,
handling replicates in a statistically robust manner.

Installation
------------

pip::

   pip install termseq-peaks

conda::

   conda install termseq-peaks --channel conda-forge --channel bioconda

source::

   git clone <repo>
   cd <repo>
   python setup.py install

Background
----------

For peak-calling, the venerable macs2 peak-caller [0] is the go-to method.
However, macs2 is designed for ChIP-seq data and all the assumptions that go
along with that (modeling peaks based on fragment size, controlling for open
chromatin, comparing IP against a background input, etc). We found that naively
applying macs2 to term-seq data resulted in suboptimal peak calls. Term-seq
signal is composed of single-bp positions of read ends and has a very high
dynamic range since it is coming from trancribed RNA that can have very high
copy numbers per cell (in contrast to ChIP-seq which we have just 1 copy per
haploid cell).

Furthermore, many peak-callers have limited support for multiple replicates.
One general solution to this is to apply the Irreproducible Discovery Rate
method (IDR), originally developed for the ENCODE project. By design, the IDR
method only takes two replicates at a time.

This tool takes two novel approaches that together yield empirically better
results for calling precise peaks in term-seq data: 1) a signal processing
approach and 2) an implementation of multi-way IDR to handle >2 replicates.

Usage
-----

Required input are normalized signal bedGraphs for each replicate, with
separate files for each strand. The bedGraphs can be gzipped if they end in
``.gz``.

One way of doing this might be with deepTools bamCoverage, here, retrieving
unique reads on the minus strand and normalizing using the RPKM method:

::

   bamCoverage \
     --bam rep1.bam \
     -o rep1_minus.bedgraph \
     --outFileFormat bedgraph \
     --binSize 1 \
     --minMappingQuality 20 \
     --samFlagInclude 16 \
     --normalizeUsing RPKM

These bedGraphs can then be provided to termseq-peaks::

   termseq-peaks rep1_minus.bedgraph rep2_minus.bedgraph peaks_minus --strand - 

By default peaks below an IDR threshold of 0.05 will be output to ``peaks_minus.bed`` in narrowPeak format. The full oracle file (the leniently-called peaks on a bedgraph that merges all provided bedgraphs) will be output as ``peaks_minus.bed.oracle.narrowPeak``. 

For more help, run::

   termseq-peaks -h


Algorithm
---------

This tool takes multiple normalized bedGraph files representing the normalized
signal for each replicate, and calls a set of consistent peaks at a provided
IDR [2] cutoff.


- Peaks are called using scipy.signal.find_peaks [2] with very lenient
  parameters to intentionally include both real peaks and noise. These peaks
  are called on each replicate.


- All bedGraphs are additionally merged together and peaks are similarly called
  on that merged signal to get the "oracle" peaks.

- The score for the peaks is the "prominence" value for each peak; see [0] for
  details.

- For each unique combination of replicates, IDR routines from [1] are run,
  resulting in an output file containing merged peaks from those two files
  along with IDR values for each. In practice the tool stores these as temp
  files. The number of peaks falling below the IDR threshold is counted. The
  minimum such number, N, across all pairwise combinations of replicates is
  used as the final number of peaks to select.

- The oracle peaks are then ranked by their score and the top N peaks are
  selected as the final peaks. The scores in the final peaks are the scores
  from the oracle peaks, that is, the peak prominences from calling peaks on
  the merged bedGraphs.

Output
------
The ``find_peaks`` function returns various metrics. Here, we retrieve the
prominence and the width. The prominence is the vertical distance between the
peak and the lowest contour line, and the width is measured at half the
prominence. See these documentation pages for a visualization of these metrics:
`prominences
<https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.peak_prominences.html>`_
and `widths
<https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.peak_widths.html>`_.

Output files are in the `narrowPeak
<https://genome.ucsc.edu/FAQ/FAQformat.html#format12>`_ format, which shows the
peak width as well as the position of the summit. We report the prominence as
the score as well as the signal value. The position of the peak is the 1-bp
position of the prominence.

Caveats
-------
The find_peaks function operates on 1-dimensional vectors, and so returns peak
positions in terms of indexes into the input vectors. Internally, we
interpolate to back-calculate the corresponding genomic coordinates and round
to integers. This may potentially have issues where two peaks that are
genomically far away have nearby indexes (for example, if the intervening
region has zero reads anywhere). Empirically we do not observe this to be an
issue, but a solution would be to pad out the vector to include zeros at every
position in the chromosome/plasmid.

The biggest downside currently is speed and RAM. This is not an issue for the
small bacterial genomes the tool was designed for; it takes about 30s to run
for E. coli data, and NumPy arrays are used to store signal. For larger
eukaryotic genomes, parallelization across chromosomes may be required and
substantial RAM may be required. This tool remains untested on larger genomes,
but has worked quite well for term-seq in several bacterial genomes.

References
----------

[0] https://github.com/macs3-project/MACS/wiki/Advanced%3A-Call-peaks-using-MACS2-subcommands
[1] https://github.com/nboley/idr
[2] https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
