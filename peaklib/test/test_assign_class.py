from peaklib import assign_class
from pathlib import Path

HERE = Path(__file__).resolve().parent

narrowPeak = HERE / "forward.unique.narrowPeak"
bw = [HERE / 'sample1.subset.bam.forward.unique.bigwig',
      HERE / 'sample2.subset.bam.forward.unique.bigwig']
cluster_length = 100

def test_curatepeaks():
    assign_class.CuratePeaks('asdf', narrowPeak, bw, cluster_length, output=HERE / 'testout',  region='chr:2277000:2287000')
