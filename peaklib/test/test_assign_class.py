from peaklib import assign_class
from pathlib import Path
import hashlib
import numpy as np

HERE = Path(__file__).resolve().parent

narrowPeak = HERE / "forward.unique.narrowPeak"
bw = [HERE / 'sample1.subset.bam.forward.unique.bigwig',
      HERE / 'sample2.subset.bam.forward.unique.bigwig']
cluster_length = 100
gtf = HERE / 'ecoli-subset.gtf'
#intergenic_min_size = 50
fasta = HERE / 'ecoli_subset.fasta.gz'
param_antisense = 50
param_down = 50
trRNA = HERE / 'trRNA.tsv'
param_upstart = 200
param_downstart = 50
param_upstop = 50
param_fasta = 50
param_downfasta = 10


def test_curatepeaks():
    assign_class.CuratePeaks('asdf', narrowPeak, bw, cluster_length,
                             output=HERE / 'testout',  region='chr:2277000:2287000')
    # possible test is to check md5sums of expected and observed outputs
    obs_md5 = hashlib.md5(open(HERE / 'testout' / 'asdf.raw.tsv','rb').read()).digest()
    exp_md5 = hashlib.md5(open(HERE / 'expected_output' / 'asdf.raw.tsv','rb').read()).digest()
    assert obs_md5 == exp_md5

def test_prepgtfs():
    proc_gtf = assign_class.PrepGtfs(gtf)

    # mean and std of start and end cols, to compare with expected
    assert proc_gtf.orf['start'].mean() == 2281973.6363636362
    assert proc_gtf.orf['start'].std() == 3125.841431446172
    assert proc_gtf.orf['end'].mean() == 2282890.1818181816
    assert proc_gtf.orf['end'].std() == 3438.519123639763
    assert proc_gtf.orf3end['start'].mean() == 2282346.4545454546
    assert proc_gtf.orf3end['start'].std() == 3289.286924658181

def test_classassign():
    peaks = assign_class.CuratePeaks('asdf', narrowPeak, bw, cluster_length,
                                     output=HERE / 'testout',  region='chr:2277000:2287000')
    peak = peaks.curated.set_index("name", drop=False)
    peakdf = peak.copy()
    proc_gtf = assign_class.PrepGtfs(gtf)

    assigned = assign_class.ClassAssign(
        peakdf=peakdf,
        gtfs=proc_gtf,
        fasta=str(fasta),
        param_antisense=param_antisense,
        param_down=param_down,
        trRNA=trRNA,
        param_upstart=param_upstart,
        param_downstart=param_downstart,
        param_upstop=param_upstop,
        param_fasta=param_fasta,
        param_downfasta=param_downfasta,
        output=HERE / 'testout',
        sample='asdf',
        offset=2277000)
    assert ','.join([str(i) for i in assigned.antisense.to_list()]) == 'A:gene_id "EG11419"; gene_name "bcr";,nan,nan,nan,nan'
    assert ','.join([str(i) for i in assigned.internal.to_list()]) == 'nan,I:gene_id "EG10885"; gene_name "rplY";,I:gene_id "EG10885"; gene_name "rplY";,I:gene_id "EG12049"; gene_name "yejM";,nan'
    assert ','.join([str(i) for i in assigned.primary['primary'].tolist()]) == 'nan,nan,nan,nan,P:gene_id "EG30067"; gene_name "proL";'
    assert ''.join([str(i) for i in assigned.details_class['classification'].tolist()]) == 'antisense, internal, internal, internal, primary, '
    assert ','.join([str(i) for i in assigned.upstart_upstop['within_10bp-downstream'].tolist()]) == 'nan,nan,nan,nan,gene_id "EG30067"; gene_name "proL";'
    assert assigned.fasta.loc['chr_2278567_2278567','upstream_3end_seq'] == 'GGTTGACATGCATCAACAAAGGAAGCCTTTTAGCTTCCTCGTTGTGCAATAGATCACCGTT'

def test_assign_fun():
    assign_class.assign(
        sample = 'asdf',
        narrowPeak = narrowPeak,
        bw = bw,
        fasta = str(fasta),
        gtf = gtf,
        trRNA = trRNA,
        output= HERE / 'testout',
        region='chr:2277000:2287000',
        offset=2277000
    )
    obs_md5 = hashlib.md5(open(HERE / 'testout' / 'all.asdf.tsv','rb').read()).digest()
    exp_md5 = hashlib.md5(open(HERE / 'expected_output' / 'all.asdf.tsv','rb').read()).digest()
    assert obs_md5 == exp_md5

def test_table_fun():
    assign_class.table_output(
        sample = 'asdf',
        assigned = HERE / 'expected_output' / 'all.asdf.tsv',
        assigned2 = HERE / 'assigned2.tsv',
        kinefold_scores = HERE / 'kinefold.tsv',
        output= HERE / 'testout'
    )
    obs_md5 = hashlib.md5(open(HERE / 'testout' / 'asdf_TableS1.tsv','rb').read()).digest()
    exp_md5 = hashlib.md5(open(HERE / 'expected_output' / 'asdf_TableS1.tsv','rb').read()).digest()
    assert obs_md5 == exp_md5
