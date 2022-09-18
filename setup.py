from distutils.core import setup

requirements = [i.strip() for i in open('requirements.txt')]
setup(
    name='termseq_peaks',
    version='0.2',
    author='Ryan Dale',
    author_email='dalerr@nih.gov',
    license='MIT',
    scripts=['bin/termseq_peaks', 'bin/class_annot', 'bin/table_output'],
    url='https://github.com/nichd-bspc/termseq-peaks',
    packages=['peaklib'],
    install_requires=requirements,
    long_description=open('README.rst').read(),
)
