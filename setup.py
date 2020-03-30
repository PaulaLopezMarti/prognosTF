from setuptools import setup, find_packages
from Cython.Build import cythonize


setup(
    name='Meta-Waffle',
    version='0.0.1',
    packages=find_packages() + ['test'],
    ext_modules = cythonize('meta_waffle/*.pyx'),
    include_package_data=True,
    package_dir={'test': 'test'},
    package_data={'test': ['*.pickle', 'data/*.bam*', 'data/*.tsv', 'data/*.bed', 'data/*.pickle']},
    scripts=['scripts/waffle-peaks.py',
             'scripts/waffle-plot.py',
             'scripts/waffle-predict.py',
             'scripts/waffle-merge.py',
             'scripts/waffle-bam2count.py'],
    license='GPLv3'
)
