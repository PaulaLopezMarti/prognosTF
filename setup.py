from setuptools import setup, find_packages
from Cython.Build import cythonize


setup(
    name='Meta-Waffle',
    version='0.0.1',
    packages=find_packages(),
    ext_modules = cythonize('meta_waffle/*.pyx'),
    scripts=['scripts/waffle-peaks.py',
             'scripts/waffle-plot.py',
             'scripts/waffle-predict.py',
             'scripts/waffle-merge.py',
             'scripts/waffle-bam2count.py'],
    license='GPLv3'
)
