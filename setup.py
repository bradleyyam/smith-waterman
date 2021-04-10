from setuptools import setup

setup(
    name='swalign',
    version='0.0.1',
    description='An implementation of the Smith-Waterman algorithm for determining optimal local alignment between two sequences given a similarity matrix.',
    url='git@github.com/bradleyyam/smith-waterman.git',
    author='Bradley Yam',
    author_email='bradley.yam@yale.edu',
    license='GPL',
    packages=['swalign'],
    zip_safe=False
)