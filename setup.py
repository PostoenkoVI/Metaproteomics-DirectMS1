#!/usr/bin/env python

'''
setup.py file for metaproteomics script
'''
from setuptools import setup

version = open('VERSION').readline().strip()


setup(
    name                 = 'MetaDirectMS1',
    version              = version,
    description          = '''Ultrafast metaproteomics for quantitative assessment of strain isolates and microbiomes''',
    long_description     = (''.join(open('README.md').readlines())),
    long_description_content_type = 'text/markdown',
    author               = 'Valeriy Postoenko, Kazakova Liza, Ivanov Mark',
    author_email         = 'pyteomics@googlegroups.com',
    install_requires     = [line.strip() for line in open('requirements.txt')],
    classifiers          = ['Intended Audience :: Science/Research',
                            'Programming Language :: Python :: 3.9',
                            'Topic :: Education',
                            'Topic :: Scientific/Engineering :: Bio-Informatics',
                            'Topic :: Scientific/Engineering :: Chemistry',
                            'Topic :: Scientific/Engineering :: Physics'],
    license              = 'License :: OSI Approved :: Apache Software License',
    packages         = ['MetaDirectMS1', ],
    entry_points         = {'console_scripts': ['MetaDirectMS1 = MetaDirectMS1.cli:run',]},
    )
