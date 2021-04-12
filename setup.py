#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.md') as history_file:
    history = history_file.read()

requirements = [ ]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]

setup(
    author="Alex Orlek",
    author_email='alex.orlek@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="figleaf_fasta applies hard/soft masking to a FASTA file or excludes/extracts sub-sequences from a FASTA file.",
    install_requires=['biopython>=1.61'],
    license="MIT license",
    long_description_content_type='text/markdown',
    include_package_data=True,
    keywords='figleaf_fasta',
    name='figleaf_fasta',
    packages=find_packages(include=['figleaf_fasta', 'figleaf_fasta.*']),
    test_suite='tests',
    url='https://github.com/AlexOrlek/figleaf_fasta',
    version='1.1.0',
    zip_safe=False,
    entry_points={
        'console_scripts': ['figleaf=figleaf_fasta.figleaf:main']
    },
)
