# Rampler

[![Build status for c++/clang++](https://travis-ci.org/isovic/rampler.svg?branch=master)](https://travis-ci.org/rvaser/rampler)

Module for sampling genomic sequences.

## Dependencies
1. gcc 4.8+ or clang 3.4+
2. cmake 3.2+

## Installation
To install Rampler run the following commands:

```bash
git clone --recursive https://github.com/rvaser/rampler.git rampler
cd rampler
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

After successful installation, an executable named `rampler` will appear in `build/bin`.

Optionally, you can run `sudo make install` to install rampler executable to your machine.

***Note***: if you omitted `--recursive` from `git clone`, run `git submodule init` and `git submodule update` before proceeding with compilation.

## Usage
Usage of rampler is as following:

    usage: rampler [options ...] <mode>

    <mode>
        subsample <sequences> <reference length> <coverage> [<coverage> ...]

            <sequences>
                input file in FASTA/FASTQ format containing sequences
                to be subsampled
            <reference length>
                integer denoting length of the reference genome (or
                assembly) from which the sequences originate
            <coverage>
                integer denoting desired coverage of the subsampled
                sequences

        split <sequences> <chunk size>

            <sequences>
                input file in FASTA/FASTQ format containing sequences
                which will be split into smaller chunks
            <chunk size>
                integer denoting the desired chunk size in bytes

    options:
        -h, --help
            prints out the help

## Contact information

For additional information, help and bug reports please send an email to one of the following: robert.vaser@fer.hr

## Acknowledgment

This work has been supported in part by Croatian Science Foundation under the project UIP-11-2013-7353.
