# MutSpec: mutation spectra analysis tool suite

- [Description](https://github.com/IARCbioinfo/mutspec/tree/modifs_MA#description)
- [Overview of MutSpec tools and workflow](https://github.com/IARCbioinfo/mutspec/tree/modifs_MA#overview-of-mutspec-tools-and-workflow)
- [Getting Started](https://github.com/IARCbioinfo/mutspec/tree/modifs_MA#getting-started)
- [MutSpec tools](https://github.com/IARCbioinfo/mutspec/tree/modifs_MA#mutspec-tools)
	- [MutSpec-Annot](https://github.com/IARCbioinfo/mutspec/blob/modifs_MA/docs/mutspec_annot.md)
	- [MutSpec-Filter](https://github.com/IARCbioinfo/mutspec/blob/modifs_MA/docs/mutspec_filter.md)
	- [MutSpec-Split](https://github.com/IARCbioinfo/mutspec/blob/modifs_MA/docs/mutspec_split.md)
	- [MutSpec-Hotspot](https://github.com/IARCbioinfo/mutspec/blob/modifs_MA/docs/mutspec_hotspot.md)
	- [MutSpec-Stat](https://github.com/IARCbioinfo/mutspec/blob/modifs_MA/docs/mutspec_stat.md)
	- [MutSpec-NMF](https://github.com/IARCbioinfo/mutspec/blob/modifs_MA/docs/mutspec_nmf.md)
	- [MutSpec-Compare](https://github.com/IARCbioinfo/mutspec/blob/modifs_MA/docs/mutspec_compare.md)
- [Download test data](https://github.com/IARCbioinfo/mutspec/tree/modifs_MA#download-test-data)
- [Galaxy installation](https://github.com/IARCbioinfo/mutspec/tree/modifs_MA#galaxy-installation)
- [Authors](https://github.com/IARCbioinfo/mutspec/tree/modifs_MA#authors)
- [Citation](https://github.com/IARCbioinfo/mutspec/tree/modifs_MA#citation)
- [License](https://github.com/IARCbioinfo/mutspec/tree/modifs_MA#license)

## Description

MutSpec includes a suite of tools for analysing and interpreting mutational signatures from next-generation sequencing technologies from human tumours and experimental systems.

MutSpec tools are designed to:
- Annotate genome variations (MutSpec-Annot).
- Filter known polymorphisms and user-defined regions (MutSpec-Filter).
- Split list of variants by sample ID (MutSpec-Split).
- Compute variant frequency on a dataset (MutSpec-HotSpot).
- Compute various statistics describing mutation spectra features (MutSpec-Stat).
- Extract mutational signatures defined by the six types of single base substitutions (SBS) in their trinucleotide sequence context (MutSpec-NMF).
- Compare the obtained signatures with published ones (MutSpec-Compare).

The tools work in a logical sequence, using as input the outputs of each preceding tool, and were developed for human, mouse and rat genomes.

The different tools can be run using command lines and are also available on Galaxy.


## Overview of MutSpec tools and workflow

![workflow](https://github.com/IARCbioinfo/mutspec/blob/modifs_MA/docs/mutspecPipeline.png)

## Getting Started

- [Prerequisites](https://github.com/IARCbioinfo/mutspec/blob/modifs_MA/docs/prerequisites.md)
- [Annovar installation](https://github.com/IARCbioinfo/mutspec/blob/modifs_MA/docs/annovar_installation.md)


## MutSpec tools:

To install MutSpec scripts, clone the github repository.

```bash
git clone https://github.com/IARCbioinfo/mutspec.git
```

Then run the different scripts from mutspec folder.

- [MutSpec-Annot](https://github.com/IARC-bioinfo/mutspec/docs/mutspec_annot.md)
- [MutSpec-Filter](https://github.com/IARC-bioinfo/mutspec/docs/mutspec_filter.md)
- [MutSpec-Split](https://github.com/IARC-bioinfo/mutspec/docs/mutspec_split.md)
- [MutSpec-Hotspot](https://github.com/IARC-bioinfo/mutspec/docs/mutspec_hotspot.md)
- [MutSpec-Stat](https://github.com/IARC-bioinfo/mutspec/docs/mutspec_stat.md)
- [MutSpec-NMF](https://github.com/IARC-bioinfo/mutspec/docs/mutspec_nmf.md)
- [MutSpec-Compare](https://github.com/IARC-bioinfo/mutspec/docs/mutspec_compare.md)


## Download test data

The data and the workflow used in [MutSpec publicaton](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1011-z) are available at 
https://usegalaxy.org/u/maude-ardin/p/mutspectestdata



## Galaxy installation

A package name **mutspec** is available on [Galaxy toolshed](https://toolshed.g2.bx.psu.edu/repository?repository_id=f5c1f75e9fb33f8e) with dedicated instruction informations.


## Authors

Ardin Maude, Cahais Vincent and Robitaille Alexis.

## Citation

Ardin M, Cahais V, Castells X, Bouaoun L, Byrnes G, Herceg Z, Zavadil J and Olivier M (2016). "MutSpec: a Galaxy toolbox for streamlined analyses of somatic mutation spectra in human and mouse cancer genomes." BMC Bioinformatics. [doi: 10.1186/s12859-016-1011-z](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1011-z).

## License

GNU public license version 2 (GPL v2)
