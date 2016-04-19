==============================
          MutSpec-Suite        
==============================

Created by Maude Ardin and Vincent Cahais (Mechanisms of Carcinogenesis Section, International Agency for Research on Cancer F69372 Lyon France, http://www.iarc.fr/)

Version 1.0

Released under GNU public license version 2 (GPL v2)

Package description: Ardin et al. - 2016 - MutSpec: a Galaxy toolbox for streamlined analyses of somatic mutation spectra in human and mouse cancer genomes - BMC Bioinformatics

Test data: https://usegalaxy.org/u/maude-ardin/p/mutspectestdata



### Requirements

	# python-dev
build-essential and python-dev packages must be installed on your machine before installing MutSpec tools:
$ sudo apt-get install build-essential python-dev


	# Annovar
If you do not have ANNOVAR installed, you can download it here: http://www.openbioinformatics.org/annovar/annovar_download_form.php

1) Once downloaded, install annovar per the installation instructions and edit the PATH variable in galaxy deamon (/etc/init.d/galaxy) to reflect the location of directory containing perl scripts.

2) Create directories for saving Annovar databases
	2-a Create a folder (annovardb) for saving all Annovar databases, e.g. hg19db
	2-b Create a subfolder (seqFolder) for saving the reference genome, e.g. hg19db/hg19_seq

3) Download the reference genome (by chromosome) from UCSC for all desired builds as follows:
$ annotate_variation.pl -buildver <build> -downdb seq <seqFolder>

where <build> can be hg18, hg19 or hg38 for the human genome or mm9, mm10 for the mouse genome.
and <seqFolder> is the location where the sequences (by chromosme) should be stored, e.g. hg19db/hg19_seq


4) Download all desired databases for all desired builds as follows:
$ annotate_variation.pl -buildver <build> [-webfrom annovar] -downdb <database> <annovardb>

/!\ At least the database refGene must be downloaded /!\

where <build> can be hg18, hg19 or hg38 for the human genome or mm9, mm10 for the mouse genome.
and <database> is the database file to download, e.g. refGene
and <annovardb> is the location where all database files should be stored, e.g. hg19db

The list of all available databases can be found here: http://annovar.openbioinformatics.org/en/latest/user-guide/download/


5) Edit the annovar_index.loc file (in the folder galaxy-dist/tool-data/toolshed/repos/iarc/mutspec/revision/) to reflect the location of annovardb folder (containing all the databases files downloaded from Annovar).
Restart galaxy instance for changes in .loc file to take effect or reload it into the admin interface.

6) Edit the file build_listAVDB.txt in the mutspec install directory to reflect the name and the type of the databases installed


### Installation

	# MutSpec-Stat and MutSpec-NMF
By default 8 CPUs are used by these tools, but you may edit mutspecStat_wrapper.sh and mutspecNmf_wrapper.sh to change this number to the maximum number of CPU available on your server.

MutSpec-Stat and MutSpec-NMF tools allow parallel computations that are time consuming.
It is recommended to use the highest number of cores available on the Galaxy server to reduce the computation time of these tools.




	# MutSpec-Annot
The maximum CPU value needs to be specified when installing MutSpec package by editing the file mutspecAnnot.pl to reflect the maximum number of CPU available on your server.

This tool may be time consuming for large files. For example, annotating a file of more than 25,000 variants takes 1 hour using 1 CPU (2.6 GHz), while annotating this file using 8 CPUs takes only 5 minutes.
We have optimized MutSpec-Annot so that the tool uses more CPUs, if available, as follows:
-files with less than 5,000 lines: 1 CPU is used
-files with more than 5,000 and less than 25,000 lines: 2 CPUs are used
-files with more than 25,000 and less than 100,000 lines: 8 (or maximum CPUs, if less than 8 CPUs are available) are used (our benchmark results didn't show any time saving using more than 8 cores for files with more than 25,000
but less than 100,000 lines)
-files with more than 100,000: maximum CPUs are used 
