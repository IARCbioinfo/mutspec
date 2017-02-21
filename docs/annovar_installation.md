## Annovar installation

1. Unpack the package and add ANNOVAR into your system executable path
	```bash
	tar xvfz annovar.latest.tar.gz
	export PATH="/dir/annovar/:$PATH"
	```

2. First you need to create directories for saving Annovar databases files used for the annotation:
	- Create a folder (genomedb) for saving all Annovar databases, e.g. hg19db
	- Create a subfolder (seqFolder) for saving the reference genome, e.g. hg19db/hg19_seq

3. Download the reference genome (by chromosomes) from UCSC for all desired builds as follows:

	```perl
	annotate_variation.pl -buildver <build> -downdb seq <seqFolder>
	```

	`build` is your prefered reference genome  
	`seqFolder` is the location where the sequences (by chromosmes) should be stored, e.g. hg19db/hg19_seq

4. Download all desired databases for your builds. **At least the database refGene must be downloaded.**  

	```perl
	annotate_variation.pl -buildver <build> -webfrom annovar -downdb <database> <genomedb>
	```
	`build` is your prefered reference genome  
	`database` is the name of database file to download, e.g. refGene  
	`genomedb` is the location where all database files should be stored, e.g. hg19db  

	The list of available Annovar databases can be found here: http://annovar.openbioinformatics.org/en/latest/user-guide/download/


5. Edit the file build_listAVDB.txt present in mutspec folder to reflect the name and the type of the databases you have installed