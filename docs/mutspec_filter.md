### MutSpec-Filter (optional)

Filter out variants that are likely neutral polymorphisms using Annovar annotations (dbSNP, ESP...) or additional databases in [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) or [BED](https://genome.ucsc.edu/FAQ/FAQformat#format1) format.  

#### Input format

Tab delimited text files generated by [MutSpec-Annot](https://github.com/IARCbioinfo/mutspec/blob/master/docs/mutspec_annot.md).

#### Output

Tab delimited text file filtered for variants considered as neutral polymorphisms.

#### Usage

```perl
perl mutspecFilter.pl --thG --filter bed1 --filter bed2 --refGenome hg19 --pathAVDBList mutspec/hg19_listAVDB.txt --outfile output_filename input_file
```

List of parameters:

| Parameter | Default value | Description          |
|-----------|---------------|----------------------|
| --dbSNP   | 0             | Column number for dbSNP (start to count from 1) |
| --segDup  | false         | Remove variants present at >= 0.9 frequency in the genomic duplicate segments database |
| --esp     | false         | Remove variants present at frequency > 0.001 in the Exome Sequencing Project database **(only valid for human genomes)** |
| --thG     | false         | Remove variants present at frequency > 0.001 in the 1000 genome database **(only valid for human genomes)** |
| --exac    | false         | Remove variants present at frequency > 0.001 in the EXome Agregate Consortium database **(only valid for human genomes)** |
| --filter  |  				| Path to one VCF or BED file to filter against (Multiples --filter are allowed) |
| input     |  				| Input file to filter |
| --outfile |  				| Path to output file |
| --help    |               | Print help message |