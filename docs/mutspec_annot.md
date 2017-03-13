### MutSpec-Annot

Provides functional annotations from [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/) software, as well as the strand transcript orientation (from refGene database) and the sequence context of variants (extrated from the reference genome selected).

#### Input format

VCF(s) (version 4.1 and 4.2) or tab-delimited text file from various variant callers.

Filenames must be &#60;= 31 characters and should not contains "." in their suffix.

Files should contain at least four columns describing for each variant: the chromosome number, the start genomic position, the reference allele and the alternate alleles. These columns can be in any order and other columns may be present.

If multiple input files are specified as input for the tool, they should be from the **same genome build** and in the **same format**.

The tool supports different column names (**names are case-sensitive**) depending on the source file as follows:

- `mutect`:     contig position ref_allele alt_allele

- `vcf`:        version [4.1](https://samtools.github.io/hts-specs/VCFv4.1.pdf) and [4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf)

- `cosmic`:     Mutation_GRCh37_chromosome_number Mutation_GRCh37_genome_position Description_Ref_Genomic Description_Alt_Genomic

- `icgc`:       chromosome chromosome_start reference_genome_allele mutated_to_allele

- `tcga`:       Chromosome Start_position Reference_Allele Tumor_Seq_Allele2

- `ionTorrent`: chr Position Ref Alt            

- `proton`:     Chrom Position Ref Variant        

- `varScan2`:   Chrom Position Ref VarAllele

- `varScan2 somatic`:   chrom position ref var

- `annovar`:    Chr Start Ref Obs                 

- `custom`:     Chromosome Start Wild_Type Mutant


For COSMIC and ICGC files, variants are reported on several transcripts. These duplicate variants need to be removed before annotating the file.


#### Output

The output is a tabular text file containing the retrieved annotations in the first columns and all columns from the original file at the end.

Only classic chromosomes are considered for the annotation, all other chromosomes are excluded from MutSpec-Annot output.
For example for human genome only chr1 to chrY are annotated.

For MuTect and MuTect2 output files, only confident calls are considered as other calls are very likely to be dubious calls or artefacts.
Variants containing the string REJECT in the judgement column or not passing MuTect2 filters are not annotated and excluded from MutSpect-Annot output. 


The minimum annotations retrieved are:

- Gene-based: RefSeqGene

- Transcript orientation

The strand annotation corresponding to transcript orientation within genic regions is recovered from RefSeqGene database.

- Sequence context

Flanking bases in both sides in 5' and 3' of the variant position are retrieved from the reference genome used.


#### Usage

```perl
perl mutspecAnnot.pl --refGenome hg19 --pathAnnovarDB path/hg19db --pathAVDBList mutspec/hg19_listAVDB.txt input
```

`--refGenome`, `--pathAnnovarDB`, `--pathAVDBList` and `input` are compulsory.  

List of parameters with default values:

| Parameter         | Default value                        | Description          |
|-------------------|--------------------------------------|----------------------|
| --refGenome       |                                      | Name of the reference genome |
| --pathAnnovarDB   |                                      | Path to annovar databases files |
| --pathAVDBList    |                                      | Path to build_listAVDB.txt |
| input             |                                      | Input folder or input file |
| --interval        | 10                                   | Number of bases to retrieve for the sequence context |     
| --outfile         | input directory                      | Output directory |
| --temp            | directory in which the script is run | Path for saving temporary files |
| --fullAnnotation  | yes                                  | Recover all Annovar annotations specified in build_listAVDB.txt |
| --max_cpu         | 8                                    | Number of CPUs to be used for the annotation |
| --help            |                                      | Print help message |


##### Recommendation for setting the parameters:

- --fullAnnotation

Set this parameter to "no" if your analysis includes millions of variants and you are just interested in having a quick overview of the mutation spectrum. Only the annotations needed for the tool [MutSpec-Stat](https://github.com/IARCbioinfo/mutspec/blob/master/docs/mutspec_stat.md) will be added (annotation from refGene, the strand orientation and the sequence context).

- --max_cpu

Annotating large files may be time consuming. For example, annotating a file of more than 25,000 variants takes 1 hour using 1 CPU (2.6 GHz),
while annotating this file using 8 CPUs takes only 5 minutes.  

Recommanded number of CPUs to use:  
-files with less than 5,000 lines: 1 CPU  
-files with more than 5,000 and less than 25,000 lines: 2 CPUs  
-files with more than 25,000 and less than 100,000 lines: 8 CPUs (our benchmark
results didn't show any time saving using more than 8 CPUs for files with more than 25,000 but less than 100,000 lines)  
-files with more than 100,000: maximum CPUs available on your machine