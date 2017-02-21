### MutSpec-Hotspot (optional)

Compute variant frequency on a defined dataset.

The dataset is defined by a collection of sample-specific lists of variants. Frequencies can be calculated on the entire collection, as well as by user-defined categories of samples.

#### Figure 1 - Count variant frequency on a define dataset
![Figure 1](https://github.com/IARCbioinfo/mutspec/blob/master/computefreq.png "Figure 1 - Count variant frequency on a define dataset") 

#### Input format

##### Data
Collection of VCFs (version 4.1 and 4.2) or tab-delimited text files. The accepted format is under the same restriction as MutSpec Annot (see below).

Filenames should not contains "." in their suffix.

Files should contain at least four columns describing for each variant: the chromosome number, the start genomic position, the reference allele and the alternate alleles. These columns can be in any order.

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

For MuTect and MuTect2 output files, only confident calls are considered as other calls are very likely to be dubious calls or artefacts.
Variants containing the string REJECT in the judgement column or not passing MuTect2 filters are not annotated and excluded from MutSpect-Annot output. 

For COSMIC and ICGC files, variants are reported on several transcripts. These duplicate variants need to be remove before annotated the file.

##### InfoFile (optional)
A tabular file associating the name of the different input files of the input collection to sample categories as follows:

###### If using a tumor-normal pair design (with or without tumor replicates)
Name the first column as Normal, the second column will be Tumor. Fill the columns with appropriate file names. You may include unmatched samples, but in this case use “NA” for missing Normal and Tumor files.

| Normal     |   Tumor    | Replicates |
|------------|------------|----------------------|
| Sample_1_N | Sample_1_T | Sample_1_TDup |
| Sample_2_N | Sample_2_T | |
|     NA     | Sample_3_T | |

> This tabular file does not support empty field for Normal and Tumor file, so please name "NA" a missing file.

> The files in the column "Replicates" will be considers as "Tumor" file for computing the variant frequency and count.

###### If using any type of categories (one or more categories)
Organize categories of samples by columns. You may use any names for columns (use preferably short names).

| Category_1 | Category_2 | Category_3 | Category_N |
|------------|------------|------------|------------|
|  Sample_1  |  Sample_2  |  Sample_3  |  Sample_4  |
|  Sample_5  |            |  Sample_7  |  Sample_8  |
|  Sample_9  |            |            |            |

>You can name your column as you want. It supports empty field.

###### Without provided an InfoFile

If you don't provide an InfoFile, you can not make a pair analysis. (see [Figure 1](https://github.com/IARCbioinfo/mutspec/blob/master/computefreq.png))

> Only one category will be consider by the tool to compute the frequency.

##### Pair (boolean)

**If you choose a pair analysis, you must provide an InfoFile.**
>The tool will add an annotation on the statut of the variant as describe in [Figure 3](https://github.com/IARCbioinfo/mutspec/blob/master/annotation.png).

>You must provide an InfoFile that begin with "Normal".

>You can not provide an InfoFile with more than 3 column (Normal - Tumor - Duplicate).

**If you don't choose a pair analysis, you don't need to provide an InfoFile.**
>Even if you provide an InfoFile that begin by "Normal", you will not have any annotation other than the variant frequency if input field paired analysis is "No".

#### Output

HotSpot tool generates two files ([Figure 2](https://github.com/IARCbioinfo/mutspec/blob/master/hotspot.png)) :

**1 - Variants_summary.vcf**:
This file contains **all unique variants detected in the dataset**, annotated with counts and frequencies in each user-defined category and with sample name in which they were found.

**2 - Annotated_dataset**:
This output is a collection that contains **all the input files annotated** with variant frequencies and counts in each user-defined category.
If working with the case scenario "Normal-Tumor-Duplicates", an additional annotation is included on variant germline or somatic status ([Figure 3](https://github.com/IARCbioinfo/mutspec/blob/master/annotation.png)).

#### Figure 2 - HotSpot workflow
![Figure 2](https://github.com/IARCbioinfo/mutspec/blob/master/hotspot.png "Figure 2 - HotSpot workflow") 
>General example of HotSpot application

#### Figure 3 - Annotation for Normal-Tumor(-Replicate) pairs
![Figure 3](https://github.com/IARCbioinfo/mutspec/blob/master/annotation.png "Figure 3 - Annotation for Normal-Tumor(-Replicate) pairs") 

>Example of annotation of variants in the case "Normal-Tumor-Replicates"

>For a variant annotated as Somatic, if it's find in the Tumor twice (Replicates also), it will be annotated as Somatic **confirmed if there is a normal sample for this individual and the variant was not found.**

>For a variant annotated as Somatic, if it's find in the Tumor once (Tumor or Replicate), it will be annotated as Somatic **not confirmed if there is a normal sample for this individual and the variant was not found.**

>For a variant annotated as Somatic, if it's find in the Tumor once (Tumor or Replicate) or twice (Tumor and Replicate), it will be annotated as Somatic **NA if there is no a normal sample for this individual**

#### Usage

```perl
perl mutspecStat.pl -i infoFile.txt -d input_dir/ -o output_dir/ -s variant_summary.vcf -p Y
```
>**All parameters are compulsory.**

List of parameters:

| Parameter | Description                             |
|-----------|-----------------------------------------|
| -i        |  input file of category or "None" if not provided  |
| -d        |  path to the directory of input VCF files          |
| -o        |  path to the directory of output VCF files         |
| -s        |  path to the output file variants_summary.vcf      |
| -p        |  Boolean for paired analysis (Y/N)                 |
| -h        |  describe the help                                 |
