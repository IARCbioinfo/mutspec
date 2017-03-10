### MutSpec-NMF

Extract mutational signatures composed of 96 SBS types (6 SBS types in their trinucleotide sequence context) using the non-negative matrix ([NMF](http://www.nature.com/nature/journal/v401/n6755/full/401788a0.html)) factorisation algorithm of Brunet with the Kullback-Leibler divergence penalty implemented in the [R package NMF](http://www.biomedcentral.com/1471-2105/11/367).

#### Input format

Input matrix created with the tool [MutSpec-Stat](https://github.com/IARCbioinfo/mutspec/blob/modifs_MA/docs/mutspec_stat.md).  
The input matrix can be found in MutSpec-Stat ouput folder: Mutational_Analysis/Figures/Input_NMF/Input_NMF_Count.txt

#### Output

- Composition of mutational signatures and contribution of each signature to each sample.
- Input matrix for the tool [MutSpec-Compare](https://github.com/IARCbioinfo/mutspec/blob/modifs_MA/docs/mutspec_compare.md).


#### Usage

```R
Rscript R/somaticSignature_Galaxy.r --input outfile_MutSpec-Stat/Mutational_Analysis/Figures/Input_NMF/Input_NMF_Count.txt --nbSignature 2 --cpu 8 --output output_dir
```

List of parameters:

| Parameter          | Description          |
|--------------------|----------------------|
| --input            | Input matrix created with the tool MutSpec-Stat |
| --nbSignature      | Number of signatures to extract (min=2) |
| --cpu              | Number of CPUs |
| --output           | Output directory |
| --html             | Path to HTML page (ONLY FOR GALAXY WRAPPER) |
| --help             | Print help message |