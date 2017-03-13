### MutSpec-NMF_Estimate_Signatures

Compute statistics for estimating the number of signatures to extract with [MutSpec-NMF](https://github.com/IARCbioinfo/mutspec/blob/modifs_MA/docs/mutspec_nmf.md).


#### Input format

Input matrix created with the tool [MutSpec-Stat](https://github.com/IARCbioinfo/mutspec/blob/modifs_MA/docs/mutspec_stat.md).  
The input matrix can be found in MutSpec-Stat ouput folder: Mutational_Analysis/Figures/Input_NMF/Input_NMF_Count.txt

#### Output

- Statistics of several approaches used for estimating the number of signatures to extract with NMF:  
	- [Brunet et al.](https://www.ncbi.nlm.nih.gov/pubmed/15016911), proposed to take the first number of signature for which the cophenetic coefficient starts decreasing,
	- [Hutchins et al.](https://www.ncbi.nlm.nih.gov/pubmed/18852176), suggested to choose the first value where the RSS curve presents an inflection point.
	- [Frigyesi et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2623306/), considered the smallest value at which the decrease in the RSS is lower than the decrease of the RSS obtained from random data.  

The estimation are based on Brunet’s algorithm computed from 50 runs for each value of signature to estimate.

The original data are shuffled for comparing the quality measures obtained with our data (Data x) and from randomized data (Data y). The curves for the actual data are in solid line, those for the randomized data are in dashed line.


#### Usage

```R
Rscript R/estimateSign_Galaxy.r --input outfile_MutSpec-Stat/Mutational_Analysis/Figures/Input_NMF/Input_NMF_Count.txt --stop 8 --cpu 8 --output estimate_signatures.png
```

List of parameters:

| Parameter   | Description          |
|-------------|----------------------|
| --input     | Input matrix created with the tool MutSpec-Stat |
| --stop      | Maximum number of signatures to compute (Selecting a number above 8 may not work on small datasets) |
| --cpu       | Number of CPUs |
| --output    | Output figure |
