### MutSpec-Split (optional)

Split list of variants by sample ID.  
**This tool need to be used if you have a single file containing all your samples.**

#### Input format

A tab delimited text file, with a column containing the samples ids.

#### Output

A folder containing tab delimited text files for each samples present in the input file.

#### Usage

```perl
perl mutspecSplit.pl --file filename --column index
```

List of parameters:

| Parameter | Description          |
|-----------|----------------------|
| --file    | Input file to split |
| --column  | Column number containing the sample ids (start to count from 1) |