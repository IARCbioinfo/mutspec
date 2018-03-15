### MutSpec-Compare

Computes the similarity between the signatures identified by MutSpec-NMF and a set of published signatures using the cosine similarity method.

#### Input format

- A matrix containing the signatures identified by MutSpec-NMF.  
The input matrix can be found in MutSpec-NMF ouput folder: NMF/Files/MatrixW-Normto100.txt  

- A matrix containing a set of published signatures.  

Two matrices of published signatures are provided:
- `COSMIC30-Hupki-Others.txt`: Contains the 30 mutational signatures operative in human cancer from [COSMIC](http://cancer.sanger.ac.uk/cosmic/signatures), plus 4 experimental signatures obtained in Hupki mouse cells after exposure to AA, AID, BaP, MNNG published by [Olivier et al.](http://www.nature.com/articles/srep04482), plus addition experimental signatures obtained in mouse models for DMBA published by [Nassar et al.](http://www.nature.com/nm/journal/v21/n8/full/nm.3878.html), and MNU and Urethane published by [Westcott et al.](http://www.nature.com/nature/journal/v517/n7535/full/nature13898.html).
- `sigProfiler60-Hupki-Others.txt`: Contains 60 mutational signatures identified in human cancer by the PCAWG project with the sigProfiler algorithm (https://fr.mathworks.com/matlabcentral/fileexchange/38724-sigprofiler), plus the experimental signatures described above.

#### Output

Heatmap showing the cosine similarity between the signatures identified by MutSpec-NMF and a set of published signatures.


#### Usage

```R
Rscript R/compareSignature_Galaxy.r Frequency-COSMICv72-Hupki.txt outfolder_MutSpec-NMF/NMF/Files/MatrixW-Inputggplot2.txt output_dir_MutSpec-NMF/NMF/
```

List of parameters:

| Parameter           | Description          |
|---------------------|----------------------|
| Published_Signature | Matrix containing a set of published signatures |
| New_Signature       | Matrix containing the signatures identified by MutSpec-NMF |
| Output_Folder       | Output directory |
