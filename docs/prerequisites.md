## Prerequisites

### Build sofwares
```bash
sudo apt-get install build-essential
```

### Build Python modules
```bash
sudo apt-get install python-dev
```

### Perl (>= v5.18.1)
```bash
sudo apt-get install perl
```
- Perl modules:  
[Parallel::ForkManager](http://search.cpan.org/~dlux/Parallel-ForkManager-0.7.5/ForkManager.pm) (>= v0.7.5)  
[Math::Round](http://search.cpan.org/dist/Math-Round/Round.pm) (>= v0.06)  
[Spreadsheet::WriteExcel](http://search.cpan.org/~jmcnamara/Spreadsheet-WriteExcel-2.40/lib/Spreadsheet/WriteExcel.pm) (v2.40)  

These modules can be installed through [CPAN](http://www.cpan.org/)

```bash
sudo apt-get install cpan
cpan
nolock_cpan[1]> install Module::Name
```

### R (v3.2.x)
```bash
sudo apt-get install r-base r-base-dev
```

- R packages  
[getopt](https://cran.r-project.org/web/packages/getopt/index.html) (v1.20.0)  
[lsa](https://cran.r-project.org/web/packages/lsa/index.html) (v0.73.1)  
[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) (v1.0.1)  
[reshape](https://cran.r-project.org/web/packages/reshape/index.html) (>= v0.8.5)  
[NMF](https://cran.r-project.org/web/packages/NMF/index.html) (v0.20.6)  
[gplots](https://cran.r-project.org/web/packages/gplots/index.html) (>= 2.17.0)  
[gtable](https://cran.r-project.org/web/packages/gtable/index.html) (>= v0.1.2)  
[scales](https://cran.r-project.org/web/packages/scales/index.html) (>= 0.2.5)  
[gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) (0.9.1)  
 

These packages can be installed in R:  

```R
install.packages(https://cran.r-project.org/src/contrib/package_1.1.tar.gz, repos=NULL, type="source")
```

### WebLogo3 (>= v3.3)

```bash
sudo pip install weblogo
```

### Annovar (>= 2016Feb01)

The latest release can be downloaded at http://www.openbioinformatics.org/annovar/annovar_download_form.php after free registration.