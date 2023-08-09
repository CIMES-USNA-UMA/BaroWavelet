# BaroWavelet

## Description

Developed by Alvaro Chao-Ecija (2022-2023 Collaboration Scholarship at the 
Department of Physiology and the Autonomic Nervous System Unit at CIMES, University of
Malaga), which started as part of my Final Degree Project, under the supervision of 
MD-PhD. Marc Stefan Dawid Milner.


## Installation

To install the package, use the following code line in R (package devtools is required):

```ruby
devtools::install_github("CIMES-USNA-UMA/BaroWavelet")
```

## Shiny application

There is a complementary shiny application available for this package. To access it,
you need to have package *shiny* already installed.

You can use the following BaroWavelet command:

```ruby
BaroWavelet::RunBaroWaveletApp()
```

Alternatively, you can use the following code line (you may need to install *ggplot2*, *ggpubr*, *gridExtra* and *BaroWavelet*):

```ruby
shiny::runGitHub("BaroWaveletApp", "CIMES-USNA-UMA", launch.browser = TRUE)
```
## Citation

An article of this project is in the process of being published and is currently
available as a Journal Pre-proof, with the following citation details (these details
will be updated once the final version of the article is published):

```ruby
citation("BaroWavelet")
#>
#> To cite BaroWavelet in publications use:
#>
#>  Chao-Ecija, A., and Dawid-Milner, M. (2023). BaroWavelet: An R-based tool for
#>  dynamic baroreflex evaluation through wavelet analysis techniques. Comput Methods
#>  Programs Biomed, 107758. doi: 10.1016/J.CMPB.2023.107758.
#>
#> A BibTeX entry for LaTeX users is
#>
#>  @Article{,
#>    title = {BaroWavelet: An R-based tool for dynamic baroreflex evaluation through wavelet analysis techniques},
#>    author = {A Chao-Ecija and MS Dawid-Milner},
#>    journal = {Computer Methods and Programs in Biomedicine},
#>    year = {2023},
#>    pages = {107758},
#>    doi = {https://doi.org/10.1016/j.cmpb.2023.107758},
#>    url = {https://www.sciencedirect.com/science/article/pii/S0169260723004248},
#>  }
```



## Issues and requests

Please access the following link to create an issue or request:

https://github.com/CIMES-USNA-UMA/BaroWavelet/issues

## Contact information

Email: alvarochaoecija.rprojects@gmail.com

ORCID: https://orcid.org/0000-0002-2691-6936
