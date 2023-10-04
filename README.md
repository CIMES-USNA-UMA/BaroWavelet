# BaroWavelet

<!-- badges: start -->
[![R-CMD-check](https://github.com/CIMES-USNA-UMA/BaroWavelet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CIMES-USNA-UMA/BaroWavelet/actions/workflows/R-CMD-check.yaml)
[![DOI:10.1016/j.cmpb.2023.107758](https://img.shields.io/badge/DOI-10.1016/j.cmpb.2023.107758-dcdcdc.svg)](https://doi.org/10.1016/j.cmpb.2023.107758)
<!-- badges: end -->

## Description

Developed by Alvaro Chao-Ecija (2022-2023 Collaboration Scholarship at the 
Department of Physiology and the Autonomic Nervous System Unit at CIMES, University of
Malaga), which started as part of my Final Degree Project, under the supervision of 
MD-PhD. Marc Stefan Dawid Milner.


## Installation

To install the package, use the following code line in R (package remotes is required):

```ruby
remotes::install_github("CIMES-USNA-UMA/BaroWavelet")
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

To cite *BaroWavelet*, use the following citation information:

```ruby
citation("BaroWavelet")
#>To cite BaroWavelet in publications use:
#>
#>  A. Chao-Ecija, M.S. Dawid-Milner, BaroWavelet: An R-based tool for dynamic baroreflex evaluation
#>  through wavelet analysis techniques, Comput Methods Programs Biomed. 242 (2023) 107758.
#>  https://doi.org/10.1016/j.cmpb.2023.107758.
#>
#>A BibTeX entry for LaTeX users is
#>
#>  @Article{CHAOECIJA2023107758,
#>    title = {BaroWavelet: An R-based tool for dynamic baroreflex evaluation through 
#>  wavelet analysis techniques},
#>    author = {A. Chao-Ecija and M.S. Dawid-Milner},
#>    journal = {Computer Methods and Programs in Biomedicine},
#>    year = {2023},
#>    volume = {242},
#>    pages = {107758},
#>    doi = {10.1016/j.cmpb.2023.107758},
#>    issn = {0169-2607},
#>    url = {https://www.sciencedirect.com/science/article/pii/S0169260723004248},
#>  }
```



## Issues and requests

Please access the following link to create an issue or request:

https://github.com/CIMES-USNA-UMA/BaroWavelet/issues

## Contact information

Email: alvarochaoecija.rprojects@gmail.com

ORCID: https://orcid.org/0000-0002-2691-6936
