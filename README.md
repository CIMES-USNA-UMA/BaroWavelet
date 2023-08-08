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

## Issues and requests

Please access the following link to create an issue or request:

https://github.com/CIMES-USNA-UMA/BaroWavelet/issues

## Contact information

Email: alvarochaoecija.rprojects@gmail.com

ORCID: https://orcid.org/0000-0002-2691-6936
