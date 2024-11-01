# diffusr <img src="https://rawgit.com/dirmeier/diffusr/master/inst/fig/diffusion.gif" align="right" width="160px"/>

[![Project Status](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)
[![Build Status](https://travis-ci.org/dirmeier/diffusr.svg?branch=master)](https://travis-ci.org/dirmeier/diffusr)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/dirmeier/diffusr?branch=master&svg=true)](https://ci.appveyor.com/project/dirmeier/diffusr)
[![codecov](https://codecov.io/gh/dirmeier/diffusr/branch/master/graph/badge.svg)](https://codecov.io/gh/dirmeier/diffusr)
[![CRAN](http://www.r-pkg.org/badges/version/diffusr?color=brightgreen)](https://cran.r-project.org/package=diffusr)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/diffusr?color=brightgreen)](https://cran.r-project.org/package=diffusr)

Network diffusion algorithms in R.

## Introduction

`diffusr` implements several algorithms for network diffusion such as *Markov random walks with restarts* and *weighted neighbor classification*. Network diffusion has been studied extensively in bioinformatics, e.g. in the field of cancer gene prioritization. Network diffusion algorithms generally spread information, e.g. encoded as node weights, along the edges of a graph to other nodes. These weights can for example be interpreted as temperature, an initial amount of water, the activation of neurons in the brain, or the location of a random surfer in the internet. The information (node weights) is iteratively propagated to other nodes until a equilibrium state or stop criterion occurs.

## Before installation

Before installation, we recommended you install Intel oneAPI Math Kernel Library (oneMKL) to optimize the computational performance of linear algebra.

Windows users can download oneMKL from [Intel's website](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) and install it in the default directory. The default directory is: `C:\Program Files (x86)\Intel\oneAPI`.

Debian users can download oneMKL using apt in the Debian non-free repo:

``` bash
# Install oneMKL version 2020.4.304-4
sudo apt install intel-mkl-full
```

Or using the Intel repo:

``` bash
# Set up the repository and signed the entry
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
| gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
# Update the package list
sudo apt update
# Install the latest oneMKL (version 2024.0)
sudo apt install intel-oneapi-mkl
```

Fedora users can download oneMKL by using dnf:

``` bash
# Create dnf repository file
tee > /tmp/oneAPI.repo << EOF
[oneAPI]
name=Intel® oneAPI repository
baseurl=https://yum.repos.intel.com/oneapi
enabled=1
gpgcheck=1
repo_gpgcheck=1
gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
EOF
sudo mv /tmp/oneAPI.repo /etc/yum.repos.d
# Install the latest oneMKL (version 2025.0)
sudo dnf install intel-oneapi-mkl intel-oneapi-mkl-devel intel-oneapi-mkl-core
```

## Installation

Install `diffusr` from GitHub:

``` r
if ("devtools" %in% installed.packages()[, "Package"])
  install.packages("devtools")
remotes::install_github("randef1ned/diffusr", upgrade = "always", build_vignettes = TRUE, build_manual = TRUE)
```

## Usage

Load the package using `library(diffusr)`. We provide a vignette for the package that can be called using: `vignette("diffusr")`.
Basically that is all you have to know.


## Author

* Simon Dirmeier <a href="simon.dirmeier@gmx.de">simon.dirmeier@gmx.de</a>
