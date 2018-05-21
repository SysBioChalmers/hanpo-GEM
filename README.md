# Aspergillus_fumigatus-GEM

- Brief Model Description

This repository contains the current genome-scale metabolic model of _Aspergillus fumigatus_. Various different updated versions can be downloaded as [releases](https://github.com/SysBioChalmers/Candida_albicans-GEM/releases).

- Abstract

_Aspergillus fumigatus_ is a fungus of the genus _Aspergillus_, and is one of the most common _Aspergillus_ species to cause disease in individuals with an immunodeficiency. (source: [Wikipedia](https://en.wikipedia.org/wiki/Aspergillus_fumigatus))

- Model KeyWords

**GEM Category:** Species; **Utilisation:** Predictive simulation; **Field:** Metabolic-network reconstruction; **Type of Model:** Reconstruction; **Model Source:** XXX; **Omic Source:** [Genomics](addURLtoGenomePaper)(empty); **Taxonomy:** _Aspergillus fumigatus_; **Metabolic System:** General Metabolism; **Strain:** StrainID; **Condition:** Complex medium;

- Reference: TBA

- Pubmed ID: TBA

- Last update: 2018-05-15

- The model:

|Taxonomy | Template Model | Reactions | Metabolites| Genes |
| ------------- |:-------------:|:-------------:|:-------------:|-----:|
|Aspergillus fumigatus | XXX | TBA | TBA | TBA |


This repository is administered by [@edkerk](https://github.com/edkerk/), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology


## Installation

### Recommended Software:
* A functional Matlab installation (MATLAB 7.3 or higher).
* [RAVEN Toolbox 2](https://github.com/SysBioChalmers/RAVEN) for MATLAB (required for contributing to development). 
* libSBML MATLAB API ([version 5.16.0](https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/MATLAB%20interface/)  is recommended).
* [Gurobi Optimizer for MATLAB](http://www.gurobi.com/registration/download-reg).
* For contributing to development: a [git wrapper](https://github.com/manur/MATLAB-git) added to the search path.

### Installation Instructions
* Clone the [master](https://github.com/SysBioChalmers/Aspergillus_fumigatus-GEM) branch from [SysBioChalmers GitHub](https://github.com/SysBioChalmers).
* Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com).

### Contribute To Development
1. Fork the repository to your own Github account
2. Create a new branch from [`devel`](https://github.com/SysBioChalmers/Aspergillus_fumigatus-GEM/tree/devel).
3. Make changes to the model
    + [RAVEN Toolbox 2](https://github.com/SysBioChalmers/RAVEN) for MATLAB is highly recommended for making changes
    + Before each commit, run in Matlab the `newCommit(model)` function from the `ComplementaryScripts` folder
    + Make a Pull Request to the `devel` folder, including changed `txt`, `yml` and `xml` files

## Contributors
* [Eduard J. Kerkhoven](https://www.chalmers.se/en/staff/Pages/Eduard-Kerkhoven.aspx) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Göteborg, Sweden
* Francisco Zorrilla ([@franciscozorrilla](https://github.com/franciscozorrilla)), Chalmers University of Technology, Göteborg, Sweden 
