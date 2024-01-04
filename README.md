# hanpo-GEM: The genome-scale metabolic model of _Hansenula polymorpha_

# Abstract

Genome-scale metabolic models (GEMs) provide a useful framework for modeling the metabolism of microorganisms. While the applications of GEMs are wide and far reaching, the reconstruction and continuous curation of such models can be perceived as a tedious and time-consuming task. Using RAVEN, a MATLAB-based toolbox designed to facilitate the reconstruction analysis of metabolic networks, this protocol practically demonstrates how researchers can create their own GEMs using a homology-based approach. To provide a complete example, a draft GEM for the industrially relevant yeast _Hansenula polymorpha_ is reconstructed.

# Description

This repository contains the current genome-scale metabolic model of _Hansenula polymorpha_, synonymously known as _Ogataea polymorpha_. _Hansenula/Ogataea polymorpha_ is a filamentous yeast from the family _Saccharomycetaceae_, and is an industrially relevant methylotrophic species. 

Clone this repo to download the model and associated data:

```
$ git clone https://github.com/SysBioChalmers/hanpo-GEM.git
```

# Citation

If you use hanpo-GEM, or the reconstruction protocol used for its generation, please cite the following textbook chapter:

  > Zorrilla, F., Kerkhoven, E.J. (2022). Reconstruction of Genome-Scale Metabolic Model for Hansenula polymorpha Using RAVEN. In: Mapelli, V., Bettiga, M. (eds) Yeast Metabolic Engineering. Methods in Molecular Biology, vol 2513. Humana, New York, NY. https://doi.org/10.1007/978-1-0716-2399-2_16

# Keywords

**Utilisation:** predictive simulation\
**Field:** metabolic-network reconstruction\
**Type of Model:** homology-based reconstruction\
**Model Source:** [hanpo-GEM](https://doi.org/10.1007/978-1-0716-2399-2_16)\
**Omic Source:** [genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001664045.1/);[protein](https://mycocosm.jgi.doe.gov/Hanpo2/Hanpo2.home.html)\
**Taxonomy:** _Hansenula polymorpha_/_Ogataea polymorpha_\
**Metabolic System:** General Metabolism\
**Strain:** NCYC 495 leu1.1\
**Condition:** Complex medium

# Model overview

|Taxonomy | Template Model | Reactions | Metabolites| Genes |
| ------------- |:-------------:|:-------------:|:-------------:|-----:|
|_Hansenula polymorpha_| [yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM) & [rhto-GEM](https://github.com/SysBioChalmers/rhto-GEM/) | 2370 | 2118 | 984 |

# Dependencies

If you want to use the model for your own model simulations, you can use **any software** that accepts SBML L3V1 FBCv3 formatted model files. This includes any of the following:
* MATLAB-based
  * [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) version 2.8.3 or later (recommended)  
  * [COBRA Toolbox](https://github.com/opencobra/cobratoolbox)

* Python-based
  * [cobrapy](https://github.com/opencobra/cobrapy)  

Please see the installation instructions for each software package.

# Contributing

Contributions are always welcome! Please read the [contributions guideline](https://github.com/SysBioChalmers/yeast-GEM/blob/main/.github/CONTRIBUTING.md) to get started.

# Contributors
* [Eduard J. Kerkhoven](https://www.chalmers.se/en/staff/Pages/Eduard-Kerkhoven.aspx) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Göteborg, Sweden
* [Francisco Zorrilla](https://www.mrc-tox.cam.ac.uk/staff/francisco-zorrilla) ([@franciscozorrilla](https://github.com/franciscozorrilla)), MRC Toxicology Unit, University of Cambridge, UK
