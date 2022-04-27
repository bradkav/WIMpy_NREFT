# WIMpy_NREFT

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/bradkav/WIMpy_NREFT/master?filepath=Examples/NREFT_example.ipynb) [![DOI](https://zenodo.org/badge/98175259.svg)](https://zenodo.org/badge/latestdoi/98175259) [![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php)

The `WIMpy_NREFT` code (also known as just `WIMpy`) allows you to calculate Dark Matter-Nucleus scattering rates in the framework of NREFT (see e.g. [arXiv:1203.3542](https://arxiv.org/abs/1203.3542)). The code is written in python, with more detailed documentation to appear soon.

The code currently supports operators <img src="https://rawgit.com/bradkav/WIMpy_NREFT/master/svgs/2330706abb8aba7916b511ca5afa9e62.svg?invert_in_darkmode" align=middle width=19.56603pt height=22.38192pt/> to <img src="https://rawgit.com/bradkav/WIMpy_NREFT/master/svgs/917244ca615745a80feccbe760feb728.svg?invert_in_darkmode" align=middle width=26.09409pt height=22.38192pt/>, as well as millicharged and magnetic dipole Dark Matter. The code can be used to generate spectra for Xenon, Argon, Carbon, Germanium, Iodine and Fluorine targets.

`WIMpy_NREFT` now includes functionality to calculate *directional* recoil spectra, as well as signals from coherent neutrino-nucleus scattering (including fluxes from the Sun, atmosphere and diffuse supernovae).

**Authors:** Bradley J Kavanagh, Tom D P Edwards.

For questions, comments or bug reports, please contact Bradley J Kavanagh (bradkav@gmail.com).

**Updates:**
* 27/04/2022: Fixed some bugs with incorrect nuclear response functions (this should only affect relatively rare parameter combinations for certain elements)
* 06/04/2022: Added new targets to Nuclei.txt.
* 29/09/2021: **Version 1.1** - Some operators were missing powers of (q/mN)^2 in the rate calculation, which has now been corrected.  

### Installation

You can install `WIMpy_NREFT` using `pip`:

```
pip install git+https://github.com/bradkav/WIMpy_NREFT
```

Requires python3 as well as [NumPy](http://www.numpy.org) and [SciPy](https://www.scipy.org).

### Usage

Most of the relevant routines are contained in the module `DMUtils.py`. Load with

```python
from WIMpy import DMUtils as DMU
```

For how to use the routines, there are a number of examples in the `Examples/` folder:

*  [`NREFT_example.ipynb`](Examples/NREFT_example.ipynb), which contains examples of how to use the different parts of the code, including calculating a range of spectra.  
* [`Spectra.ipynb`](Examples/Spectra.ipynb), which can be used to generate plots of spectra for all NREFT operators and a range of experiments.  
* [`Directional.ipynb`](Examples/Directional.ipynb), which demonstrates how to calculate *directional* recoil spectra, as well as how to transform into different coordinate systems and account for *time-integrated* directionality.  
* [`Neutrinos.ipynb`](Examples/Neutrinos.ipynb), which shows how to calculate neutrino-nucleus scattering spectra.


### Notes + Caveats

The code is still a work in progress, so be aware of the following:

- The code assumes a spin-1/2 DM particle.
- The code can only be used for NREFT operators up to <img src="https://rawgit.com/bradkav/WIMpy_NREFT/master/svgs/917244ca615745a80feccbe760feb728.svg?invert_in_darkmode" align=middle width=26.09409pt height=22.38192pt/>.

### Citation

If you use the WIMpy code, please cite it as
```
B. J. Kavanagh and T. D. P. Edwards, WIMpy NREFT v1.0 [Computer Software], doi:10.5281/zenodo.1230503. Available at https://github.com/bradkav/WIMpy_NREFT, (2018)
```
The corresponding bibtex is:
```
@misc{WIMpy-code,
author = {Kavanagh, Bradley J. and Edwards, Thomas D. P.},
title = {\textnormal{WIMpy\_NREFT v1.1 [Computer Software]}, \href{https://doi.org/10.5281/zenodo.1230503}{\textnormal{doi:10.5281/zenodo.1230503}}\textnormal{. Available at }\url{https://github.com/bradkav/WIMpy_NREFT}},
year = {2018}
}
```

### Publications

The code has been used in a number of publications, including:

- *Results on photon-mediated dark matter-nucleus interactions from the PICO-60 C3F8 bubble chamber*, Ali et al. (PICO Collaboration, 2022), [arXiv:2204.10340](https://arxiv.org/abs/2204.10340)  
- *Constraints on dark matter-nucleon effective couplings in the presence of kinematically distinct halo substructures using the DEAP-3600 detector*, Adhikari et al., (DEAP-3600 Collaboration, 2020), [arXiv:2005.14667](https://arxiv.org/abs/2005.14667)  
- *Digging for dark matter: Spectral analysis and discovery potential of paleo-detectors*, Edwards et al. (2019), [arXiv:1811.10549](https://arxiv.org/abs/1811.10549)  
- *Dark Matter Model or Mass, but Not Both: Assessing Near-Future Direct Searches with Benchmark-free Forecasting*, Edwards, Kavanagh & Weniger (2018), [arXiv:1805.04117](https://arxiv.org/abs/1805.04117)

