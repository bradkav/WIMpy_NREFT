# WIMpy_NREFT

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/bradkav/WIMpy_NREFT/master?filepath=Examples/NREFT_example.ipynb)

The `WIMpy_NREFT` code allows you to calculate Dark Matter-Nucleus scattering rates in the framework of NREFT (see e.g. [arXiv:1203.3542](https://arxiv.org/abs/1203.3542)). The code is written in python, with more detailed documentation to appear soon.

The code currently supports operators <img src="https://rawgit.com/bradkav/WIMpy_NREFT/master/svgs/2330706abb8aba7916b511ca5afa9e62.svg?invert_in_darkmode" align=middle width=19.56603pt height=22.38192pt/> to <img src="https://rawgit.com/bradkav/WIMpy_NREFT/master/svgs/917244ca615745a80feccbe760feb728.svg?invert_in_darkmode" align=middle width=26.09409pt height=22.38192pt/> and contains nuclear form factor information for Xenon, Carbon, Fluorine and Iodine.

For questions, comments or bug reports, please contact Bradley J Kavanagh (bradkav@gmail.com).

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
-  [`NREFT_example.ipynb`](Examples/NREFT_example.ipynb), which contains examples of how to use the different parts of the code, including calculating a range of spectra.  
- [`Spectra.ipynb`](Examples/Spectra.ipynb), which can be used to generate plots of spectra for all NREFT operators and a range of experiments.

More detailed documentation should be coming soon...


### Notes + Caveats

The code is still a work in progress, so be aware of the following:

- The code assumes a spin-1/2 DM particle.
- The code can only be used for NREFT operators up to <img src="https://rawgit.com/bradkav/WIMpy_NREFT/master/svgs/917244ca615745a80feccbe760feb728.svg?invert_in_darkmode" align=middle width=26.09409pt height=22.38192pt/>.
