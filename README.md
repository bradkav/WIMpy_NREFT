# WIMpy_NREFT

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/bradkav/WIMpy_NREFT/master?filepath=NREFT_example.ipynb)

The `WIMpy_NREFT` code allows you to calculate Dark Matter-Nucleus scattering rates in the framework of NREFT (see e.g. [arXiv:1203.3542](https://arxiv.org/abs/1203.3542)). The code is written in python and is a work in progress. It will hopefully be updated soon.

The code currently supports operators <img src="https://rawgit.com/bradkav/WIMpy_NREFT/master/svgs/2330706abb8aba7916b511ca5afa9e62.svg?invert_in_darkmode" align=middle width=19.56603pt height=22.38192pt/> to <img src="https://rawgit.com/bradkav/WIMpy_NREFT/master/svgs/917244ca615745a80feccbe760feb728.svg?invert_in_darkmode" align=middle width=26.09409pt height=22.38192pt/> and contains nuclear form factor information for Xenon, Carbon, Fluorine and Iodine.

For questions, comments or bug reports, please contact Bradley J Kavanagh (bradkav@gmail.com).

### Usage

Requires python 2.7, as well as [NumPy](http://www.numpy.org) and [SciPy](https://www.scipy.org).

To get started, check out the [`NREFT_example.ipynb`](NREFT_example.ipynb), which contains examples of how to use the different parts of the code. The same code is also available as a pure python script in `NREFT_example.py`.

More detailed documentation should be coming soon...


### Notes + Caveats

The code is still a work in progress, so be aware of the following:

- The code assumes a spin-1/2 DM particle.
- The code can only be used for NREFT operators up to <img src="https://rawgit.com/bradkav/WIMpy_NREFT/master/svgs/917244ca615745a80feccbe760feb728.svg?invert_in_darkmode" align=middle width=26.09409pt height=22.38192pt/>.
- The code is *not* optimised for speed.
