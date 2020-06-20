# 21cmMG

| 21cmMG:    | A suite for probing modified gravity with 21cm cosmology  |
|------------|-----------------------------------------------------------|
| Author:    | Evan J. Arena                                             |

* (c) Evan J. Arena (Drexel University Department of Physics), 2020.
* For questions, please email `evan.james.arena@drexel.edu.`
* If you wish to use 21cmMG, please cite [arXiv:1810.09572](https://arxiv.org/abs/1810.09572).

## Required Packages

* numpy
* scipy
* matplotlib
* h5py
* pickle
* SNR data cube from Cosmic Visions (optional, see below)
* hi_class Boltzmann code (provided, for convenience, in `21cmMG_hi_class`:
   * [hi_class code and its python wrapper classy](http://miguelzuma.github.io/hi_class_public/)

## Packages and modules

* `21cmMG_FisherForecast`: A Fisher forecast for a general 21cm experiment with a modified gravity parameter-space extension. See [Fisher21cm](https://github.com/evanjarena/Fisher21cm) for specific module details.

* `21cmMG_fsigma8`: examines the effects that Horndeski models of modified gravity/dark energy have on the growth rate of structure fsigma8(z).

* `21cmMG_CMB`: examines the effects that Horndeski models of modified gravity/dark energy have on the CMB lensing potential.

* `21cmMG_fsigma8_and_CMB`: examines the effects that Horndeski models of modified gravity/dark energy have on fsigma8(z), sigma8(z), the CMB lensing potential, and the CMB temperature power spectrum.


