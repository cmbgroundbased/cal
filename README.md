The CAL library is born from the needs of the CMB ground-based telescope simulation framework to have a module that can take into account the atmospheric effects. Different experiments are lead by different people that choose different programming solutions to implement the characteristic of their instruments.

In general, this approach is not particularly problematic even if some guys working for multiple frameworks that are written not in the same code. In this case, the porting of parts of code, that they do the same things, is not easy and the result is a colossal waste of time spent, for first to get used to with the framework, maybe with the language and finally to re-implement the same algorithm within it.

One of these tasks is represented by atmospheric simulations. The atmospheric time evolution and emission proprieties are the same for the same locations and are not up to the instrumental feature.

Out of this consideration, we decided to extract the code that does the atmospheric simulation from the TOAST framework. This bunch of code is released by (Who actually do release this?) under the <LICENSE> license.

The maintenance and future improvements actions of the CAL library are up to the CAL authors, that are listed in the AUTHORS file in this repository.

## Install

## Contribution guidelines

## Continuous build status

### Official Builds

### Community Supported Builds

## Resources

## License
