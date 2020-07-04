The CMB Atmospheric Library (CAL) project is born from the needs of the CMB ground-based telescope simulation framework to have a module that can take into account the atmospheric effects. Different experiments are led by various people that choose different programming solutions to implement their instruments' characteristics.

In general, this approach could be problematic for some guys that are working for multiple frameworks (that are written not in the same languages) ant they have to implement the same part of code. 

In this case, the porting of parts of code, that they do the same things, is not easy and the result is a colossal waste of time spent, for first, to get used to with the framework, maybe with the language and finally to re-implement the same algorithm within it.

Atmospheric simulations represent one of these tasks. The atmospheric time evolution and emission proprieties are the same for the same locations and are not up to the instrumental features.

Out of this consideration, we decided to make independent the code that does the atmospheric evolution and observation, and it is already provided by the TOAST framework. This bunch of code is released by  under the BSD/MIT license.

The atmospheric model is based on the S. Church 1995 paper (https://doi.org/10.1093/mnras/272.3.551). The mathematical description is available on @ref church_model

## Build status

| Distribution      | Status |
|:-----------:|:------:|
| Linux/Amd64 |  ![C/C++ CI](https://github.com/cmbgroundbased/libcal/workflows/C/C++%20CI/badge.svg) |   |   |   |


## Intall the C++ library and Python bindings

Work in progress ...

## AUTHOR

HPC4CMB/TOAST Author <br />
Theodore Kisner <work@theodorekisner.com> <br />
Reijo Keskitalo <reijo.keskitalo@gmail.com> <br />
Andrea Zonca <zonca@sdsc.edu> <br />
Giuseppe Puglisi <giuse.puglisi@gmail.com> <br />

CAL Maintainer <br />
Stefano Mandelli <stefano.mandelli@unimi.it>


## License


Time Ordered Astrophysics Scalable Tools (TOAST)

Copyright (c) 2015-2018, The Regents of the University of California
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

