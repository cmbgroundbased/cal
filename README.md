# CMB Atmspheric Library (CAL)

The CAL code is extracted from the famous [TOAST framework](https://github.com/hpc4cmb/toast) with the goal to make indipendend the atmospheric effects simulation for CMB ground-based telescope.

The project is born from the needs of the CMB ground-based telescope simulation framework to have a module that can take into account the atmospheric effects. Different experiments are led by various people that choose different programming solutions to implement their instruments' characteristics that, more often that not, are particulary exotic!

In general, this approach could be problematic for those DEVs that are working for multiple frameworks (that are written not in the same languages) and they have to implement the same part of code.

The atmospheric simulations represent one of these tasks. The atmospheric time evolution and emission proprieties are the same for the same locations and are not up to the instrumental features.

Out of this consideration, we decided to make independent the code that does the atmospheric evolution and observation, and it is already provided by the TOAST framework, with python binding (in order to make easy the implementation in your instrumental python framework) and for the future we are planning to release wrappers also for julia (>1.0) language.

The atmospheric model is based on the S. Church 1995 paper (https://doi.org/10.1093/mnras/272.3.551). The mathematical description is available on @ref church_model 

# Build status and informations

<center>

<table>
  <tr><th>Build status</th><th>Documentation</th><th>Code coverage</th></tr>
<tr><td>
  
| Distribution| Status |  
|:-----------:|:------:| 

| Singularity Build | ![Singularity Build (docker)](https://github.com/cmbgroundbased/cal/workflows/Singularity%20Build%20(docker)/badge.svg?branch=master) | 

</td><td>

| Version      | Status |
|:-----------:|:------:|
| 0.9 |  [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cmbgroundbased.github.io/cal/)|


</td><td>
  
 |             |             |
 |:-----------:|:-----------:|
 |             |             |

</tr> </table>



</center>




# Install the C++ library and Python bindings

Clone the repository and the submodules

`git https://github.com/cmbgroundbased/cal` <br/>
`cd cal`<br/>
`git submodule update --init --recursive`<br/>

Build the project

`mkdir build`<br />

If(you don't need `Julia` wraps)<br/>
	&nbsp;&nbsp;`cd build; cmake ../`<br />
else<br/>
	&nbsp;&nbsp;`cd build; cmake -DJULIA=Y ../`<br/> 

For Python user install:<br/>
        &nbsp;&nbsp;`cd build; cmake -DPYTHON_USER_INSTALL=Y ../`<br />

`make -J <N>`<br />
`make install`<br />

## CAL requirements:

`openompi (>= 4.0.0)` <br/>
`libaatm`<br/>
`SuiteSparse`<br/>
`LAPACK`<br/>
`python3`

The `Julia` and `Python` wrappers are provided with `libcxxwraps-julia` and `pybind11`, respectively. These two packages has just been added as `git submodule` and to build the wrappers you have to fetch all the submodule in this way `git submodule update --init --recursive`



# AUTHOR

HPC4CMB/TOAST Author <br />
Theodore Kisner <work@theodorekisner.com> <br />
Reijo Keskitalo <reijo.keskitalo@gmail.com> <br />
Andrea Zonca <zonca@sdsc.edu> <br />
Giuseppe Puglisi <giuse.puglisi@gmail.com> <br />

CAL Maintainer <br />
Stefano Mandelli <stefano.mandelli@unimi.it>


# License


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

