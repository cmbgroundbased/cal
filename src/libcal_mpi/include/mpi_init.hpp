
// Copyright (c) 2015-2020 by the parties listed in the AUTHORS file.
// All rights reserved.  Use of this source code is governed by
// a BSD-style license that can be found in the LICENSE file.

#ifndef CAL_MPI_HPP
#define CAL_MPI_HPP

#include <mpi.h>
#include <CAL_MPI_AtmSim.hpp>

namespace cal
{

void mpi_init(int argc, char * argv[]);

void mpi_finalize();

}

#endif // ifndef CAL_MPI_HPP
