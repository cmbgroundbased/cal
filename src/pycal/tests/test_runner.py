from ..mpi import MPI, use_mpi

import os
import sys
import unittest

from .mpi import MPITestRunner
from ..vis import set_backend

# if atm_available: # where is this variable?
#     from . import ops_sim_atm as testopsatm

from . import ops_sim_atm as testopsatm
from . import single_proc_atm as testsingleatm


def test(name=None, verbosity=2):
    # We run tests with COMM_WORLD if available
    comm = None
    rank = 0
    if use_mpi:
        comm = MPI.COMM_WORLD
        rank = comm.rank

    set_backend()

    outdir = "cal_test_output"

    if rank == 0:
        outdir = os.path.abspath(outdir)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)

    if comm is not None:
        outdir = comm.bcast(outdir, root=0)

    # Run python tests.

    loader = unittest.TestLoader()
    mpirunner = MPITestRunner(comm, verbosity=verbosity, warnings="ignore")
    suite = unittest.TestSuite()

    # suite.addTest(loader.loadTestsFromModule(testopsatm))
    suite.addTest(loader.loadTestsFromModule(testsingleatm))

    ret = 0
    _ret = mpirunner.run(suite)
    
    if not _ret.wasSuccessful():
        ret += 1

    if ret > 0:
        sys.exit(ret)

    return ret
