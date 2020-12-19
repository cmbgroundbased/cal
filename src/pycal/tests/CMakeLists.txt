
# Install the python files

install(FILES
    __init__.py
    _helpers.py
    test_runner.py
    mpi.py
    ops_sim_atm.py
    single_proc_atm.py
    DESTINATION ${PYTHON_SITE}/pycal/tests
)
