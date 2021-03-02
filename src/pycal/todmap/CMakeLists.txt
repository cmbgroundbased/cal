
# Install the python files

install(FILES
    __init__.py
    atm.py
    pointing_math.py
    sim_det_atm.py
    sim_tod.py
    pointing.py
    classes.py
    DESTINATION ${PYTHON_SITE}/pycal/todmap
)
