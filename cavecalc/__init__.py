"""The Cavecalc module provides a model for cave dripwater and carbonate
chemistry.

Author  :   Robert Owen
email   :   robert.owen@earth.ox.ac.uk

Package structure:
    /__init__.py
    /caves.py           - Core geochemical code & IPhreeqc IO
    /forward_models.py  - High-level scripting of Cavecalc models.
    /setter.py          - Utility code to generate model settings.
    /util.py            - Utility code largely relating to file IO.
    /analyse.py         - Analysis and processing of model output.
    /gui/
        /__init__.py
        /gui.py         - GUI Programming
        /layout.py      - Layout data for GUI
        /mapping.py     - Name mapping for gui.py and setter.py
    /data/
        /__init__.py
        /defaults.py           - default model input parameters
        /phreeqc_templates.py  - PHREEQC format strings for caves.py
        /types_and_limits.py   - Model parameter allowable values for setter.py
        /oxotope.dat           - PHREEQC database file
"""

import cavecalc.caves
import cavecalc.setter
import cavecalc.util
import cavecalc.forward_models
import cavecalc.analyse