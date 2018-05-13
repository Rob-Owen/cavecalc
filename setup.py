import os.path
import sys
import subprocess

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# check system meets basic requirements
assert sys.platform in ('darwin', 'linux', 'win32'), \
"Platform not supported: %s" % (sys.platform)
assert sys.version_info >= (3,0), \
"Python version not supported. Python 3.5+ is recommended."

if sys.argv[1] == 'install':
    # update pip if necessary
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--upgrade', 'pip'])
    # install compiled dependencies with pip
    # for some reason explicitly downloading with pip works better / quicker.
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pandas'])
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'numpy'])
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'scipy'])
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'matplotlib'])
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'phreeqpy'])

setup(  name =          'cavecalc',
        description =   'Cave Geochemical Modelling',
        author =        'Robert Owen',
        author_email =  'robert.owen@earth.ox.ac.uk',
		url =			'https://www.earth.ox.ac.uk/people/robert-owen/',
        version =       '1.1',
        packages =      ['cavecalc', 'cavecalc.data', 'cavecalc.gui'],
        package_data =  {'cavecalc.data' : ['*.dat']},
        scripts =       ['scripts/cc_input_gui.py', 'scripts/cc_output_gui.py'],
        install_requires = [ 'scipy', 'numpy', 'matplotlib', 'inflection',
                             'seaborn']   )
        
# install phreeqpy separately for non-Windows users.
# phreeqpy doesn't work properly in install_requires
if sys.argv[1] == 'install' and sys.platform.lower() != 'win32':
    import phreeqpy
    print("Done\nAttempting to patch phreeqpy using 2to3...")
    loc = os.path.dirname(os.path.abspath(phreeqpy.__file__))
    com_file = os.path.join(loc, 'iphreeqc', 'phreeqc_com.py')
    subprocess.check_call(["2to3", "-w", com_file])
else:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pywin32'])

print("Cavecalc installation complete. Run example1.py to test.")