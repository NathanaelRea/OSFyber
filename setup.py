import sys
from setuptools import setup

if sys.version_info[0] == 3 and sys.version_info[1] < 5:
    sys.exit('Sorry, Python < 3.5 is not supported')

setup(
    name='OSFyber',
    version='0.1.0',
    description='2D Moment Curvature',
    author='Nathanael Rea',
    author_email='nrea@eng.ucsd.edu',
    download_url='https://github.com/NathanaelRea/OSFyber',
    license='GPL-3.0',
    packages=['osfyber'],
    include_package_data=True,
    install_requires=[
        'matplotlib',
        'numpy',
        'pybind11',
        'meshpy'
    ],
)