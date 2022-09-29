from setuptools import setup, find_packages
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory /"../README.md").read_text()
setup(
    name='glidetools',
    version='0.0.4',
    description='DSD, GLIDE and GLIDER tools',
    author='Kapil Devkota',
    license_files = {'LICENSE.txt',},
    long_description=long_description,
    long_description_content_type="text/markdown",
    author_email='kapil.devkota@tufts.edu',
    url='https://github.com/kap-devkota/GLIDER',
    packages=find_packages(),
    install_requires=["pandas", "numpy", "scipy", "networkx"],
    entry_points = {
        'console_scripts': [
            'glide-compute = glidetools.commands.glide_main:main'
        ]}
)
