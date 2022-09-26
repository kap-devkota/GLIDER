from setuptools import setup, find_packages

setup(
    name='dsd-gtools',
    version='1',
    description='DSD, GLIDE and GLIDER tools',
    author='Kapil Devkota',
    author_email='kapil.devkota@tufts.edu',
    url='https://github.com/kap-devkota/GLIDER',
    packages=find_packages(exclude=('tests', 'docs', 'results', 'data'))
)
