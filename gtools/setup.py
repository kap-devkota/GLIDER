from setuptools import setup, find_packages

setup(
    name='glidetools',
    version='0.0.1',
    description='DSD, GLIDE and GLIDER tools',
    author='Kapil Devkota',
    author_email='kapil.devkota@tufts.edu',
    url='https://github.com/kap-devkota/GLIDER',
    packages=find_packages(),
    install_requires=["pandas", "numpy", "scipy", "networkx"],
    entry_points = {
        'console_scripts': [
            'glide-compute = commands.glide_main:main'
        ]}
)
