from setuptools import setup

setup(
    name='WIMpy',
    version='1.2',
    description='Calculates the DM direct detection signal according to the NREFT',
    author='Bradley Kavanagh',
    author_email='kavanagh@ifca.unican.es',
    url = 'https://github.com/bradkav/WIMpy_NREFT',
    packages=['WIMpy'],
    include_package_data=True,
    package_data={'WIMpy': ['Nuclei.txt', 'nu_spectra/*'], }
    #data_files=['nuclei.txt'],
    #data_files=[('WIMpy',['nuclei.txt'])]
)
