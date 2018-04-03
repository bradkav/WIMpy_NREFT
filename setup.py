from setuptools import setup

setup(
    name='WIMpy',
    version='1.0.0',
    description='Calculates the DM direct detection signal according to the NREFT',
    author='Bradley Kavanagh',
    author_email='b.j.kavanagh@uva.nl',
    url = 'https://github.com/bradkav/WIMpy_NREFT',
    packages=['WIMpy'],
    package_data={'WIMpy': [] },
    long_description=open('README.md').read(),
)
