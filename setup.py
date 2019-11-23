from setuptools import setup, find_packages

with open('README.txt', 'r') as fh: 
    long_description = fh.read()

# noinspection PyPackageRequirements
setup(
    name = 'microcircuit',
    version = '0.0.0',
    packages = find_packages(),
    author = 'HJJ',
    author_email = 'h.jiang@fz-juelich.de',
    description = 'microcircuit package',
    long_description = long_description,
    license = '',
    url = 'git@github.com:HanjiaJiang/microcircuit.git',    
)

