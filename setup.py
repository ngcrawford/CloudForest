from distutils.core import setup

setup(
    name='CloudForest',
    version='0.0.1',
    author='Nicholas G. Crawford',
    author_email='ngcrawford@gmail.com',
    packages=['cloudforest', 'cloudforest.test'], # Not tests written yet...
    scripts=['bin/stowe-towels.py','bin/wash-towels.py'],
    url='http://pypi.python.org/pypi/CloudForest/',
    license='LICENSE.txt',
    description='Generate gene and species trees from thousands of loci using Map-Reduce.',
    long_description=open('README.txt').read(),
    install_requires=[
        "mrjob >= 0.2.8",
        "flask >= 0.8",
        "flask-wtf >= 0.5.2",
        "numpy >= 1.6.1"
    ],
)