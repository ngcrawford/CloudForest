from distutils.core import setup

setup(
    name='CloudForest',
    version='0.0.1',
    author='Nicholas G. Crawford',
    author_email='ngcrawford@gmail.com',
    packages=['cloudforest', 'cloudforest.test'], # Not tests written yet...
    url='http://pypi.python.org/pypi/CloudForest/',
    license='LICENSE.txt',
    description='Generate gene and species trees from thousands of loci using Map-Reduce.',
    long_description=open('README.rst').read(),
    install_requires=[
        "mrjob >= 0.2.8",
        "flask >= 0.8",
        "flask-wtf >= 0.5.2",
        "numpy >= 1.6.1"
    ],
    package_data = {
            '':['*.txt'],
            'cloudforest':[
                'bin/*',
                'gzips/*',
                'test/*',
                'test/alignments/*',
                'test/alignments/nexus_primates/*',
                'test/alignments/phylip_primates/*',
            ]
        },
    include_package_data = True,
)
