import distribute_setup
distribute_setup.use_setuptools()
from setuptools import setup

setup(
    name='CloudForest',
    version='0.0.1',
    author='Nicholas G. Crawford',
    author_email='ngcrawford@gmail.com',
    url='http://pypi.python.org/pypi/CloudForest/',
    license='LICENSE.txt',
    description='Generate gene and species trees from thousands of loci using Map-Reduce.',
    long_description=open('README.rst').read(),
    test_suite='tests',
    install_requires=[
        "mrjob >= 0.2.8",
        "flask >= 0.8",
        "flask-wtf >= 0.5.2",
        "numpy >= 1.6.1"
    ],
    packages=[
            'cloudforest',
            'cloudforest.test',
            'cloudforest.webapp',
            ],
    package_data = {
            '':['*.txt'],
            'cloudforest':[
                'test/alignments/*.oneliners',
                'test/alignments/*.phylip',
                'test/alignments/nexus_primates/*',
                'test/alignments/phylip_primates/*',
                'webapp/static/images/*',
                'webapp/static/javascripts/*',
                'webapp/static/stylesheets/*',
                'webapp/templates/*'
            ]
        },
    include_package_data = True,
)
