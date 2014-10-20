# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

long_desc = '''
Python parts of the NGS pipeline at the Charite CompBio group.
'''

requires = []

setup(
    name='cbpipeline',
    version='0.1.0',
    url='http://github.com/holtgrewe/cbqc',
    download_url='http://github.com/holtgrewe/cbqc',
    license='MIT',
    author='Manuel Holtgrewe',
    author_email='manuel.holtgrewe@charite.de',
    description='Python parts of NGS pipeline at the Charite CompBio group',
    long_description=long_desc,
    zip_safe=False,
    platforms='any',
    packages=find_packages(),
    include_package_data=True,
    install_requires=requires,
    namespace_packages=['cbpipeline'],
    scripts=['qc_report.sh'],
)
