#!/usr/bin/env python

from setuptools import setup

requires = ['numpy', 'netCDF4']

setup(
    name='wrfncxnj',
    platforms=['GNU/Linux'],
    version='0.1',
    description="WRF netCDF Extract and Join",
    install_requires=requires,
    packages=['wrfncxnj', ],
    include_package_data=True,
    license='GNU General Public License',
    scripts=['bin/wrfncxnj']
)
