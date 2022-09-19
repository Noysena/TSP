#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Jun 25 15:02:09 2019

@author: K. NOYSENA
"""

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='Transient search source',
	version = '0.2',
	description = 'Transient search algorithm',
	long_description = 'This search algortihm is designed to locate unknown candidate in a very large field of view image in TAROT telescopes.',
	classifiers = [
		'Development Status :: 1 - Plaining',
		'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
         'Natural Language :: English',
         'Operating System :: POSIX :: Linux'
		'Programming Language :: Python :: 3.7+',
		'Topic :: Utilities',
		],
	keywords = 'Transient search algorithim',
	url = 'https://github.com/Noysena/TSP',
	author = 'K. Noysena',
	author_email = 'noysena@astronu.com',
	license = 'GPL-3.0',
	packages = ['tarot'],
	install_requires=[
            'scipy',
            'numpy',
            'beautifulsoup4',
            'h5py',
            'pyyaml',
            'matplotlib',
            'pytz',
            'scikit-image',
            'pandas',
            'objgraph',
            'setuptools',
            'mock',
            'astropy',
            'astroquery',
            'healpy',
            'cryptography',
            'scikit-learn',
            'numba',
	    'MOCPy',
            'pyds9',
            'lalsuite',
	    'joblib',
            'ligo.skymap',],
	include_package_data=True,
	zip_safe = False,
	scripts = ['bin/transient-search'],
	test_suite = 'nose.collector',
	tests_require = ['nose'],
)
