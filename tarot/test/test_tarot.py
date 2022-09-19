#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 16:40:17 2018

@author: kanthanakorn
"""

from unittest import TestCase

try:
    basestring
except NameError:
    basestring = str

import tarot
class TestTarot(TestCase):
    def test_is_string(self):
        stg = tarot.TAROT_PIP.MAIN()
        self.assertTrue(isinstance(stg, basestring))