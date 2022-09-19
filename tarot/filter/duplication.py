#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 17:54:08 2018

@author: kanthanakorn
"""
import numpy as np

#Input idx (referenct to catalog's idx) but output is in data's idx
def Duplicate_Match(idx):
    same = {}; keep_dup_data = []; data_dup = [] # keep_same = {};

    for i in idx:
        same[i] = same.get(i, 0) + 1
        
    print(same)
    
    for i, j in same.items():
        if j >= 2:
            keep_dup_data.extend(np.where(idx == i)) #provide data index

    for i in keep_dup_data:
        data_dup.extend(i)
        
    print(data_dup)
    
    Duplicate_idx = sorted(data_dup)
    return Duplicate_idx

# =============================================================================
# Test duplication in python3
# =============================================================================
if False:
    num = np.random.randint(5, size=(1,10))
    print(num)
    
    a = Duplicate_Match(num.flatten())
    print(len(a))
    print(a)
else:
    pass