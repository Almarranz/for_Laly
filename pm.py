#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 17:45:11 2024

@author: amartinez
"""
import numpy as np
import matplotlib.pyplot as plt
from compare_lists import compare_lists 
from astropy.table import Table

# Pixel scale = 0.0135"/pixel

folder_ls = '/Users/amartinez/Desktop/for_people/for_Laly/field_6/'
l_a = 'IB224_calib_field6_2011.txt'
l_b = 'IB224_calib_field6_2013.txt'

l1 = Table.read(folder_ls +l_a, format = 'ascii') 
l2 = Table.read(folder_ls +l_b, format = 'ascii') 

l1_xy = np.array([l1['x'],l1['y']]).T
l2_xy = np.array([l2['x'],l2['y']]).T

comp = compare_lists(l1_xy,l2_xy,1)

l1_com = l1[comp[:,4].astype(int)]
l2_com = l2[comp[:,5].astype(int)]
diff_mag = l1_com['IB224_diff'] - l2_com['IB224_diff'] 
fig, ax = plt.subplots(1,1)
ax.hist(diff_mag)
