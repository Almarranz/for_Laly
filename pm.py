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
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
import matplotlib.colors as colors_plt
from skimage import data
from skimage import transform
import astroalign as aa
from astroquery.gaia import Gaia
import skimage as ski
import sys
from astropy.stats import sigma_clip


# Pixel scale = 0.0135"/pixel



field = 2

folder_ls = '/Users/amartinez/Desktop/for_people/for_Laly/fields/field_%s/'%(field)
l_a = 'IB_calib_field2_2011.txt'
l_b = 'IB_calib_field2_2013.txt'

l1 = Table.read(folder_ls +l_a, format = 'ascii') 
l2 = Table.read(folder_ls +l_b, format = 'ascii') 

l1_xy = np.array([l1['x'],l1['y']]).T
l2_xy = np.array([l2['x'],l2['y']]).T


d_m = 1
deg = 1
sig_cl = 2



# comp_t = Table(comp, names =('l1_x','l1_y','l2_x','l2_y','ind_1','ind_2','dist'))

for loop in range(2):
    
    comp = compare_lists(l1_xy,l2_xy,d_m)
    print(f'Common in loop {loop} = %s'%(len(comp['ind_1'])))
    l1_com = l1[comp['ind_1']]
    l2_com = l2[comp['ind_2']]
    
    diff_mag = l1_com['IB224_diff'] - l2_com['IB224_diff'] 
    diff_mag1 = l1_com['IB230_diff'] - l2_com['IB230_diff'] 
    
    
    mask_m, l_lim,h_lim = sigma_clip(diff_mag, sigma=sig_cl, masked = True, return_bounds= True)
    
    
    
    l2_clip = l2_com
    l1_clip = l1_com
    
# =============================================================================
#     l1_clip = l1_com[~mask_m.mask]
#     l2_clip = l2_com[~mask_m.mask]
#     
#     fig, (ax,ax1) = plt.subplots(1,2)
#     ax.set_xlabel('$\Delta$ IB224')
#     ax1.set_xlabel('$\Delta$ IB230')
#     
#     ax.hist(diff_mag, label = 'matches = %s\ndist = %.2f arcsec'%(len(comp['ind_1']), d_m*0.0135))
#     ax1.hist(diff_mag1)
#     
#     ax.axvline(l_lim, color = 'red', ls = 'dashed', label = '$\pm$%s$\sigma$'%(sig_cl))
#     ax.axvline(h_lim, color = 'red', ls = 'dashed')
#     
#     ax.legend()
# =============================================================================
    # ax1.legend()
    
    
    
    
    xy_1c = np.array([l1_clip['x'],l1_clip['y']]).T
    xy_2c = np.array([l2_clip['x'],l2_clip['y']]).T
    p = ski.transform.estimate_transform('polynomial',
                                          xy_1c, 
                                          xy_2c, order = deg)
    # p = ski.transform.estimate_transform('similarity',
    #                                       xy_1c, 
    #                                       xy_2c)
    
    l1_xy = p(l1_xy)
    
  



