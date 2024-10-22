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

import Polywarp as pw


# Pixel scale = 0.0135"/pixel



field = 2

folder_ls = '/Users/amartinez/Desktop/for_people/for_Laly/fields/field_%s/'%(field)
pruebas = '/Users/amartinez/Desktop/for_people/for_Laly/pruebas/'

l_a = 'IB_calib_field2_2011.txt'
l_b = 'IB_calib_field2_2013.txt'
#%%
l1 = Table.read(folder_ls +l_a, format = 'ascii') 
l2 = Table.read(folder_ls +l_b, format = 'ascii') 
sig_cl_ls = np.arange(1,3,0.1)

# with open(pruebas + 'sig_and_com.txt', 'w') as file:
#     file.write('# sig_cl num of mathches\n')


    

sig_cl =3
d_m = 1
deg = 1
print(f'Degree = {deg} sig = {sig_cl}')
for loop in range(10):
    l1_xy = np.array([l1['x'],l1['y']]).T
    l2_xy = np.array([l2['x'],l2['y']]).T
    comp = compare_lists(l1_xy,l2_xy,d_m)
    print(f'Common in loop {loop} = %s'%(len(comp['ind_1'])))
    # if loop == 1:
    #     with open(pruebas + 'sig_and_com.txt', 'a') as file:
    #         file.write('%.1f %.0f\n'%(sig_cl, len(comp['ind_1'])))

    l1_com = l1[comp['ind_1']]
    l2_com = l2[comp['ind_2']]
    
    diff_mag = l1_com['IB224_diff'] - l2_com['IB224_diff'] 
    diff_mag1 = l1_com['IB230_diff'] - l2_com['IB230_diff'] 
    diff_x =  l2_com['x'] - l1_com['x'] 
    diff_y =  l2_com['y'] - l1_com['y'] 
    diff_xy = (diff_x**2 + diff_y**2)**0.5
    mask_m, l_lim,h_lim = sigma_clip(diff_mag, sigma=sig_cl, masked = True, return_bounds= True)
    
    
    
    # l2_clip = l2_com
    # l1_clip = l1_com
    
    l1_clip = l1_com[~mask_m.mask]
    l2_clip = l2_com[~mask_m.mask]
    
    fig, (ax,ax1) = plt.subplots(1,2)
    ax.set_xlabel('$\Delta$ IB224')
    ax.hist(diff_mag, label = 'matches = %s\ndist = %.2f arcsec'%(len(comp['ind_1']), d_m*0.0135))
    ax.axvline(l_lim, color = 'red', ls = 'dashed', label = '$\pm$%s$\sigma$'%(sig_cl))
    ax.axvline(h_lim, color = 'red', ls = 'dashed')
    ax.legend()
    
    ax1.hist(diff_x, label = 'x', histtype = 'step')
    ax1.hist(diff_y, label = 'y',histtype = 'step')
    ax1.hist(diff_xy, label = 'xy',color = 'grey', alpha = 0.3)
    
    # ax1.hist(comp['dist'], color = 'pink', alpha = 0.5)
    ax1.hist(comp['dist'], color = 'red',  histtype = 'step', lw=2)
    ax1.set_xlabel('$\Delta$ pixel')
    # ax1.set_xlabel('$\Delta$ IB230')
    # ax1.hist(diff_mag1)
    ax1.legend()
    
    
    
    
    xy_1c = np.array([l1_clip['x'],l1_clip['y']]).T
    xy_2c = np.array([l2_clip['x'],l2_clip['y']]).T
    
    # Kx,Ky=pw.polywarp(xy_1c[:,0],xy_1c[:,1],xy_2c[:,0],xy_2c[:,1],degree=1)
    Kx,Ky=pw.polywarp(xy_2c[:,0],xy_2c[:,1],xy_1c[:,0],xy_1c[:,1],degree=deg)
    
    xi=np.zeros(len(l1))
    yi=np.zeros(len(l1))
    
    for k in range(deg+1):
                for m in range(deg+1):
                    xi=xi+Kx[k,m]*l1['x']**k*l1['y']**m
                    yi=yi+Ky[k,m]*l1['x']**k*l1['y']**m
    
    l1['x'] = xi
    l1['y'] = yi
    
    
    # p = ski.transform.estimate_transform('polynomial',
    #                                       xy_1c, 
    #                                       xy_2c, order = deg)
    # print(p.params)
    # sys.exit()
    # p = ski.transform.estimate_transform('affine',
    #                                       xy_1c, 
    #                                       xy_2c)
    
    # l1_xy = p(l1_xy)
    
    # d_m += 0.51



