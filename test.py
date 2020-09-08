from channel_extr_classdef import ChanImage
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import numpy as np

im = ChanImage('Sunderbands_NDVI.tif')
im.preprocess(clahe = False)
im.applyGaborFilter(n = 24, kernel_ra = (2,2), w = 1)
im.extractChannels(value_treshold = 65, size_treshold = 0)
p = PatchCollection(im.patches, match_original = True, alpha = 0.5, edgecolor = 'Darkgreen', facecolor = 'Maroon')

im2 = ChanImage('Sunderbands_NDVI.tif')
im2.applyGaborFilter(n = 24, kernel_ra = (2,2), w = 1)
im2.extractChannels(value_treshold = 65, size_treshold = 0)
p2 = PatchCollection(im2.patches, match_original = True, alpha = 0.5, edgecolor = 'Darkgreen', facecolor = 'Maroon')

fig, ax = plt.subplots(ncols = 4, figsize = (100,25))
ax[0].imshow(im.arr_orig)
ax[0].set_title('original', fontsize = 50)
ax[1].imshow(im.im_prepro)
ax[1].set_title('preprocessed', fontsize = 50)
ax[2].imshow(im.im)
#ax[2].add_collection(p)
ax[2].set_title('preprocessed and Gabor', fontsize = 50)
ax[3].imshow(im2.im)
#ax[3].add_collection(p2)
ax[3].set_title('only Gabor', fontsize = 50)
fig.savefig('test.png')

import os
os.system('nomacs test.png')
