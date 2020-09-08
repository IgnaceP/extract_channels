""" Script to define a class to load an image and extract channels
    largely based on Yang et al. 2015
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from shapely.geometry import Polygon as shPol
from osgeo import gdal
import cv2

class ChanImage:
    def __init__(self, path):
        self.path = path
        self.geo, self.bands, self.arr = self.importRaster(path)
        self.rows, self.cols = np.shape(self.arr)

        self.arr_orig = self.arr.copy()

        self.im = np.array(255*(self.arr + 1)/2, dtype = np.uint8)

    @staticmethod
    def importRaster(file, band = 1):
        """
        Method to load a geotiff and translate to gdal geo object, a band and its array

        Args:
            - band: (int, defaults to 1) which band to extract
        """
        Raster = gdal.Open(file)
        Band= Raster.GetRasterBand(band)
        Array= Band.ReadAsArray()
        return(Raster, Band, Array)

    def preprocess(self, movingmean = True, clahe = True, sharpen = False, bilateral = True):
        """
        Method to apply the following pre-processing steps onto the image:
            - moving mean
        """

        if movingmean:
            arr_mean = self.movingmean(3)
            self.im = np.array(255*(arr_mean + 1)/2, dtype = np.uint8)
        else:
            self.im = np.array(255*(self.arr + 1)/2, dtype = np.uint8)

        if clahe:
            self.applyCLAHE()

        if sharpen:
            self.applySharpening()

        if bilateral:
            self.im = cv2.bilateralFilter(self.im, 5 ,25, 25)

        self.im_prepro = self.im.copy()

    def movingmean(self, kernel_size = 3):
        """
        Method to calculate a moving mean with a square kernel

        Args:
            - kernel_size: (int, defaults to 3) total width or height of square kernel
        """

        arr = self.arr
        arr_cp = arr.copy()

        # initiate the mean array
        arr_mean = np.zeros([self.rows-(kernel_size-1), self.cols-(kernel_size-1)])

        # nested loop to loop over all the windowed arrays
        for i in range(0,kernel_size):
            for j in range(0, kernel_size):
                arr_mean += arr[i:i+self.rows-(kernel_size-1),j:j+self.cols-(kernel_size-1)]
        # the total of windowed arrays = (i+1)/(j+1)
        arr_mean /= (i+1)*(j+1)

        return arr_mean

    def applyCLAHE(self):
        """
        Method to apply a Contrast limited Adaptive Histogram Equalization
        """

        # create a CLAHE object and apply it
        clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize = (8,8))
        self.arr_clahe = clahe.apply(self.im)

    def applySharpening(self, smoothing_param = 10):

        #self.im = cv2.addWeighted(self.im, cv2.GaussianBlur(self.im, (smoothing_param, smoothing_param)))
        pass

    def applyGaborFilter(self, w = 2, kernel_ra = (2,2), n = 12):
        """
        Method to apply GABOR filter to enhance the curvilinear characteristics of channels
        """

        angles = np.linspace(-np.pi/2, np.pi/2, n)
        gabor_resp = np.zeros([self.rows, self.cols, n])

        for i in range(n):

            theta = angles[i]

            # bounding box
            x,y = np.meshgrid(np.arange(-kernel_ra[0], kernel_ra[0]+1), np.arange(-kernel_ra[1], kernel_ra[1]+1))

            # define parameters based on desired detectable channel width
            W = 2*w + 1
            sigmax = sigmay = W/(2*(2*np.log(2))**0.5)
            f0 = 1/W

            # rotation of coordinates
            x_ = x*np.cos(theta) + y*np.sin(theta)
            y_ = y*np.cos(theta) + x*np.sin(theta)

            # real Gabor filter
            gb = (1/(2*sigmax*sigmay))*np.exp(-0.5*(x_**2+y_**2)/(sigmax*sigmay))*np.cos(2*np.pi*f0*x_)

            try: gabor_resp[1:-1,1:-1,i] = cv2.filter2D(self.im,-1,gb)
            except: gabor_resp[:,:,i] = cv2.filter2D(self.im,-1,gb)

        gb = np.max(gabor_resp, axis = 2)
        self.im = gb

    def extractChannels(self, value_treshold = 0, size_treshold = 0):

        # make a binary image
        (binary_treshold, binary) = cv2.threshold(src = self.im,
                                                  thresh = value_treshold,
                                                  maxval = 255,
                                                  type = cv2.THRESH_BINARY)

        # convert to proper data type
        binary = np.uint8(binary)

        # contour
        contours = cv2.findContours(image = binary, mode = cv2.RETR_LIST, method = cv2.CHAIN_APPROX_SIMPLE)
        contours = contours[0]
        contours_len = []

        self.patches = []

        # loop over all found contours
        for C in contours:
            if len(C)> size_treshold:
                # adapt the dimensions
                C = C[:,0,:]

                # create a polygon
                pol = Polygon(C, closed = True, facecolor = 'red', edgecolor = 'red', alpha = 0.05) # plt polygon
                self.patches.append(pol)
