from PIL import Image
import numpy as np
import os
import xml.etree.ElementTree as ET

class ReadImageCube():
    def __init__(self,image_dir,band_list,rotate_angle):
        """
        :param image_dir: directory with tif image of palimpsest
        :param band_list: list of bands
        :param coord: (left, upper, right, lower) tuple of bounding box coordinates
        """
        self.image_dir = image_dir
        if len(band_list)==0:
            for fname in os.listdir(self.image_dir):
                if ".tif" in fname:
                    band_list.append(fname[:-4])
        self.band_list = sorted(band_list)
        self.rotate_angle = rotate_angle
        self.max_values = self.max_of_band_im(25)
        self.stretch_contrast = True

    #open file as multipage tif
    #close file
    def thumbnail(self,band_idx,scale_ratio):
        """Opens image thumbnail"""
        path = os.path.join(self.image_dir,self.band_list[band_idx]+".tif")
        im = self.read_image(path,False,scale_ratio,[0,0,0,0],self.max_values[band_idx])
        return im

    def max_of_band_im(self,scale_ratio):
        """Reads image of multiple bands"""
        self.stretch_contrast = False
        max_vals = np.zeros((len(self.band_list)))
        for idx, band_name in enumerate(self.band_list):
            path = os.path.join(self.image_dir, self.band_list[idx] + ".tif")
            im = self.read_image(path, False, scale_ratio, [0, 0, 0, 0],1)
            max_vals[idx] = find_max(im,np.amax(im))
        return max_vals

    def read_image(self,path,crop,scale_ratio,coord,max_val):
        """Reads image of one band
        Args:
            coords - [left, upper, right, lower]
            """

        with Image.open(path) as im:
            if crop:
                "coords - left, upper, right, lower"
                im = im.crop(coord)
            if scale_ratio>1:
                width, height = im.size
                new_size = (width//scale_ratio,height//scale_ratio)
                im = im.resize(new_size)
            if self.rotate_angle>0:
                rotation = eval("Image.ROTATE_{}".format(self.rotate_angle))
                im = im.transpose(rotation)

            im = np.array(im)
            if self.stretch_contrast:
                im = self.strech_contrast_fun(im,max_val)

        return im

    def read_msi_image(self,coord):
        """Reads image of multiple bands
        coord - [x1,y1,x2,y2]"""
        size = (coord[3]-coord[1],coord[2]-coord[0])
        ims = np.zeros(((len(self.band_list),)+size))
        for idx, band_name in enumerate(self.band_list):
            fpath = os.path.join(self.image_dir,band_name+".tif")
            im = self.read_image(fpath,True,1,coord,self.max_values[idx])
            ims[idx]= im
        return ims

    def strech_contrast_fun(self,im,max_val):
        """Strech im by max value without oversaturated pixels
        max_val - int, bit depth"""
        im = np.clip(im,a_max=max_val,a_min=0.)
        im = im/max_val
        return im


def find_max(im,max_val):
    """Find max value without oversaturated pixels
    max_val - scalar, bit depth"""
    if max_val == 1:
        bin_width = 1 / 256.0
    elif max_val == 256:
        bin_width = 1
    else:
        bin_width = 10
    bins = np.arange(0, max_val, bin_width)
    hist, bins = np.histogram(im, bins=bins)
    hist = hist/np.sum(hist)
    for idx in reversed(range(len(hist))):
        if hist[idx] < 0.001 and hist[idx - 1] == 0:
        # print("Hist now {}, hist before {}".format(hist[-idx],hist[-idx-1]))
            max_val = bins[idx - 1]
        else:
            break
    return max_val



