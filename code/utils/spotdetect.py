#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Script to identify spot locations and generate QC output.
Python script to identify spot locations using Hough Transform from high-res 
images output by SpaceRanger. Also outputs various QC plots.
Args:
    infolder: path to folder containing outs folder output by SpaceRanger
    -o --outfolder: path to output folder for location csv and QC images
    -s --samplename: sample name used for generating output folders and labeling figure
Output:
    spots_detected.csv: csv of positions of detected spots. Three columns, X position, Y position, radius. No header
    tissue_hires_spotdetect.png: image of spots detected overlaid on tissue_hires_image.png (SpaceRanger output)
    detected_tissue_spotdetect.png: image of spots detected overlaid on detected_tissue_image.jpg (SpaceRanger output)
    high_contrast_spotdetect.png: image of spots detected overlaid on high-res image after contrast manipulation
    high_contrast_centers_spotdetect.png: image of spot centers detected overlaid on high-res image after contrast manipulation
    high_contrast_centers_tenx.png: image of spot centers determined by SpaceRanger overlaid on high-res image after contrast manipulation
    high_contrast_centers_both.png: image of spot centers determined by SpaceRanger (white) and spot centers ID'ed by spotdetect (black) overlaid on high-res image after contrast manipulation
Usage:
spotdetect.py /home/spot_checking/Human_FFPE_PDAC_E2547_ROI4_s1_D -o /home/spot_checking/outfolder -s Human_FFPE_PDAC_E2547_ROI4_s1_D 
@author: Hannah Pliner
Created on Mon Feb 20 16:17:41 2023

"""

import os
import sys
import argparse
import cv2 as cv
import numpy as np
import pandas as pd
import json


def draw_circles(image, circles, color_outer = (0,255,0), color_center = (0,0,255)):
    for i in circles[0,:]:
        cv.circle(image,(i[0],i[1]),i[2],color_outer,1)
        cv.circle(image,(i[0],i[1]),2,color_center,1)
    return image

def main(arguments):

    parser = argparse.ArgumentParser()
    parser.add_argument('infolder', help="Path to input folder")
    parser.add_argument('-o', '--outfolder', help="Path to output folder")
    parser.add_argument('-s', '--samplename', help="Sample ID for output files")
    args = parser.parse_args(arguments)

    # Make output folder if doesn't exist
    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    hires_img = cv.imread(os.path.join(args.infolder,'spatial/tissue_hires_image.png'))

    # Grayscale
    gray = cv.cvtColor(hires_img, cv.COLOR_BGR2GRAY)

    # Equalize histogram to increase contrast
    gray = cv.equalizeHist(gray)

    # ID circles
    circles = cv.HoughCircles(gray,cv.HOUGH_GRADIENT,1,20,
                              param1=200,param2=12,minRadius=5,maxRadius=10) #param2=15
    circles = np.uint16(np.around(circles))

    # Output csv
    np.savetxt(os.path.join(args.outfolder, args.samplename + "_spots_detected.csv"), circles[0,:], delimiter=",")

    # Output images
    hires_img = draw_circles(hires_img, circles)
    cv.imwrite(os.path.join(args.outfolder, args.samplename + "_tissue_hires_spotdetect.png"), hires_img)

    det_img = cv.imread(os.path.join(args.infolder, 'spatial/detected_tissue_image.jpg'))
    det_img = draw_circles(det_img, circles)
    cv.imwrite(os.path.join(args.outfolder, args.samplename + "_detected_tissue_spotdetect.png"), det_img)

    gray_cent = gray
    for i in circles[0,:]:
      cv.drawMarker(gray_cent, (i[0], i[1]), (0,0,0), markerSize=4, thickness=1, line_type=cv.LINE_AA)
      
    cv.imwrite(os.path.join(args.outfolder, args.samplename + "_high_contrast_centers_spotdetect.png"), gray_cent)
    
    gray_c = draw_circles(gray, circles, color_outer = (0,0,255))
    cv.imwrite(os.path.join(args.outfolder, args.samplename + "_high_contrast_spotdetect.png"), gray_c)
    
    hires_img = cv.imread(os.path.join(args.infolder,'spatial/tissue_hires_image.png'))

    # Grayscale
    gray = cv.cvtColor(hires_img, cv.COLOR_BGR2GRAY)

    # Equalize histogram to increase contrast
    gray = cv.equalizeHist(gray)
    
    with open(os.path.join(args.infolder,'spatial/scalefactors_json.json'), 'r') as f:
      scales = json.load(f)

    try:
        tenx_pos = pd.read_csv(os.path.join(args.infolder,'spatial/tissue_positions_list.csv'), header=None)
        df1 = tenx_pos[[4,5]].mul(scales['tissue_hires_scalef']).round(0).astype(int)
    except:
        tenx_pos = pd.read_csv(os.path.join(args.infolder,'spatial/tissue_positions.csv'))
        df1 = tenx_pos[['pxl_row_in_fullres','pxl_col_in_fullres']].mul(scales['tissue_hires_scalef']).round(0).astype(int)
    
    arr = df1.to_numpy()
    
    for item in arr:
      cv.drawMarker(gray, (item[1], item[0]), (255,0,0), markerSize=4, thickness=1, line_type=cv.LINE_AA)
    
    cv.imwrite(os.path.join(args.outfolder, args.samplename + "_high_contrast_centers_tenx.png"), gray)
    
    for i in circles[0,:]:
      cv.drawMarker(gray, (i[0], i[1]), (0,0,0), markerSize=4, thickness=1, line_type=cv.LINE_AA)
      
    cv.imwrite(os.path.join(args.outfolder, args.samplename + "_high_contrast_centers_both.png"), gray)

    return


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

