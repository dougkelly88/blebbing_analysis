# Port program written by Dr. Satoru Okuda for curvature analysis to ImageJ (from C/OpenCV)
# Motivation:	let ImageJ do heavy lifting on UI, image processing, and file IO side, leaving code 
#				that does the clever stuff less obscured. Also makes cross-platform deployment easier, 
#				promotes use by biologists and makes looping over multiple inputs easier. 
# D.J. Kelly, 2018-09-26, douglas.kelly@riken.jp

# python (jython) imports
import os, sys, math
from datetime import datetime

# java imports - aim to have UI components entirely swing, listeners and layouts awt
from java.awt import Dimension, GridBagLayout, GridBagConstraints, GridLayout
import javax.swing as swing
import javax.swing.table.TableModel

# imagej imports
from ij import IJ, WindowManager
from ij.gui import NonBlockingGenericDialog, WaitForUserDialog
from ij.io import OpenDialog, DirectoryChooser
from ij.plugin import ChannelSplitter
from loci.plugins import BF as bf

# breakout file chooser UI to enable faster debug
def file_location_chooser(default_directory):
	# input
	od = OpenDialog('Choose original file...', 
				 default_directory, 
				 '*.tif');
	file_path = od.getPath();

	# output
	DirectoryChooser.setDefaultDirectory(default_directory);
	dc = DirectoryChooser('Select the root folder for saving output');
	output_root = dc.getDirectory();
	print(output_root);
	if output_root is None:
		print("no output root");
	return file_path, output_root;

# move user-defined anchor points onto automatically-segmented membrane
# can do this more efficiently but probably not worth optimising...
# use set to implicitly avoid degeneracy
def fix_anchors_to_membrane(anchors_list, membrane_roi):
	outline = membrane_roi.getContainedFloatPoints();
	xs = outline.xpoints;
	ys = outline.ypoints;
	fixed_anchors_set = set();
	for anchor_idx, anchor in enumerate(anchors_list):
		# debug...
		#print("manual anchor " + str(anchor_idx) + ", unfixed:");
		#print("(" + str(anchor.x) + ", " + str(anchor.y) + ")");
		last_dsq = 100000;
		for idx in xrange(1,  outline.npoints):

			d2 = math.pow((xs[idx] - anchor.x), 2) + math.pow((ys[idx] - anchor.y), 2);
			if d2 < last_dsq:
				last_dsq = d2;
				fixed_anchor = (xs[idx], ys[idx]);	
		fixed_anchors_set.add(fixed_anchor);
		# debug...
		#print("fixed anchor " + str(anchor_idx) + ":");
		#print("(" + str(fixed_anchor[0]) + ", " + str(fixed_anchor[1]) + ")");
	if (len(fixed_anchors_set) < 3) :
		raise ValueError('degeneracy between anchor points!');
	return fixed_anchors_set;

def main():
	#print (sys.version_info) # debug
	#print(sys.path) # debug
	file_path = "D:\\data\\Inverse blebbing\\MAX_2dpf marcksl1b-EGFP inj_TgLifeact-mCh_movie e4_split-bleb1.tif" # debug
	output_root = "D:\\data\\Inverse blebbing\\output" # debug

	## prompt user for file locations
	#default_directory = "D:\\data\\Inverse blebbing\\";
	#file_path, output_root = file_location_chooser(default_directory);
	#if ((file_path is None) or (output_root is None)): 
	#	return;

	# get image file
	imps = bf.openImagePlus(file_path);
	imp = imps[0];
	imp.show();
	stack = imp.getStack();
	
	## prepare output folders
	timestamp = datetime.strftime(datetime.now(), '%Y-%m-%d %H-%M-%S');
	output_folder = os.path.join(output_root, (timestamp + ' output'));
	os.mkdir(output_folder);
	
	
	# handle zooming to reasonable size, assuming reasonable screen resolution...
	h = imp.height;
	w = imp.width;
	d = imp.getNSlices();
	c = imp.getNChannels();
	f = imp.getNFrames();
	zoom_factor = 100* math.floor(1080 / (2 * h));
	IJ.run("Set... ", "zoom=" + str(zoom_factor) + " x=" + str(math.floor(w/2)) + " y=" + str(math.floor(h/2)));
	IJ.run("Scale to Fit", "");

	# prompt user to select ROI
	IJ.setTool("rect");
	WaitForUserDialog("Crop", "If desired, select a rectangular ROI to crop to...").show();
	roi = imp.getRoi();
	if roi is not None:
		IJ.run(imp, "Crop", "");
		IJ.run("Set... ", "zoom=" + str(zoom_factor) + " x=" + str(math.floor(w/2)) + " y=" + str(math.floor(h/2)));
		h = imp.height;
		w = imp.width;

	# binarise/segment
	IJ.setTool("multipoint");
	selected_points = 0;
	while (selected_points != 3):
		WaitForUserDialog("Select channel", "Select the membrane-label channel, and position exactly three points at extremes of membrane and in the middle...").show();
		roi = imp.getRoi();
		if roi is not None:
			selected_points = len(roi.getContainedPoints());
		if ((roi is None) or (selected_points != 3)):
			WaitForUserDialog("Error!", "Wrong number of points selected! Please try again...").show();

	anchor1 = roi.getContainedPoints()[0];
	anchor2 = roi.getContainedPoints()[1];
	anchor = roi.getContainedPoints();
	
	membrane_channel = imp.getChannel();
	split_channels = ChannelSplitter.split(imp);
	membrane_channel_imp = split_channels[membrane_channel-1];
	membrane_channel_imp.show();
	# perform binary manipulations
	IJ.run(membrane_channel_imp, "Make Binary", "method=Moments background=Dark calculate");
	IJ.run(membrane_channel_imp, "Fill Holes", "stack");
	IJ.run(membrane_channel_imp, "Open", "stack");
	IJ.run(membrane_channel_imp, "Close-", "stack");
	
	# generate edge
	# 	fix anchors:
	IJ.run(membrane_channel_imp, "Create Selection", "stack");
	roi = membrane_channel_imp.getRoi();
	anchors = fix_anchors_to_membrane(anchor, roi);
	print(anchors);
	
	
	# generate curvature
	# output colormapped images and kymographs 

# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main()