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
from ij.gui import Roi, PointRoi, PolygonRoi, NonBlockingGenericDialog, WaitForUserDialog
from ij.io import OpenDialog, DirectoryChooser
from ij.plugin import ChannelSplitter
from ij.process import FloatPolygon
from loci.plugins import BF as bf


# prompt the user to provide a number of points on the image
def prompt_for_points(imp, title, message, n_points):
	imp.killRoi();
	if (n_points == 1): 
		IJ.setTool("point");
	else:
		IJ.setTool("multipoint");
	selected_points = 0;
	while (selected_points != n_points):
		WaitForUserDialog(title, message).show();
		roi = imp.getRoi();
		if roi is not None:
			selected_points = len(roi.getContainedPoints());
		if ((roi is None) or (selected_points != n_points)):
			WaitForUserDialog("Error!", "Wrong number of points selected! Please try again...").show();
	return roi.getContainedPoints();

# move user-defined anchor points onto automatically-segmented membrane
# can do this more efficiently but probably not worth optimising...
# use set to implicitly avoid degeneracy, convert back to list. again, 
# inefficient, but OK for small set
def fix_anchors_to_membrane(anchors_list, membrane_roi):
	outline = membrane_roi.getPolygon();
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
	if (len(fixed_anchors_set) < (anchor_idx+1)):
		raise ValueError('degeneracy between anchor points!');
	return list(fixed_anchors_set);

# figure out which edge of the roi is the membrane, since IJ might start the roi
# from anywhere along the perimeter w.r.t. the user defined anchors
def get_membrane_edge(roi, fixed_anchors, fixed_midpoint):
	smth_polygon = roi.getPolygon();
	xs = smth_polygon.xpoints;
	ys = smth_polygon.ypoints;
	started = False;
	e1 = FloatPolygon();
	e2 = FloatPolygon();
	for idx,(x,y) in enumerate(zip(xs,ys)):
		if (((x, y) == (fixed_anchors[0])) or ((x, y) == fixed_anchors[1])) and not started:
			started = True;
		elif (((x, y) == (fixed_anchors[0])) or ((x, y) == fixed_anchors[1])) and started:
			started = False;
		if started:
			e1.addPoint(x, y);
		else:
			e2.addPoint(x, y);
	if e1.contains(fixed_midpoint[0][0], fixed_midpoint[0][1]):
		return PolygonRoi(e1, Roi.POLYLINE);
	else:
		return PolygonRoi(e2, Roi.POLYLINE);

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
	#timestamp = datetime.strftime(datetime.now(), '%Y-%m-%d %H-%M-%S');
	#output_folder = os.path.join(output_root, (timestamp + ' output'));
	#os.mkdir(output_folder);
	
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
	anchor = prompt_for_points(imp, 
						"Select channel + extrema", 
						"Select the membrane-label channel, and position \n" + 
						"exactly TWO points at extremes of membrane", 
						2);

	midpoint = prompt_for_points(imp, 
								"Choose midpoint", 
								"Now select a point halfway between the extremes, along the membrane", 
								1);
	
	membrane_channel = imp.getChannel();
	split_channels = ChannelSplitter.split(imp);
	membrane_channel_imp = split_channels[membrane_channel-1];
	#membrane_channel_imp.show(); # debug
	
	# perform binary manipulations
	IJ.run(membrane_channel_imp, "Make Binary", "method=Moments background=Dark calculate");
	IJ.run(membrane_channel_imp, "Fill Holes", "stack");
	IJ.run(membrane_channel_imp, "Open", "stack");
	IJ.run(membrane_channel_imp, "Close-", "stack");
	
	# generate edge
	# 	fix anchors:
	IJ.run(membrane_channel_imp, "Create Selection", "");
	roi = membrane_channel_imp.getRoi();
	anchors = fix_anchors_to_membrane(anchor, roi);
	midpoint = fix_anchors_to_membrane(midpoint, roi);

	#	identify which side of the segmented roi to use and perform interpolation/smoothing:
	membrane_edge = get_membrane_edge(roi, anchors, midpoint);
	imp.setRoi(membrane_edge);
	imp.show();
	IJ.run(imp, "Interpolate", "interval=0.5 smooth");
	IJ.run(imp, "Fit Spline", "");
	
	# generate curvature
	# output colormapped images and kymographs 

# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main()