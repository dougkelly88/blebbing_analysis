# Port program written by Dr. Satoru Okuda for curvature analysis to ImageJ (from C/OpenCV)
# Motivation:	let ImageJ do heavy lifting on UI, image processing, and file IO side, leaving code 
#				that does the clever stuff less obscured. Also makes cross-platform deployment easier, 
#				promotes use by biologists and makes looping over multiple inputs easier. 
# D.J. Kelly, 2018-09-26, douglas.kelly@riken.jp

# python (jython) imports
import os, sys, math, csv
from datetime import datetime

# java imports - aim to have UI components entirely swing, listeners and layouts awt
from java.awt import Dimension, GridBagLayout, GridBagConstraints, GridLayout
import javax.swing as swing
import javax.swing.table.TableModel

# imagej imports
from ij import IJ, WindowManager, ImagePlus, ImageStack
from ij.gui import Roi, PointRoi, PolygonRoi, GenericDialog, WaitForUserDialog
from ij.io import OpenDialog, DirectoryChooser, FileSaver
from ij.plugin import ChannelSplitter, Straightener
from ij.process import FloatPolygon
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
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
	fixed_anchors_set = set();
	for anchor_idx, anchor in enumerate(anchors_list):
		fixed_anchor = anchor;
		last_dsq = 100000;
		for (x,y) in zip(outline.xpoints,outline.ypoints):
			d2 = math.pow((x - anchor.x), 2) + math.pow((y - anchor.y), 2);
			if d2 < last_dsq:
				last_dsq = d2;
				fixed_anchor = (x, y);
		fixed_anchors_set.add(fixed_anchor);
	if (len(fixed_anchors_set) < (anchor_idx+1)):
		raise ValueError('degeneracy between anchor points!');
	return list(fixed_anchors_set);

# remove all blobs other than the largest by area
def keep_largest_blob(imp):
	rt = ResultsTable();
	mxsz = imp.width * imp.height;
	#pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER | ParticleAnalyzer.DOES_STACKS, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE, rt, 0, mx);
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE, rt, 0, mxsz);

	roim = RoiManager();
	for idx in range(1, imp.getImageStackSize()+1):
		roim.reset();
		rt.reset();
		imp.setPosition(idx);
		pa.analyze(imp);
		#rt_slices = [int(r) for r in rt.getColumn(rt.getColumnIndex("Slice")).tolist()]
		rt_areas = rt.getColumn(rt.getColumnIndex("Area")).tolist();
		#rois = roim.getRoisAsArray().tolist();
		mx_ind = rt_areas.index(max(rt_areas))
		indices_to_remove = [a for a in range(0,len(rt_areas)) if a != mx_ind]
		print(indices_to_remove);
		for rem_idx in indices_to_remove:
			roim.select(imp, rem_idx);
			roim.runCommand(imp, "Fill");
	roim.close();

# return angle between two lines
def angle_between_vecs(u_start, u_end, v_start, v_end):
    u = (u_end[0] - u_start[0], u_end[1] - u_start[1]);
    v = (v_end[0] - v_start[0], v_end[1] - v_start[1]);
    return math.atan2(v[1], v[0]) - math.atan2(u[1], u[0]);

# return vector length
def vector_length(start, end):
	return math.sqrt(math.pow((start[0] - end[0]),2) + math.pow((start[1] - end[1]),2));

# figure out which edge of the roi is the membrane, since IJ might start the roi
# from anywhere along the perimeter w.r.t. the user defined anchors
def get_membrane_edge(roi, fixed_anchors, fixed_midpoint):
	poly = roi.getPolygon();
	started = False;
	e1 = FloatPolygon();
	e2 = FloatPolygon();
	for idx,(x,y) in enumerate(zip(poly.xpoints,poly.ypoints)):
		if (((x, y) == (fixed_anchors[0])) or ((x, y) == fixed_anchors[1])) and not started:
			e2.addPoint(x, y);
			started = True;
		elif (((x, y) == (fixed_anchors[0])) or ((x, y) == fixed_anchors[1])) and started:
			e1.addPoint(x, y);
			started = False;
		if started:
			e1.addPoint(x, y);
		else:
			e2.addPoint(x, y);
	anchors_midpoint = (fixed_anchors[1][0] - fixed_anchors[0][0], 
						fixed_anchors[1][1] - fixed_anchors[0][1]);
	e1_mean = (sum(e1.xpoints)/e1.npoints, sum(e1.ypoints)/e1.npoints);
	e2_mean = (sum(e2.xpoints)/e2.npoints, sum(e2.ypoints)/e2.npoints);
	theta_e1 = angle_between_vecs(fixed_anchors[0], fixed_anchors[1], fixed_anchors[0], e1_mean);
	theta_e2 = angle_between_vecs(fixed_anchors[0], fixed_anchors[1], fixed_anchors[0], e2_mean);
	sign = lambda x: (1, -1)[x < 0]
	if sign(theta_e1) is not sign(theta_e2):
		theta_midpoint = angle_between_vecs(fixed_anchors[0], fixed_anchors[1], fixed_anchors[0], fixed_midpoint);
		use_edge = (e1, e2)[sign(theta_midpoint) == sign(theta_e2)];
	else:
		use_edge = (e1, e2)[vector_length(anchors_midpoint, e1_mean) < vector_length(anchors_midpoint, e2_mean)]
	return 	PolygonRoi(use_edge, Roi.POLYLINE);
	
# return a line profile taking the maximum value over n pixels perpendicular to roi line
def maximum_line_profile(imp, roi, pixel_width):
	imp.setRoi(roi);
	ip = Straightener().straightenLine(imp, pixel_width);
	width = ip.getWidth();
	height = ip.getHeight();
	max_profile = [];
	for x in range(0, width):
		pix = ip.getLine(x, 0, x, height);
		max_profile.append(max(pix));
	#print("Length of maximum profile = " + str(len(max_profile)));
 # debug
	return max_profile;

# generate arrays of points along the membrane that are separated by path length l - currently in pixels
def generate_l_spaced_points(roi, l): 
	poly = roi.getPolygon();
	p1, p2, p3 = ([] for i in range(3));
	for idx,(x,y) in enumerate(zip(poly.xpoints,poly.ypoints)):
		if ((idx > 0) and (idx < poly.npoints - 1)): # ignore first and last points: by definition these won't have anything useful on either side
			# look backwards and calculate pathlength at successive points
			db = 0;
			iidx = idx - 1;
			while ((iidx >= 0) and (db < l)):
				dbnew = db + math.sqrt(math.pow((poly.xpoints[iidx] - poly.xpoints[iidx+1]), 2)  + 
									math.pow((poly.ypoints[iidx] - poly.ypoints[iidx+1]), 2));
				if (dbnew >= l):
					xx = poly.xpoints[iidx+1] + ((l - db)/(dbnew - db))*(poly.xpoints[iidx] - poly.xpoints[iidx+1]);
					yy = poly.ypoints[iidx+1] - ((l - db)/(dbnew - db))*(poly.ypoints[iidx] - poly.ypoints[iidx+1]);
					dbnew = db + math.sqrt(math.pow(((l - db)/(dbnew - db))*(poly.xpoints[iidx] - poly.xpoints[iidx+1]), 2)  + 
									math.pow(((l - db)/(dbnew - db))*(poly.ypoints[iidx] - poly.ypoints[iidx+1]), 2));
				else:
					iidx = iidx-1;
				db = dbnew;
			if (db == l):
				pp1 = (xx, yy);
				pp2 = (x, y);
				# then look forwards ONLY IF backwards search was successful...
				iidx = idx + 1;
				df = 0;
				while ((iidx < poly.npoints - 1) and (df < l)):
					dfnew = df + math.sqrt(math.pow((poly.xpoints[iidx] - poly.xpoints[iidx-1]), 2)  + 
									math.pow((poly.ypoints[iidx] - poly.ypoints[iidx-1]), 2));
					if (dfnew >= l):
						xx = poly.xpoints[iidx-1] + ((l - df)/(dfnew - df))*(poly.xpoints[iidx] - poly.xpoints[iidx-1]);
						yy = poly.ypoints[iidx-1] + ((l - df)/(dfnew - df))*(poly.ypoints[iidx] - poly.ypoints[iidx-1]);
						dfnew = df + math.sqrt(math.pow(((l - df)/(dfnew - df))*(poly.xpoints[iidx] - poly.xpoints[iidx-1]), 2)  + 
										math.pow(((l - df)/(dfnew - df))*(poly.ypoints[iidx] - poly.ypoints[iidx-1]), 2));
					else:
						iidx = iidx+1;
					df = dfnew;
				if (df == l):
					p1.append(pp1);
					p2.append(pp2);
					p3.append((xx, yy));
	return p1, p2, p3;

# generate a line profile of local curvatures using three-point method and SSS theorem
# (see http://mathworld.wolfram.com/SSSTheorem.html)
def calculate_curvature_profile(centre_curv_points, curv_points1, curv_points2, remove_negative_curvatures):
	curvature_profile = [];
	for (cp, p1, p2) in zip(centre_curv_points, curv_points1, curv_points2):
		a = math.sqrt( math.pow((cp[0] - p1[0]), 2) + math.pow((cp[1] - p1[1]), 2) );
		b = math.sqrt( math.pow((cp[0] - p2[0]), 2) + math.pow((cp[1] - p2[1]), 2) );
		c = math.sqrt( math.pow((p2[0] - p1[0]), 2) + math.pow((p2[1] - p1[1]), 2) );
		s = 0.5 * (a + b + c);
		K = math.sqrt(s * (s - a) * (s - b) * (s - c));
		if (K == 0):
			curv = 0;
		else:
			R = (a * b * c)/(4 * K);
			c_1 = tuple(r1-rc for r1, rc in zip(p1, cp));
			c_2 = tuple(r2-rc for r2, rc in zip(p2, cp));
			sign = round((c_1[0]*c_2[1] - c_2[0]*c_1[1]) / (abs(c_1[0]*c_2[1] - c_2[0]*c_1[1]) + 1e-10));
			curv = sign * 1/R;
		if (remove_negative_curvatures and (curv < 0)):
			curv = 0;
		curvature_profile.append((cp, curv));
	return curvature_profile;

# overlay curvature pixels on membrane image
def overlay_curvatures(imp, curvature_stack, curvature_profiles, membrane_channel, limits):
	overlay_base_imp = imp.clone();
	overlay_imp = ImagePlus("Curvature stack", curvature_stack);
	IJ.run(overlay_imp, "Orange Hot", "");
	IJ.setMinAndMax(overlay_imp, int(limits[0]), int(limits[1]));
	IJ.run(overlay_imp, "RGB Color", "");
	overlaid_stack = ImageStack(overlay_imp.width, overlay_imp.height);
	for fridx in range(1, curvature_stack.getSize()+1):
		raw_idx = overlay_base_imp.getStackIndex(membrane_channel, 1, fridx);
		ip = overlay_base_imp.getStack().getProcessor(raw_idx).convertToRGB();
		pix = overlay_imp.getStack().getProcessor(fridx).getPixels();
		base_pix = ip.getPixels();
		for ((x, y), c) in curvature_profiles[fridx-1]:
			if (c > 0):
				base_pix[y * imp.width + x] = pix[y * imp.width + x];
		overlaid_stack.addSlice(ip);
	return ImagePlus("Overlaid curvatures", overlaid_stack);

# generate overlays to display curvature
def generate_curvature_overlays(curvature_profile, curvature_stack): 
	w = curvature_stack.getWidth();
	overlay = IJ.createImage("Curvature", "16-bit", w, curvature_stack.getHeight(), 1);
	pix = overlay.getProcessor().getPixels();
	mx = max([c for ((x, y), c) in curvature_profile]);
	for ((x, y), c) in curvature_profile:
		pix[y * w + x] = int(1000 * c);
	curvature_stack.addSlice(overlay.getProcessor());
	return curvature_stack;

# save curvature
def save_curvature_as_csv(curvature_profiles, file_path):
	f = open(file_path, 'wb');
	writer = csv.writer(f);
	writer.writerow(['Frame', 'x', 'y', 'curvature']);
	for idx, curvature_profile in enumerate(curvature_profiles):
		for ((x, y), c) in curvature_profile:
			writer.writerow([idx, x, y, c]);
	f.close();

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
	timestamp = datetime.strftime(datetime.now(), '%Y-%m-%d %H-%M-%S');
	output_folder = os.path.join(output_root, (timestamp + ' output'));
	os.mkdir(output_folder);
	
	# handle zooming to reasonable size, assuming reasonable screen resolution...
	h = imp.height;
	w = imp.width;
	d = imp.getNSlices();
	n_channels = imp.getNChannels();
	n_frames = imp.getNFrames();
	zoom_factor = 100* math.floor(1080 / (2 * h));
	IJ.run("Set... ", "zoom=" + str(zoom_factor) + " x=" + str(math.floor(w/2)) + " y=" + str(math.floor(h/2)));
	IJ.run("Scale to Fit", "");

	## prompt user to select ROI
	#IJ.setTool("rect");
	#WaitForUserDialog("Crop", "If desired, select a rectangular ROI to crop to...").show();
	#roi = imp.getRoi();
	#if roi is not None:
	#	IJ.run(imp, "Crop", "");
	#	IJ.run("Set... ", "zoom=" + str(zoom_factor) + " x=" + str(math.floor(w/2)) + " y=" + str(math.floor(h/2)));
	#	h = imp.height;
	#	w = imp.width;

	# binarise/segment
	anchors = prompt_for_points(imp, 
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
	keep_largest_blob(membrane_channel_imp);
	
	# generate edge - this needs to be looped over slices
	curvature_stack = ImageStack(w, h);
	actin_profiles = [];
	curvature_profiles = [];
	curv_limits = (10E4, 0);
	for fridx in range(0, n_frames):
		imp.setPosition(membrane_channel, 1, fridx+1);
		membrane_channel_imp.setPosition(fridx+1);
		# 	fix anchors:
		IJ.run(membrane_channel_imp, "Create Selection", "");
		roi = membrane_channel_imp.getRoi();
		fixed_anchors = fix_anchors_to_membrane(anchors, roi);
		fixed_midpoint = fix_anchors_to_membrane(midpoint, roi);
		fixed_midpoint = (midpoint[0].x, midpoint[0].y);

		#	identify which side of the segmented roi to use and perform interpolation/smoothing:
		membrane_edge = get_membrane_edge(roi, fixed_anchors, fixed_midpoint);
		imp.setRoi(membrane_edge);
		imp.show();
		IJ.run(imp, "Interpolate", "interval=0.5 smooth");
		IJ.run(imp, "Fit Spline", "");
		membrane_edge = imp.getRoi();
		
		# generate curvature - this needs to be looped over slices
		curv_points1, centre_curv_points, curv_points2 = generate_l_spaced_points(membrane_edge, 5.0);
		remove_negative_curvatures = True;
		curvature_profiles.append(calculate_curvature_profile(centre_curv_points, 
														curv_points1, 
														curv_points2, 
														remove_negative_curvatures));
		curv_limits = (min(curv_limits[0], min([c[1] for c in curvature_profiles[-1]])), 
						max(curv_limits[1], max([c[1] for c in curvature_profiles[-1]])));
		curvature_stack = generate_curvature_overlays(curvature_profiles[-1], curvature_stack)
		
		# generate actin-channel line profile - assume 2-channel image...
		actin_channel = (membrane_channel + 1) % n_channels;
		actin_channel_imp = split_channels[actin_channel-1];
		actin_profiles.append(maximum_line_profile(actin_channel_imp, membrane_edge, 3));
	
	# output colormapped images and kymographs 
	curvature_stack_imp = overlay_curvatures(imp, curvature_stack, curvature_profiles, membrane_channel, curv_limits);	
	curvature_stack_imp.show();
	FileSaver(curvature_stack_imp).saveAsTiffStack(os.path.join(output_folder, "curvature_stack.tif"));
	save_curvature_as_csv(curvature_profiles, os.path.join(output_folder, "curvatures.csv"))
	FileSaver(membrane_channel_imp).saveAsTiffStack(os.path.join(output_folder, "binary_membrane_stack.tif"));
	IJ.setTool("zoom");

# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main()