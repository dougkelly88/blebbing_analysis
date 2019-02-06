from ij import IJ, ImagePlus, Prefs
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable, Measurements
from ij.plugin import ChannelSplitter
from ij.gui import PolygonRoi, WaitForUserDialog

import math
import os
import sys

# ensure consistent preference settings
Prefs.blackBackground = False;

release = False;

if not release:
	script_path = os.path.dirname(os.path.realpath(__file__));
else: 
	script_path = os.getcwd();
if "Fiji.app" in script_path:
	ss = script_path.split("Fiji.app");
	final_folder = "blebbing analysis";
	script_path = os.path.join(ss[0], "Fiji.app", "plugins", "Scripts", "Plugins", final_folder);
sys.path.insert(0, os.path.join(script_path, 'modules'));
sys.path.insert(0, os.path.join(script_path, 'classes'));

import membraneBlebbingFileio as mbio;
import membraneBlebbingUi as mbui;
import membraneBlebbingEngine as mb;
import membraneBlebbingFigures as mbfig;
import membraneBlebbingHighLevelFunctions as mbfs;

from Parameters import Parameters

def get_background_stdev(intensity_imp, mask_imp, params, membrane_edges, dilations=5, threshold_method=None, membrane_imp=None):
	"""get intensity standard deviation from a background region identified by auto-thresholded image and position of midpoint anchor"""
	if mask_imp is None and membrane_imp is not None:
		# do thresholding using either previous method if threhsold_method is None or using (less conservative?) threshold method
		if (threshold_method is None or not (threshold_method in params.listThresholdMethods())):
			mask_imp = mb.make_and_clean_binary(membrane_imp, params.threshold_method);
			print("using thresholding method " + str(params.threshold_method));
		else:
			mask_imp = mb.make_and_clean_binary(membrane_imp, threshold_method);
			print("using thresholding method " + str(threshold_method));
			
	std_devs = [];
	rois = [];
	IJ.setForegroundColor(0, 0, 0);
	roim = RoiManager(True);
	rt = ResultsTable();

	for fridx in range(imp.getNFrames()):
		mask_imp.setT(fridx+1);
		intensity_imp.setT(fridx+1);
		# add extra bit to binary mask from loaded membrane in case user refined edges...
		# flip midpoint anchor across the line joining the two extremes of the membrane, 
		# and fill in the triangle made by this new point and those extremes
		poly = membrane_edges[fridx].getPolygon();
		l1 = (poly.xpoints[0], poly.ypoints[0]);
		l2 = (poly.xpoints[-1], poly.ypoints[-1]);
		M = (0.5*(l1[0] + l2[0]), 0.5*(l1[1] + l2[1]));
		Mp1 = (params.manual_anchor_midpoint[0][0] - M[0], params.manual_anchor_midpoint[0][1] - M[1]);
		p2 = (M[0] - Mp1[0], M[1] - Mp1[1]);
		new_poly_x = list(poly.xpoints);
		new_poly_x.append(p2[0]);
		new_poly_y = list(poly.ypoints);
		new_poly_y.append(p2[1]);
		mask_imp.setRoi(PolygonRoi(new_poly_x, new_poly_y, PolygonRoi.POLYGON));
		IJ.run(mask_imp, "Fill", "slice");
		mask_imp.killRoi();
		
		# now dilate the masked image and identify the unmasked region closest to the midpoint anchor
		ip = mask_imp.getProcessor();
		dilations = 5;
		for d in range(dilations):
			ip.dilate();
		mask_imp.setProcessor(ip);
		IJ.run(mask_imp, "Create Selection", "");
		ip = mask_imp.getProcessor();
		ip.invert();
		mxsz = mask_imp.getWidth() * mask_imp.getHeight();
		pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER | ParticleAnalyzer.SHOW_PROGRESS, ParticleAnalyzer.CENTROID, rt, 0, mxsz);
		pa.setRoiManager(roim);
		pa.analyze(mask_imp);
		ds_to_anchor = [math.sqrt((x - params.manual_anchor_midpoint[0][0])**2 + (y - params.manual_anchor_midpoint[0][1])**2) 
				  for x, y in zip(rt.getColumn(rt.getColumnIndex("X")).tolist(), rt.getColumn(rt.getColumnIndex("Y")).tolist())];
		if len(ds_to_anchor)>0:
			roi = roim.getRoi(ds_to_anchor.index(min(ds_to_anchor)));
			rois.append(roi);
			intensity_imp.setRoi(roi);
			
			# TODO: add user intervention step to check the success of assigning background regions!
	
			# the standard deviation within this region will be the background SD
			std_dev = intensity_imp.getStatistics().stdDev;
			print("original image frame " + str(fridx+1) +  " standard deviation = " + str(std_dev));
			std_devs.append(std_dev);
			intensity_imp.killRoi();
		else:
			std_devs.append(0);
			rois.append(None);
		roim.reset();
		rt.reset();
	roim.close();
	return std_devs, rois;


mask_path = "C:\\Users\\dougk\\Desktop\\output\\2019-01-29 17-10-59 output\\Varying curvature length parameter\\curvature_length_um = 1.0\\binary_membrane_stack.tif";
original_img_path = "D:\\source\\vascular_morphogenesis_ij\\blebbing analysis\\test data\\Simple data.tif";
params_path = "C:\\Users\\dougk\\Desktop\\output\\2019-01-29 17-10-59 output\\Varying curvature length parameter\\curvature_length_um = 1.0\\parameters used.json";
edges_path = "C:\\Users\\dougk\\Desktop\\output\\2019-01-29 17-10-59 output\\Varying curvature length parameter\\user_defined_edges.zip";

params = Parameters();
params.loadParametersFromJson(params_path);
membrane_edges = mbio.load_qcd_edges2(edges_path);
mask_imp = IJ.openImage(mask_path);
mask_imp.show();
mask_imp.killRoi();
imp = IJ.openImage(original_img_path);
split_channels = ChannelSplitter.split(imp);
actin_imp = split_channels[0];
membrane_imp = split_channels[1];
membrane_imp.show();
actin_imp.show();

#print(get_background_stdev(actin_imp, mask_imp, params, membrane_edges));
print(get_background_stdev(actin_imp, None, params, membrane_edges, threshold_method='MinError', membrane_imp=membrane_imp));