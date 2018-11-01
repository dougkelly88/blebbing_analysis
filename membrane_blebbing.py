# Port program written by Dr. Satoru Okuda for curvature analysis to ImageJ (from C/OpenCV)
# Motivation:	let ImageJ do heavy lifting on UI, image processing, and file IO side, leaving code 
#				that does the clever stuff less obscured. Also makes cross-platform deployment easier, 
#				promotes use by biologists and makes looping over multiple inputs easier. 
# D.J. Kelly, 2018-09-26, douglas.kelly@riken.jp

# imports
import os, sys, math
from ij import IJ, ImageStack
from ij.io import FileSaver
from ij.gui import WaitForUserDialog
from ij.plugin import ChannelSplitter, Duplicator
from loci.formats import ImageReader
from loci.plugins import BF as bf

script_path = os.path.dirname(os.path.realpath(__file__));
if "Fiji.app" in script_path:
	ss = script_path.split("Fiji.app");
	final_folder = os.path.basename(script_path);
	script_path = os.path.join(ss[0], "Fiji.app", "plugins", "Scripts", "Plugins", final_folder);
sys.path.insert(0, os.path.join(script_path, 'modules'));
sys.path.insert(0, os.path.join(script_path, 'classes'));

import membraneBlebbingFileio as mbio;
import membraneBlebbingUi as mbui;
import membraneBlebbingEngine as mb;
import membraneBlebbingFigures as mbfig;
from Parameters import Parameters

def main():
	### debug ###############################################################
	#print (sys.version_info);
	#print(sys.path);
	#file_path = "D:\\data\\Inverse blebbing\\MAX_2dpf marcksl1b-EGFP inj_TgLifeact-mCh_movie e4_split-bleb1.tif";
	#output_root = "D:\\data\\Inverse blebbing\\output";
	#from datetime import datetime
	#timestamp = datetime.strftime(datetime.now(), '%Y-%m-%d %H-%M-%S');
	#output_folder = os.path.join(output_root, (timestamp + ' output')); 
	#os.mkdir(output_folder); 
	########################################################################
	
	# prompt user for input parameters
	params = mbui.analysis_parameters_gui();

	# prompt user for file locations
	file_path, output_folder = mbio.file_location_chooser(params.input_image_path);
	params.setInputImagePath(file_path);

	# get image file
	imps = bf.openImagePlus(file_path);
	imp = imps[0];
	imp.show();
	reader = ImageReader();
	reader.setId(file_path);
	params.setFrameInterval(reader.getMetadataValue("Frame Interval").value());
	params.setIntervalUnit(reader.getMetadataValue("Frame Interval").unit().getSymbol())
	params.setPixelPhysicalSize(1/reader.getMetadataValue("YResolution"));
	params.setPixelSizeUnit(reader.getMetadataValue("Unit"))
	mbui.autoset_zoom(imp);

	# prompt user to select ROI
	if params.perform_spatial_crop:
		original_imp, crop_params = mbui.crop_to_ROI(imp);
		if crop_params is not None:
			print("cropped");
			params.setSpatialCrop(crop_params.toString());
			mbui.autoset_zoom(imp);

	# prompt user to do time cropping
	if params.perform_time_crop:
		imp, start_end_tuple = mbui.time_crop(imp)
		params.setTimeCropStartEnd(start_end_tuple);

	h = imp.height;
	w = imp.width;
	d = imp.getNSlices();
	n_channels = imp.getNChannels();
	n_frames = imp.getNFrames();

	# binarise/segment
	anchors = mbui.prompt_for_points(imp, 
						"Select channel + extrema", 
						"Select the membrane-label channel, and position \n" + 
						"exactly TWO points at extremes of membrane", 
						2);
	midpoint = mbui.prompt_for_points(imp, 
								"Choose midpoint", 
								"Now select a point halfway between the extremes, along the membrane", 
								1);
	params.setManualAnchorMidpoint(midpoint);
	params.setManualAnchorPositions(anchors);
	
	membrane_channel = imp.getChannel();
	split_channels = ChannelSplitter.split(imp);
	membrane_channel_imp = split_channels[membrane_channel-1];
	
	# perform binary manipulations
	membrane_channel_imp = mb.make_and_clean_binary(membrane_channel_imp, params.threshold_method);
	
	# generate edge - this needs to be looped over slices
	curvature_stack = ImageStack(w, h);
	actin_profiles = [];
	curvature_profiles = [];
	membrane_edges = [];
	curv_limits = (10E4, 0);
	for fridx in range(0, n_frames):
		imp.setPosition(membrane_channel, 1, fridx+1);
		membrane_channel_imp.setPosition(fridx+1);
		# 	fix anchors:
		IJ.run(membrane_channel_imp, "Create Selection", "");
		roi = membrane_channel_imp.getRoi();
		fixed_anchors = mb.fix_anchors_to_membrane(anchors, roi);
		fixed_midpoint = midpoint[0];

		#	identify which side of the segmented roi to use and perform interpolation/smoothing:
		membrane_edge = mb.get_membrane_edge(roi, fixed_anchors, fixed_midpoint);
		imp.setRoi(membrane_edge);
		imp.show();
		IJ.run(imp, "Interpolate", "interval=0.5 smooth");
		IJ.run(imp, "Fit Spline", "");
		membrane_edge = imp.getRoi();
		membrane_edges.append(membrane_edge);
		
		# generate curvature - this needs to be looped over slices
		curv_points = mb.generate_l_spaced_points(membrane_edge, params.curvature_length_pix);
		curvature_profiles.append(mb.calculate_curvature_profile(curv_points,
																membrane_edge, 
																params.filter_negative_curvatures));
		curv_limits = (min(curv_limits[0], min([c[1] for c in curvature_profiles[-1]])), 
						max(curv_limits[1], max([c[1] for c in curvature_profiles[-1]])));
		curvature_stack = mbfig.generate_curvature_overlays(curvature_profiles[-1], curvature_stack)
		
		# generate actin-channel line profile - assume 2-channel image...
		actin_channel = (membrane_channel + 1) % n_channels;
		actin_channel_imp = split_channels[actin_channel-1];
		actin_profiles.append(mb.maximum_line_profile(actin_channel_imp, membrane_edge, 3));
	
	# output colormapped images and kymographs 
	norm_curv_kym = mbfig.generate_kymograph(curvature_profiles, params.curvature_kymograph_lut_string, "Curvature kymograph - distal point at middle");
	curv_kym = mbfig.generate_plain_kymograph(curvature_profiles, params.curvature_kymograph_lut_string, "Curvature kymograph");
	norm_actin_kym = mbfig.generate_kymograph(actin_profiles, params.actin_kymograph_lut_string, (params.labeled_species + " intensity - distal point at middle"));
	actin_kym = mbfig.generate_plain_kymograph(actin_profiles, params.actin_kymograph_lut_string, (params.labeled_species + " intensity"));
	FileSaver(actin_kym).saveAsTiff(os.path.join(output_folder, "normalised position " + params.labeled_species + " kymograph.tif"));
	FileSaver(actin_kym).saveAsTiff(os.path.join(output_folder, "raw " + params.labeled_species + " kymograph.tif"));
	FileSaver(norm_curv_kym).saveAsTiff(os.path.join(output_folder, "normalised position curvature kymograph.tif"));
	FileSaver(curv_kym).saveAsTiff(os.path.join(output_folder, "raw curvature kymograph.tif"));
	overlaid_curvature_imp, raw_curvature_imp = mbfig.overlay_curvatures(imp, curvature_stack, curvature_profiles, membrane_channel, curv_limits, params);	
	overlaid_curvature_imp.show();
	mbui.autoset_zoom(overlaid_curvature_imp);
	FileSaver(overlaid_curvature_imp).saveAsTiffStack(os.path.join(output_folder, "overlaid curvature.tif"));
	FileSaver(raw_curvature_imp).saveAsTiffStack(os.path.join(output_folder, "raw curvature.tif"));
	mbio.save_profile_as_csv(curvature_profiles, os.path.join(output_folder, "curvatures.csv"), "curvature")
	mbio.save_profile_as_csv(actin_profiles, os.path.join(output_folder, (params.labeled_species + " intensities.csv")), (params.labeled_species + " intensity"))
	FileSaver(membrane_channel_imp).saveAsTiffStack(os.path.join(output_folder, "binary_membrane_stack.tif"));
	mrg_imp = mbfig.merge_kymographs(norm_actin_kym, norm_curv_kym, params);
	bleb_len_imp, bleb_ls = mbfig.plot_bleb_evolution([t * params.frame_interval for t in range(0, len(membrane_edges))], 
											[mb.roi_length(medge) * params.pixel_physical_size for medge in membrane_edges], 
											"Edge length (" + params.pixel_unit + ")");
	bleb_a_imp, bleb_as = mbfig.plot_bleb_evolution([t * params.frame_interval for t in range(0, len(membrane_edges))], 
											[mb.bleb_area(medge)[0] * math.pow(params.pixel_physical_size,2) for medge in membrane_edges], 
											"Bleb area (" + params.pixel_unit + "^2)");
	FileSaver(bleb_len_imp).saveAsTiff(os.path.join(output_folder, "bleb perimeter length.tif"));
	FileSaver(bleb_a_imp).saveAsTiff(os.path.join(output_folder, "bleb perimeter area.tif"));
	mbio.save_1d_profile_as_csv(bleb_ls, os.path.join(output_folder, "bleb perimeter length.csv"), [("Time, " + params.interval_unit), "Length, " + params.pixel_unit]);
	mbio.save_1d_profile_as_csv(bleb_as, os.path.join(output_folder, "bleb area.csv"), [("Time, " + params.interval_unit), "Area, " + params.pixel_unit + "^2"]);
	FileSaver(mrg_imp).saveAsTiff(os.path.join(output_folder, "merged intensity and curvature kymograph.tif"));
	#mbfig.generate_intensity_weighted_curvature(raw_curvature_imp, curvature_profiles, actin_channel_imp, "physics");
	params.saveParametersToJson(os.path.join(output_folder, "parameters used.json"));
	imp.changes = False;
	IJ.setTool("zoom");
	if params.close_on_completion:
		imp.close();
		mrg_imp.close();
		bleb_a_imp.close();
		bleb_len_imp.close();
		raw_curvature_imp.close();
		overlaid_curvature_imp.close();
		norm_curv_kym.close();
		curv_kym.close();
		norm_actin_kym.close();
		actin_kym.close();


# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main()