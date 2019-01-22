# Port program written by Dr. Satoru Okuda for curvature analysis to ImageJ (from C/OpenCV)
# Motivation:	let ImageJ do heavy lifting on UI, image processing, and file IO side, leaving code 
#				that does the clever stuff less obscured. Also makes cross-platform deployment easier, 
#				promotes use by biologists and makes looping over multiple inputs easier. 
# D.J. Kelly, 2018-09-26, douglas.kelly@riken.jp

# imports
import os, sys, math
from ij import IJ, ImageStack
from ij.gui import WaitForUserDialog
from ij.plugin import ZProjector, ChannelSplitter, Duplicator
from ij import Prefs
from loci.plugins import BF as bf

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
from InnerOuterComparisonData import InnerOuterComparisonData
from CalculatedObjects import CalculatedObjects

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
	
	# ensure consistent preference settings
	Prefs.blackBackground = False;

	# prompt user for input parameters
	params = mbui.analysis_parameters_gui();

	# prompt user for file locations
	file_path, output_folder = mbio.file_location_chooser(params.input_image_path);
	params.setInputImagePath(file_path);
	params.setOutputPath(output_folder);

	# get image file
	import_opts, params = mbui.choose_series(file_path, params);
	imps = bf.openImagePlus(import_opts);
	imp = imps[0];

	# catch unexpected image dimensions - for now, project in Z...:
	if  imp.getNSlices() > 1:
		mbui.warning_dialog(["More than one Z plane detected.", 
							"I will do a maximum projection before proceeding", 
							"Continue?"]);
		imp = ZProjector.run(imp,"max all");

	params = mbio.get_metadata(params);
	params.setCurvatureLengthUm(round(params.curvature_length_um / params.pixel_physical_size) * params.pixel_physical_size);
	params.persistParameters();
	IJ.run(imp, "Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
	imp.show();
	if imp.getNChannels() > 1:
		imp.setPosition(params.membrane_channel_number, 1, 1);
	mbui.autoset_zoom(imp);
	IJ.run("Enhance Contrast", "saturated=0.35");

	# prompt user to select ROI
	original_imp, crop_params = mbui.crop_to_ROI(imp, params);
	if crop_params is not None:
		params.perform_spatial_crop = True;
		mbui.autoset_zoom(imp);

	# prompt user to do time cropping
	imp, start_end_tuple = mbui.time_crop(imp, params);
	params.setTimeCropStartEnd(start_end_tuple);

	repeats = 1;
	inner_outer_comparisons = None;
	if params.inner_outer_comparison:
		inner_outer_comparisons = InnerOuterComparisonData();
		repeats = 2;

	IJ.selectWindow(imp.getTitle());
	IJ.run("Enhance Contrast", "saturated=0.35");
	for r in range(0, repeats):
		calculated_objects = CalculatedObjects();
		if params.time_crop_start_end[0] is not None:
			calculated_objects.timelist = [idx * params.frame_interval for idx in range(params.time_crop_start_end[0], params.time_crop_start_end[1]+1)]
		else:
			calculated_objects.timelist = [idx * params.frame_interval for idx in range(imp.getNFrames())];
		calculated_objects, params, split_channels = mbfs.generate_edges(imp, params, calculated_objects, (r+1)/repeats);
		membrane_channel_imp = split_channels[0];
		actin_channel_imp = split_channels[1];
		segmentation_channel_imp = split_channels[2];

		calculated_objects = mbfs.calculate_outputs(params,
													calculated_objects, 
													split_channels, 
													repeat_fraction=(r+1)/repeats);
						
		# output colormapped images and kymographs 
		fig_imp_list = mbfs.generate_and_save_figures(imp, calculated_objects, params, membrane_channel_imp, segmentation_channel_imp);
		mbfs.save_csvs(calculated_objects, params);
	
		params.saveParametersToJson(os.path.join(params.output_path, "parameters used.json"));
		imp.changes = False;
		IJ.setTool("zoom");
		if params.close_on_completion or (params.inner_outer_comparison and r!=(repeats-1)):
			for fig_imp in fig_imp_list:
				fig_imp.close();
		elif params.close_on_completion and (r==(repeats-1)):
			imp.close();
	
	if params.inner_outer_comparison:
		tname = "Time, " + params.interval_unit;
		output_folder = os.path.dirname(params.output_path);
		profile = [[((inner, outer), (float(outer)/inner)) for inner, outer in zip(calculated_objects.inner_outer_data.inner_means, calculated_objects.inner_outer_data.outer_means)]];
		mbio.save_profile_as_csv(profile, 
								os.path.join(output_folder, "Intensity ratios.csv"), 
								"outer/inner", 
								xname="inner", 
								yname="outer",
								tname=tname, 
								time_list=calculated_objects.timelist);
		sd_profile = [[((inner, outer), (float(outer)/inner)) for inner, outer in zip(calculated_objects.inner_outer_data.inner_sds, calculated_objects.inner_outer_data.outer_sds)]];
		mbio.save_profile_as_csv(profile, 
								os.path.join(output_folder, "Intensity standard deviation ratios.csv"), 
								"outer sd/inner sd", 
								xname="inner sd", 
								yname="outer sd",
								tname=tname, 
								time_list=calculated_objects.timelist);
	return;
	
# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main()