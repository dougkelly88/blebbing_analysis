# Code to perform looped analysis of membrane blebbing, changing a single parameter to predefined values to examine effect
# Motivation:	it may not always be obvious what the best parameter values are - in this case, 
#				looping over values and looking at the output might be useful
# D.J. Kelly, 2019-01-08, douglas.kelly@riken.jp

# intended operation: have user prepare txt/csv file in format "{"parameter name (copied from parameters used.json)": 
#					[val1, val2, val3]}"; when looped analysis is run, prompt for the location of this file and save
#					to output folders named for parameter name and value

# imports
import os, sys, json
from ij import IJ
from ij.plugin import ZProjector
from ij import Prefs
from loci.plugins import BF as bf

release = True;

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
from CalculatedObjects import CalculatedObjects

def generate_params_list(loop_file_path, params):
	"""given a text file containing the name and list of values of a parameter to loop, generate a list of Parameter objects"""
	params_list = [];
	params_dict = params.__dict__;
	f = open(loop_file_path, 'r');
	try:
		loop_dct = json.loads(f.read());
	except Exception as e:
		print("Error reading loop definition file.")
		print("Please ensure that the loop file conforms to the following format:");
		print("{\"parameter name\": [value1, value2, value3]}");
		raise e;

	parameter_name = list(loop_dct.keys())[0];
	parameter_vals = list(loop_dct.values())[0];
	for v in parameter_vals:
		params_dict[parameter_name] = v;
		np = Parameters();
		np.populate_parameters_from_dict(params_dict);
		params_list.append(np);
	return params_list, parameter_name, parameter_vals;

def main():
	# ensure consistent preference settings
	Prefs.blackBackground = False;
	
	# get locations for previous and new outputs
	params = Parameters();
	output_folder_old, output_folder = mbio.rerun_location_chooser(params.input_image_path);
	params.loadParametersFromJson(os.path.join(output_folder_old, 'parameters used.json'));
	params.setOutputPath(os.path.join(output_folder, "looped analysis"));
	os.makedirs(params.output_path);

	# prompt user for path to text file describing the looped analysis
	loop_file_path = mbio.input_file_location_chooser(params.output_path, filt='*.txt', message="Please choose a text file containing the looping data...");

	# get original image file (for overlays etc.) 
	if not os.path.isfile(params.input_image_path):
		mbui.warning_dialog(["The original data can't be found at the location specified in saved parameters. ", 
							"(Possibly something as simple as a change in network drive mapping has occurred)",
							"Please specify the location of the original image file..."]);
		params.setInputImagePath(mbio.input_file_location_chooser(params.output_path));

	import_opts, params = mbui.choose_series(params.input_image_path, params);
	imps = bf.openImagePlus(import_opts);
	imp = imps[0];

	if imp.getNSlices() > 1:
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

	# import edges
	membrane_edges = mbio.load_qcd_edges(os.path.join(output_folder_old, "user_defined_edges.json"));
	mbio.save_qcd_edges(membrane_edges, params.output_path);
	
	# generate list of Parameters objects to loop over
	params_list, loop_param_name, loop_param_vals = generate_params_list(loop_file_path, params);
	dummy_output_path = params.output_path;

	for params, loop_param_val in zip(params_list, loop_param_vals):
		params.setOutputPath(os.path.join(dummy_output_path, (loop_param_name + " = " + str(loop_param_val))));
		params.close_on_completion = True;
		params.setCurvatureLengthUm(round(params.curvature_length_um / params.pixel_physical_size) * params.pixel_physical_size);
		os.makedirs(params.output_path);
		calculated_objects = CalculatedObjects();
		calculated_objects.membrane_edges = membrane_edges;
		if params.time_crop_start_end[0] is not None:
			calculated_objects.timelist = [idx * params.frame_interval for idx in range(params.time_crop_start_end[0], params.time_crop_start_end[1]+1)]
		else:
			calculated_objects.timelist = [idx * params.frame_interval for idx in range(imp.getNFrames())];

		split_channels = mbfs.split_image_plus(imp, params)
		membrane_channel_imp = split_channels[0];
		actin_channel_imp = split_channels[1];
		segmentation_channel_imp = split_channels[-1];
		if params.photobleaching_correction:
			segmentation_channel_imp = mb.make_and_clean_binary(segmentation_channel_imp, params.threshold_method)

		calculated_objects = mbfs.calculate_outputs(params,
													calculated_objects, 
													split_channels);
						
		# output colormapped images and kymographs 
		fig_imp_list = mbfs.generate_and_save_figures(imp, calculated_objects, params, membrane_channel_imp, segmentation_channel_imp);
		mbfs.save_csvs(calculated_objects, params);
	
		params.saveParametersToJson(os.path.join(params.output_path, "parameters used.json"));
		imp.changes = False;
		IJ.setTool("zoom");
		if params.close_on_completion:
			for fig_imp in fig_imp_list:
				fig_imp.close();
	if params.close_on_completion:
		imp.close();
		
	return;

# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main()
