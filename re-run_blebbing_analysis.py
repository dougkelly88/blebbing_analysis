# Code to re-run membrane blebbing analysis on previous data without going through manual edge drawing again
# Motivation:	calculation methods have changed slightly as version number has increased, leading to 
#				difficulties in comparing analysed datasets. To avoid having to redraw membranes, it's useful
#				to have a tool to re-run the analysis, loading previously determined edges and skipping 
#				edge-determination steps
# D.J. Kelly, 2019-01-04, douglas.kelly@riken.jp

# imports
import os, sys, math
from ij import IJ, ImageStack
from ij.io import FileSaver
from ij.gui import WaitForUserDialog
from ij.plugin import ChannelSplitter, Duplicator, ZProjector
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
from Parameters import Parameters

def main():
	# ensure consistent preference settings
	Prefs.blackBackground = False;
	
	# get locations for previous and new outputs
	params = Parameters();
	output_folder_old, output_folder = mbio.rerun_location_chooser(params.input_image_path);
	params.loadParametersFromJson(os.path.join(output_folder_old, 'parameters used.json'));
	params.setOutputPath(output_folder);

	# present user with the familiar setup UI, with options that don't make sense/aren't yet implemented disabled
	params = mbui.analysis_parameters_gui(rerun_analysis=True, params=params);

	# get original image file (for overlays etc.) 
	if not os.path.isfile(params.input_image_path):
		mbui.warning_dialog(["The original data can't be found at the location specified in saved parameters. ", 
							"(Possibly something as simple as a change in network drive mapping has occurred)",
							"Please specify the location of the original image file..."]);
		params.setInputImagePath(mbio.input_file_location_chooser(params.output_path));

	import_opts, params = mbui.choose_series(params.input_image_path, params);
	imps = bf.openImagePlus(import_opts);
	imp = imps[0];
	h = imp.height;
	w = imp.width;
	d = imp.getNSlices();
	n_channels = imp.getNChannels();
	n_frames = imp.getNFrames();
	n_slices = imp.getNSlices();

	if n_slices > 1:
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

	original_imp, crop_params = mbui.crop_to_ROI(imp, params);
	if crop_params is not None:
		params.perform_spatial_crop = True;
		mbui.autoset_zoom(imp);

	timelist = [idx * params.frame_interval for idx in range(n_frames)];
	tname = "Time, " + params.interval_unit
	imp, start_end_tuple = mbui.time_crop(imp, params);
	params.setTimeCropStartEnd(start_end_tuple);
	if start_end_tuple[0] is not None:
		timelist = [idx * params.frame_interval for idx in range(start_end_tuple[0], start_end_tuple[1]+1)];
	h = imp.height;
	w = imp.width;
	d = imp.getNSlices();
	n_channels = imp.getNChannels();
	n_frames = imp.getNFrames();
	n_slices = imp.getNSlices();
	
	split_channels = ChannelSplitter.split(imp);
	membrane_channel_imp = split_channels[params.membrane_channel_number-1];
	if n_channels >= 2:
		actin_channel = (params.membrane_channel_number + 1) % n_channels;
		actin_channel_imp = split_channels[actin_channel-1];
		t0_actin_mean = None;

	# import edges
	membrane_edges = mbio.load_qcd_edges(os.path.join(output_folder_old, "user_defined_edges.json"));
	mbio.save_qcd_edges(membrane_edges, params.output_path);
	
	# do calculations independent of source of edges
	curvature_stack = ImageStack(w, h);
	actin_profiles = [];
	curvature_profiles = [];
	curv_limits = (10E4, 0);
	lengths_areas_and_arearois = [mb.bleb_area(medge) for medge in membrane_edges];
	mb_lengths = [laaroi[0] * params.pixel_physical_size for laaroi in lengths_areas_and_arearois]
	mb_areas =[laaroi[1] * math.pow(params.pixel_physical_size,2) for laaroi in lengths_areas_and_arearois];
	mb_area_rois =[laaroi[2] for laaroi in lengths_areas_and_arearois];
	
	for fridx in range(0, n_frames):
		# generate curvature - this needs to be looped over slices
		membrane_edge = membrane_edges[fridx];
		curv_points = mb.generate_l_spaced_points(membrane_edge, int(round(params.curvature_length_um / params.pixel_physical_size)));
		curv_profile = mb.calculate_curvature_profile(curv_points,
														membrane_edge, 
														params.filter_negative_curvatures);
		curvature_profiles.append(curv_profile);
		curv_limits = (min(curv_limits[0], min([c[1] for c in curvature_profiles[-1]])), 
						max(curv_limits[1], max([c[1] for c in curvature_profiles[-1]])));
		curvature_stack = mbfig.generate_curvature_overlays(curvature_profiles[-1], curvature_stack)
		
		# generate actin-channel line profile if actin channel present
		if n_channels >= 2:
			actin_channel_imp.setPosition(fridx+1);
			actin_channel_imp, t0_actin_mean = mb.apply_photobleach_correction_framewise(params, 
																						actin_channel_imp, 
																						membrane_channel_imp, 
																						t0_value=t0_actin_mean);
			actin_profile = mb.maximum_line_profile(actin_channel_imp, 
													membrane_edge, 
													int(round(params.intensity_profile_width_um / params.pixel_physical_size)));
			actin_profiles.append(actin_profile);
			if params.inner_outer_comparison:
				mean_intensity = float(sum([a for ((x, y), a) in actin_profile]))/len(actin_profile);
				outer_intensity_mean.append(mean_intensity) if r==(repeats-1) else inner_intensity_mean.append(mean_intensity);
				sd = math.sqrt(sum((x - mean_intensity)**2 for ((x, y), a) in actin_profile) / len(actin_profile));
				outer_intensity_sd.append(sd) if r==(repeats-1) else inner_intensity_sd.append(sd);
					
	# output colormapped images and kymographs 
	norm_curv_kym = mbfig.generate_kymograph(curvature_profiles, params.curvature_kymograph_lut_string, "Curvature kymograph - distal point at middle");
	curv_kym = mbfig.generate_plain_kymograph(curvature_profiles, params.curvature_kymograph_lut_string, "Curvature kymograph");
	FileSaver(norm_curv_kym).saveAsTiff(os.path.join(params.output_path, "normalised position curvature kymograph.tif"));
	FileSaver(curv_kym).saveAsTiff(os.path.join(params.output_path, "raw curvature kymograph.tif"));
	overlaid_curvature_imp, raw_curvature_imp = mbfig.overlay_curvatures(imp, curvature_stack, curvature_profiles, params.membrane_channel_number, params, annotate=True);
	mbio.save_profile_as_csv(curvature_profiles, os.path.join(params.output_path, "curvatures.csv"), "curvature")
	bleb_len_imp, bleb_ls = mbfig.plot_bleb_evolution([t * params.frame_interval for t in range(0, len(membrane_edges))], 
											mb_lengths, 
											"Edge length (" + params.pixel_unit + ")");
	bleb_a_imp, bleb_as = mbfig.plot_bleb_evolution([t * params.frame_interval for t in range(0, len(membrane_edges))], 
											mb_areas, 
											"Bleb area (" + params.pixel_unit + "^2)");
	FileSaver(bleb_len_imp).saveAsTiff(os.path.join(params.output_path, "bleb perimeter length.tif"));
	FileSaver(bleb_a_imp).saveAsTiff(os.path.join(params.output_path, "bleb area.tif"));
	mbio.save_1d_profile_as_csv(bleb_ls, os.path.join(params.output_path, "bleb perimeter length.csv"), [("Time, " + params.interval_unit), "Length, " + params.pixel_unit]);
	mbio.save_1d_profile_as_csv(bleb_as, os.path.join(params.output_path, "bleb area.csv"), [("Time, " + params.interval_unit), "Area, " + params.pixel_unit + "^2"]);
	
	# actin channel
	if n_channels >= 2:
		norm_actin_kym = mbfig.generate_kymograph(actin_profiles, params.actin_kymograph_lut_string, (params.labeled_species + " intensity - distal point at middle"));
		actin_kym = mbfig.generate_plain_kymograph(actin_profiles, params.actin_kymograph_lut_string, (params.labeled_species + " intensity"));
		FileSaver(norm_actin_kym).saveAsTiff(os.path.join(params.output_path, "normalised position " + params.labeled_species + " kymograph.tif"));
		FileSaver(actin_kym).saveAsTiff(os.path.join(params.output_path, "raw " + params.labeled_species + " kymograph.tif"));
		mbio.save_profile_as_csv(actin_profiles, 
								os.path.join(params.output_path, 
								(params.labeled_species + " intensities.csv")), 
								(params.labeled_species + " intensity"), 
								tname=tname, 
								time_list=timelist)
		mrg_imp = mbfig.merge_kymographs(norm_actin_kym, norm_curv_kym, params);
		FileSaver(mrg_imp).saveAsTiff(os.path.join(params.output_path, "merged intensity and curvature kymograph.tif"));

	params.saveParametersToJson(os.path.join(params.output_path, "parameters used.json"));
	imp.changes = False;
	IJ.setTool("zoom");
	if params.close_on_completion or (params.inner_outer_comparison and r!=(repeats-1)):
		bleb_a_imp.close();
		bleb_len_imp.close();
		raw_curvature_imp.close();
		overlaid_curvature_imp.close();
		norm_curv_kym.close();
		curv_kym.close();
		if n_channels >= 2:
			mrg_imp.close();
			norm_actin_kym.close();
			actin_kym.close();
	elif params.close_on_completion and (r==(repeats-1)):
		imp.close();
	return;

# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main()
