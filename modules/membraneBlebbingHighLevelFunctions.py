# intermediate layer between main function and lower-level functions, 
# used for historical reasons to avoid a complete refactor
#
# D. J. Kelly, 2019-01-08, douglas.kelly@riken.jp

# imports
import os, math
from ij import IJ
from ij.io import FileSaver
from ij.plugin import ChannelSplitter, Duplicator

import membraneBlebbingEngine as mb
import membraneBlebbingFigures as mbfig
import membraneBlebbingFileio as mbio
import membraneBlebbingUi as mbui

def split_image_plus(imp, params):
	"""split original ImagePlus by channel, and assign an image to segment on"""
	split_channels = ChannelSplitter.split(imp);
	membrane_channel_imp = split_channels[params.membrane_channel_number-1];
	segmentation_channel_imp = Duplicator().run(membrane_channel_imp);
	if params.use_single_channel:
		actin_channel = params.membrane_channel_number;
		actin_channel_imp = Duplicator().run(membrane_channel_imp);
	else:
		if imp.getNChannels() >= 2:
			actin_channel = (params.membrane_channel_number + 1) % imp.getNChannels();
			actin_channel_imp = split_channels[actin_channel-1];
		else:
			actin_channel = params.membrane_channel_number;
			actin_channel_imp = Duplicator().run(membrane_channel_imp);
	split_channels = [membrane_channel_imp, actin_channel_imp, segmentation_channel_imp];
	return split_channels;

def generate_edges(imp, params, calculated_objects, repeat_fraction=1):
	"""generate edges automatically, then clean up with user step if necessary"""
	# binarise/segment
	extra_prompt_str = "";
	if params.inner_outer_comparison:
		if repeat_fraction == 1:
			extra_prompt_str = " OUTER";
			params.setOutputPath(os.path.join(os.path.dirname(params.output_path), "outer"));
			os.mkdir(params.output_path);
		else:
			extra_prompt_str = " INNER";
			params.setOutputPath(os.path.join(params.output_path, "inner"));
			os.mkdir(params.output_path);
	anchors = mbui.prompt_for_points(imp, 
						"Select channel + extrema", 
						"Select the membrane-label channel, and position \n" + 
						"exactly TWO points at extremes of" + extra_prompt_str + " membrane", 
						2);
	midpoint = mbui.prompt_for_points(imp, 
								"Choose midpoint", 
								"Now select a point halfway between the extremes, distant from \n the membrane in the direction of bleb formation. ", 
								1);
	membrane_channel = imp.getChannel();
	params.setMembraneChannelNumber(membrane_channel);
	params.persistParameters();
	params.setManualAnchorMidpoint(midpoint);
	anchors = mb.order_anchors(anchors, midpoint);
	params.setManualAnchorPositions(anchors);
	split_channels = split_image_plus(imp, params);
	membrane_test_channel_imp = Duplicator().run(split_channels[0]);
	segmentation_channel_imp = split_channels[-1];
	
	# perform binary manipulations
	segmentation_channel_imp = mb.make_and_clean_binary(segmentation_channel_imp, params.threshold_method);
		
	# generate edge - this needs to be looped over slices
	membrane_edges = [];
	alternate_edges = [];
	fixed_anchors_list = [];
	previous_anchors = [];
	for fridx in range(0, imp.getNFrames()):
		imp.setPosition(membrane_channel, 1, fridx+1);
		segmentation_channel_imp.setPosition(fridx+1);
		# 	fix anchors:
		IJ.run(segmentation_channel_imp, "Create Selection", "");
		roi = segmentation_channel_imp.getRoi();
		fixed_anchors = mb.fix_anchors_to_membrane(anchors, roi, params);
		fixed_midpoint = midpoint[0];
		# evolve anchors...
		if not params.inner_outer_comparison and not params.constrain_anchors:
			previous_anchors, anchors = mb.evolve_anchors(previous_anchors, fixed_anchors);
		# identify which side of the segmented roi to use and perform interpolation/smoothing:
		membrane_edge, alternate_edge = mb.get_membrane_edge(roi, fixed_anchors, fixed_midpoint);
		fixed_anchors = mb.order_anchors(fixed_anchors, midpoint);
		fixed_anchors_list.append(fixed_anchors);
		membrane_edge = mb.check_edge_order(fixed_anchors, membrane_edge);
		alternate_edge = mb.check_edge_order(fixed_anchors, alternate_edge);
		imp.show();
		imp.setRoi(membrane_edge);
		IJ.run(imp, "Interpolate", "interval=1.0 smooth adjust");
		IJ.run(imp, "Fit Spline", "");
		#membrane_edge = mb.selectionInterpolateAndFitSpline(membrane_edge);
		membrane_edge = imp.getRoi();
		imp.setRoi(alternate_edge);
		IJ.run(imp, "Interpolate", "interval=1.0 smooth adjust");
		IJ.run(imp, "Fit Spline", "");
		alternate_edge = imp.getRoi();
		#alternate_edge = mb.selectionInterpolateAndFitSpline(alternate_edge);
		membrane_edges.append(membrane_edge);
		alternate_edges.append(alternate_edge);
	
	# perform user QC before saving anything
	if params.perform_user_qc:
		imp.hide();
		if imp.getNFrames() > 1:
			membrane_edges, fixed_anchors_list = mbui.perform_user_qc(membrane_test_channel_imp, membrane_edges, alternate_edges, fixed_anchors_list, params);
		else:
			membrane_edges, fixed_anchors_list = mbui.perform_user_qc(imp, membrane_edges, alternate_edges, fixed_anchors_list, params);
		imp.show();
	else:
		mbio.save_qcd_edges2(membrane_edges, params.output_path);	

	calculated_objects.membrane_edges = membrane_edges;
	calculated_objects.fixed_anchors_list = fixed_anchors_list;

	return calculated_objects, params, split_channels;

def calculate_outputs(params, calculated_objects, split_channels, inner_outer_intensity_data=None, repeat_fraction=1):
	"""generate curvatures, lengths, areas etc."""
	membrane_channel_imp = split_channels[0];
	actin_channel_imp = split_channels[1];
	segmentation_channel_imp = split_channels[2];
	# do calculations independent of source of edges
	actin_profiles = [];
	curvature_profiles = [];
	lengths_areas_and_arearois = [mb.bleb_area(medge, params.manual_anchor_midpoint[0]) for medge in calculated_objects.membrane_edges];
	mb_lengths = [laaroi[0] * params.pixel_physical_size for laaroi in lengths_areas_and_arearois]
	full_membrane_lengths = [params.pixel_physical_size * edge.getLength() for edge in calculated_objects.membrane_edges];
	euclidean_membrane_lengths = [params.pixel_physical_size * mb.vector_length((edge.getPolygon().xpoints[0], edge.getPolygon().ypoints[0]), 
																			 (edge.getPolygon().xpoints[-1], edge.getPolygon().ypoints[-1])) for edge in calculated_objects.membrane_edges];
	mb_areas =[laaroi[1] * math.pow(params.pixel_physical_size,2) for laaroi in lengths_areas_and_arearois];
	mb_area_rois =[laaroi[2] for laaroi in lengths_areas_and_arearois];
	# save membrane channel with original anchors, fixed anchors and membrane edge for assessment of performance
	mbfig.save_membrane_edge_image(membrane_channel_imp, calculated_objects.fixed_anchors_list, calculated_objects.membrane_edges, mb_area_rois, params);
	
	# perform analysis of background regions
	bg_rois = mb.generate_background_rois(None, params, calculated_objects.membrane_edges, threshold_method='MinError', membrane_imp=membrane_channel_imp); # resegment, allowing more conservative guess at bg region, or...
	#bg_rois = mb.generate_background_rois(segmentation_channel_imp, params, calculated_objects.membrane_edges); # use existing segmentation
	# do qc
	if params.qc_background_rois:
		bg_rois = mbui.qc_background_regions(actin_channel_imp, bg_rois)
	mbio.save_qcd_edges2(bg_rois, params.output_path, "background regions.zip")
	calculated_objects.background_sd_profile = mb.get_stddevs_by_frame_and_region(actin_channel_imp, bg_rois);

	t0_actin_mean = None;
	for fridx in range(membrane_channel_imp.getNFrames()):
		# generate curvature - this needs to be looped over slices
		membrane_edge = calculated_objects.membrane_edges[fridx];
		curv_profile = mb.calculate_curvature_profile(membrane_edge, params);
		curvature_profiles.append(curv_profile);
			
		# generate actin-channel line profile if actin channel present
		if actin_channel_imp is None:
			actin_channel_imp = Duplicator().run(membrane_channel_imp);
		actin_channel_imp.setPosition(fridx+1);
		actin_channel_imp, t0_actin_mean = mb.apply_photobleach_correction_framewise(params, 
																					actin_channel_imp, 
																					segmentation_channel_imp, 
																					t0_value=t0_actin_mean);
		actin_profile = mb.maximum_line_profile(actin_channel_imp, 
												membrane_edge, 
												int(round(params.intensity_profile_width_um / params.pixel_physical_size)));
		actin_profiles.append(actin_profile);

		if params.inner_outer_comparison:
			mean_intensity = float(sum([a for ((x, y), a) in actin_profile]))/len(actin_profile);
			calculated_objects.inner_outer_data.outer_means.append(mean_intensity) if repeat_fraction==1 else calculated_objects.inner_outer_data.inner_means.append(mean_intensity);
			sd = math.sqrt(sum((x - mean_intensity)**2 for ((x, y), a) in actin_profile) / len(actin_profile));
			calculated_objects.inner_outer_data.outer_sds.append(sd) if repeat_fraction==1 else calculated_objects.inner_outer_data.inner_sds.append(sd);

	calculated_objects.curvature_profiles = curvature_profiles;
	calculated_objects.actin_profiles = actin_profiles;
	calculated_objects.bleb_perimeter_lengths = mb_lengths;
	calculated_objects.bleb_areas = mb_areas;
	calculated_objects.full_membrane_lengths = full_membrane_lengths;
	calculated_objects.euclidean_membrane_lengths = euclidean_membrane_lengths;
	params.setPhysicalCurvatureUnit(params.pixel_unit + u'\u02C9' + u'\u00B9');
	
	return calculated_objects;

def generate_and_save_figures(imp, calculated_objects, params, membrane_channel_imp, segmentation_channel_imp):
	fig_imp_list = [];
	norm_curv_kym = mbfig.generate_kymograph(calculated_objects.curvature_profiles, params.curvature_kymograph_lut_string, "Curvature kymograph - distal point at middle");
	curv_kym = mbfig.generate_plain_kymograph(calculated_objects.curvature_profiles, params.curvature_kymograph_lut_string, "Curvature kymograph");
	fig_imp_list.append(norm_curv_kym);
	fig_imp_list.append(curv_kym);
	FileSaver(norm_curv_kym).saveAsTiff(os.path.join(params.output_path, "normalised position curvature kymograph.tif"));
	FileSaver(curv_kym).saveAsTiff(os.path.join(params.output_path, "raw curvature kymograph.tif"));
	#overlaid_curvature_imp, raw_curvature_imp = mbfig.overlay_curvatures(imp, curvature_stack, curvature_profiles, params.membrane_channel_number, params, annotate=True);
	overlaid_curvature_imp, raw_curvature_imp = mbfig.overlay_curvatures(imp, calculated_objects.curvature_profiles, params, annotate=True);
	fig_imp_list.append(overlaid_curvature_imp);
	fig_imp_list.append(raw_curvature_imp);
	if segmentation_channel_imp is not None:
		FileSaver(segmentation_channel_imp).saveAsTiff(os.path.join(params.output_path, "binary_membrane_stack.tif"));
	bleb_len_imp, bleb_ls = mbfig.plot_bleb_evolution(calculated_objects.timelist, 
											calculated_objects.bleb_perimeter_lengths, 
											"Edge length (" + params.pixel_unit + ")");
	bleb_a_imp, bleb_as = mbfig.plot_bleb_evolution(calculated_objects.timelist, 
											calculated_objects.bleb_areas, 
											"Bleb area (" + params.pixel_unit + "^2)");
	fig_imp_list.append(bleb_len_imp);
	fig_imp_list.append(bleb_a_imp);
	FileSaver(bleb_len_imp).saveAsTiff(os.path.join(params.output_path, "bleb perimeter length.tif"));
	FileSaver(bleb_a_imp).saveAsTiff(os.path.join(params.output_path, "bleb area.tif"));
	mbio.save_1d_profile_as_csv([(t, l) for t, l in zip(calculated_objects.timelist, calculated_objects.bleb_perimeter_lengths)], 
							 os.path.join(params.output_path, "bleb perimeter length.csv"), 
							 [("Time, " + params.interval_unit), "Length, " + params.pixel_unit]);
	mbio.save_1d_profile_as_csv([(t, l) for t, l in zip(calculated_objects.timelist, calculated_objects.full_membrane_lengths)], 
							 os.path.join(params.output_path, "full membrane length.csv"), 
							 [("Time, " + params.interval_unit), "Length, " + params.pixel_unit]);
	mbio.save_1d_profile_as_csv([(t, a) for t, a in zip(calculated_objects.timelist, calculated_objects.bleb_areas)], 
							 os.path.join(params.output_path, "bleb area.csv"), 
							 [("Time, " + params.interval_unit), "Area, " + params.pixel_unit + "^2"]);
	mbio.save_1d_profile_as_csv([(t, euc) for t, euc in zip(calculated_objects.timelist, calculated_objects.euclidean_membrane_lengths)], 
								 os.path.join(params.output_path, "full membrane euclidean length.csv"), 
								 ["Time, " + params.interval_unit, "Membrane euclidean length, " + params.pixel_unit])
		
	# actin channel
	if calculated_objects.actin_profiles is not None:
		if len(calculated_objects.actin_profiles) > 0:
			norm_actin_kym = mbfig.generate_kymograph(calculated_objects.actin_profiles, params.actin_kymograph_lut_string, (params.labeled_species + " intensity - distal point at middle"));
			actin_kym = mbfig.generate_plain_kymograph(calculated_objects.actin_profiles, params.actin_kymograph_lut_string, (params.labeled_species + " intensity"));
			fig_imp_list.append(norm_actin_kym);
			fig_imp_list.append(actin_kym);
			FileSaver(norm_actin_kym).saveAsTiff(os.path.join(params.output_path, "normalised position " + params.labeled_species + " kymograph.tif"));
			FileSaver(actin_kym).saveAsTiff(os.path.join(params.output_path, "raw " + params.labeled_species + " kymograph.tif"));
			
			mrg_imp = mbfig.merge_kymographs(norm_actin_kym, norm_curv_kym, params);
			fig_imp_list.append(mrg_imp);
			FileSaver(mrg_imp).saveAsTiff(os.path.join(params.output_path, "merged intensity and curvature kymograph.tif"));
			#mbfig.generate_intensity_weighted_curvature(raw_curvature_imp, curvature_profiles, actin_channel_imp, "physics");
	return fig_imp_list;

def save_csvs(calculated_objects, params):
	tname = "Time, " + params.interval_unit;
	mbio.save_profile_as_csv(calculated_objects.curvature_profiles, os.path.join(params.output_path, "curvatures.csv"), "curvature");
	mbio.save_profile_as_csv(calculated_objects.actin_profiles, 
									os.path.join(params.output_path, 
									(params.labeled_species + " intensities.csv")), 
									(params.labeled_species + " intensity"), 
									tname=tname, 
									time_list=calculated_objects.timelist);
	pixel_unit_for_encoding = "um" if u'\xb5m' in params.pixel_unit else params.pixel_unit;
	mbio.save_1d_profile_as_csv([(t, std) for t, std in zip(calculated_objects.timelist, calculated_objects.background_sd_profile)], 
								 os.path.join(params.output_path, (params.labeled_species + " channel background standard deviations.csv")), 
								 ["Time, " + params.interval_unit, params.labeled_species + " bg std"]);
	#mbio.save_1d_profile_as_csv([(t, flen) for t, flen in zip(calculated_objects.timelist, calculated_objects.full_membrane_lengths)], 
	#							 os.path.join(params.output_path, "full membrane perimeter lengths.csv"), 
	#							 ["Time, " + params.interval_unit, "Membrane perimeter length, {}".format(pixel_unit_for_encoding)])
	return;