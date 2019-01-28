# functions for handling user interface for membrane blebbing analysis
#
# D. J. Kelly, 2018-10-15, douglas.kelly@riken.jp

# imports
import math
from ij import IJ, ImageStack, ImagePlus
from ij.gui import GenericDialog, NonBlockingGenericDialog, Roi, PolygonRoi, PointRoi
from ij.plugin import Duplicator
from java.awt import GraphicsEnvironment, Panel, Dimension, Checkbox, CheckboxGroup, Color
from javax.swing import Box
from loci.formats import ImageReader, MetadataTools
from loci.formats.gui import BufferedImageReader
from loci.plugins.in import ImporterOptions, ThumbLoader

from Parameters import Parameters
from MyControlDefinition import MyControlDefinition
from UpdateRoiImageListener import UpdateRoiImageListener
import membraneBlebbingEngine as mb
import membraneBlebbingFileio as mbio

def autoset_zoom(imp):
	"""set the zoom of the current imageplus to give a reasonable window size,  based on reasonable guess at screen resolution"""
	h = imp.getHeight();
	w = imp.getWidth();
	try:
		gd = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice().getDisplayMode();
		screen_h = gd.getHeight();
		screen_w = gd.getWidth();
	except:
		screen_h = 1080;
		screen_w = 1920;
	zy = 100* math.floor(screen_h / (2 * h));
	zx = 100* math.floor(screen_w / (2 * w));
	zoom_factor = min([zx, zy]);
	IJ.run("Set... ", "zoom=" + str(zoom_factor) + 
						" x=" + str(math.floor(w/2)) + 
						" y=" + str(math.floor(h/2)));
	IJ.run("Scale to Fit", "");

def crop_to_ROI(imp, params):
	"""prompt user to select ROI for subsequent analysis"""
	roi = None;
	if params.perform_spatial_crop:
		if params.spatial_crop is not None:
			roi = params.parse_roistr_to_roi();
			imp.setRoi(roi);
		else:
			IJ.setTool("rect");
			MyWaitForUser("Crop", "If desired, select an area ROI to crop to...");
			roi = imp.getRoi();
	crop_params = None;
	original_imp = Duplicator().run(imp);
	if roi is not None:
		if not roi.isArea():
			raise TypeError("selected ROI should be an area");
		if roi.getType():
			crop_params = str([(x, y) for x, y in 
							zip(roi.getPolygon().xpoints, roi.getPolygon().ypoints)]);
			IJ.run(imp, "Make Inverse", "");
			roi = imp.getRoi();
			fill_val = mb.calculate_percentile(imp, roi, 25);
			IJ.run(imp, "Set...", "value=" + str(round(fill_val)) + " stack");
			IJ.run(imp, "Make Inverse", "");
		else:
			crop_params = roi.getBounds().toString();
		IJ.run(imp, "Crop", "");
		imp.killRoi();
		params.setSpatialCrop(crop_params)
	else:
		params.setSpatialCrop(None);
	return original_imp, crop_params;

def prompt_for_points(imp, title, message, n_points):
	"""prompt the user to provide a number of points on the image"""
	imp.killRoi();
	if (n_points == 1): 
		IJ.setTool("point");
	else:
		IJ.setTool("multipoint");
	selected_points = 0;
	while (selected_points != n_points):
		MyWaitForUser(title, message);
		roi = imp.getRoi();
		if roi is not None:
			selected_points = len(roi.getContainedPoints());
		if ((roi is None) or (selected_points != n_points)):
			MyWaitForUser("Error!", "Wrong number of points selected! Please try again...");
			imp.killRoi();
	return [(p.x, p.y) for p in roi.getContainedPoints()];

def time_crop(imp, params):
	"""trim a time series based on interactively-defined start and end points"""
	start_frame = None;
	end_frame = None;
	if params.perform_time_crop:
		if params.time_crop_start_end is not None:
			start_frame = params.time_crop_start_end[0];
			end_frame = params.time_crop_start_end[1];
		else:
			MyWaitForUser("First T...", "Choose first time frame and click OK...");
			start_frame = imp.getT();
			MyWaitForUser("Last T...", "Now choose last time frame and click OK...");
			end_frame = imp.getT();
	slices = imp.getNSlices();
	channels = imp.getNChannels();
	if start_frame is not None:
		dupimp = Duplicator().run(imp, 1, channels, 1, slices, start_frame, end_frame);
		imp.changes = False;
		imp.close();
		dupimp.show()
		autoset_zoom(dupimp);
		return dupimp, (start_frame, end_frame);
	else:
		return imp, (start_frame, end_frame);

def warning_dialog(message):
	"""show a warning dialog with a user-defined message"""
	dialog = GenericDialog("Warning!");
	dialog.setCancelLabel("Cancel analysis");
	if type(message) is list:
		for line in message:
			dialog.addMessage(line);
	else:
		dialog.addMessage(message);
	dialog.showDialog();
	if dialog.wasCanceled():
		raise KeyboardInterrupt("Run canceled");
	return;

def MyWaitForUser(title, message):
	"""non-modal dialog with option to cancel the analysis"""
	dialog = NonBlockingGenericDialog(title);
	dialog.setCancelLabel("Cancel analysis");
	if type(message) is list:
		for line in message:
			dialog.addMessage(line);
	else:
		dialog.addMessage(message);
	dialog.showDialog();
	if dialog.wasCanceled():
		raise KeyboardInterrupt("Run canceled");
	return;

def perform_user_qc(in_imp, edges, alt_edges, fixed_anchors_list, params):
	"""allow the user to intervene to fix erroneously identified membrane edges"""
	output_folder = params.output_path;
	current_edges = edges;
	rgbstack = ImageStack(in_imp.getWidth(), in_imp.getHeight());
	for tidx in range(in_imp.getNFrames()): 
		in_imp.setT(tidx+1);
		ip = in_imp.getProcessor();
		rgbip = ip.convertToRGB();
		rgbstack.addSlice(rgbip);
	imp = ImagePlus(("RGB " + in_imp.getTitle()), rgbstack);
	IJ.run("Colors...", "foreground=red background=white selection=yellow");
	for tidx in range(imp.getNSlices()):
		imp.setSlice(tidx+1);
		for anchor in params.manual_anchor_positions:
			imp.setRoi(PointRoi(anchor[0], anchor[1]));
			IJ.run(imp, "Draw", "slice");
	imp.show();
	autoset_zoom(imp);
	imp.setPosition(1);
	imp.setRoi(current_edges[0]);
	listener = UpdateRoiImageListener(current_edges);
	imp.addImageListener(listener);
	IJ.setTool("freeline");
	do_flip = True;
	while do_flip:
		dialog = NonBlockingGenericDialog("User quality control");
		dialog.enableYesNoCancel("Continue", "Flip all edges");
		dialog.setCancelLabel("Cancel analysis");
		dialog.addMessage("Please redraw the membrane edges as necessary, \n" + 
						"making sure to draw beyond anchor points at either end...\n" + 
						"Click OK when done. ");
		dialog.showDialog();
		if dialog.wasCanceled():
			raise KeyboardInterrupt("Run canceled");
		elif dialog.wasOKed():
			do_flip = False;
		else:
			print("flip edges");
			do_flip = True;
			imp.removeImageListener(listener);
			current_edges = alt_edges if (current_edges == edges) else edges;
			imp.setPosition(1);
			imp.setRoi(current_edges[0]);
			listener = UpdateRoiImageListener(current_edges);
			imp.addImageListener(listener);

	last_roi = imp.getRoi();
	qcd_edges = listener.getRoiList();
	if imp.getNFrames() > imp.getNSlices():
		qcd_edges[imp.getT() - 1] = last_roi;
	else:
		qcd_edges[imp.getZ() - 1] = last_roi;
	imp.removeImageListener(listener);
	mbio.save_qcd_edges2(qcd_edges, output_folder);
	nframes = imp.getNFrames() if imp.getNFrames()>imp.getNSlices() else imp.getNSlices();
	for fridx in range(0, nframes):
		if (qcd_edges[fridx].getType()==Roi.FREELINE) or (qcd_edges[fridx].getType()==Roi.POLYLINE):
			if (fridx == 0) or params.constrain_anchors:
				anchors = params.manual_anchor_positions;
			else:
				anchors = fixed_anchors_list[fridx - 1];
			fixed_anchors = mb.fix_anchors_to_membrane(anchors, qcd_edges[fridx], params);
			fixed_anchors = mb.order_anchors(fixed_anchors, params.manual_anchor_midpoint);
			fixed_anchors_list[fridx] = fixed_anchors;
			poly =  qcd_edges[fridx].getInterpolatedPolygon(0.25, False);
			polypoints = [(x,y) for x,y in zip(poly.xpoints, poly.ypoints)];
			idx = [polypoints.index(fixed_anchors[0]), polypoints.index(fixed_anchors[1])];
			idx.sort();
			polypoints = polypoints[idx[0]:idx[1]];
			newedge = PolygonRoi([x for (x,y) in polypoints], 
									[y for (x,y) in polypoints], 
									Roi.POLYLINE);
			newedge = mb.check_edge_order(anchors, newedge);
			imp.setPosition(fridx + 1);
			imp.setRoi(newedge);
			IJ.run(imp, "Interpolate", "interval=1.0 smooth adjust");
			IJ.run(imp, "Fit Spline", "");
			qcd_edges[fridx] = imp.getRoi();
	mbio.save_qcd_edges2(qcd_edges, output_folder);
	imp.changes = False;
	imp.close();
	return qcd_edges, fixed_anchors_list;

def analysis_parameters_gui(rerun_analysis=False, params=None):
	"""GUI for setting analysis parameters at the start of a run. TODO: more effectively separate model and view"""
	if params is None:
		params = Parameters(load_last_params = True);
	dialog = GenericDialog("Analysis parameters");
	controls = [];
	if rerun_analysis:
		params.setUseSingleChannel(False);
		params.togglePerformUserQC(False);
		params.setDoInnerOuterComparison(False);

	controls.append(MyControlDefinition("Curvature length parameter (um):", 
									 MyControlDefinition.Numeric, 
									 round(params.curvature_length_um, 2), 
									 params.setCurvatureLengthUm));
	controls.append(MyControlDefinition("Width of region for intensity analysis (um):", 
									MyControlDefinition.Numeric, 
									round(params.intensity_profile_width_um, 2), 
									params.setIntensityProfileWidthUm))
	controls.append(MyControlDefinition("Threshold method: ",
									 MyControlDefinition.Choice, 
									 params.threshold_method, 
									 params.setThresholdMethod, 
									 choices=params.listThresholdMethods(), 
									 enabled=(not rerun_analysis)));
	controls.append(MyControlDefinition("Curvature overlay LUT: ", 
									 MyControlDefinition.Choice, 
									 params.curvature_overlay_lut_string, 
									 params.setCurvatureOverlayLUT, 
									 choices=IJ.getLuts()));
	controls.append(MyControlDefinition("Curvature kymograph LUT: ",
									 MyControlDefinition.Choice, 
									 params.curvature_kymograph_lut_string, 
									 params.setCurvatureKymographLUT, 
									 choices=IJ.getLuts()));
	controls.append(MyControlDefinition("Labelled species kymograph LUT: ", 
									 MyControlDefinition.Choice, 
									 params.actin_kymograph_lut_string, 
									 params.setActinKymographLUT, 
									 choices=IJ.getLuts()));
	controls.append(MyControlDefinition("Labelled species for intensity analysis: ", 
									 MyControlDefinition.String, 
									 params.labeled_species, 
									 params.setLabeledSpecies));
	controls.append(MyControlDefinition("Use intensity channel for segmentation too?", 
									 MyControlDefinition.Checkbox, 
									 params.use_single_channel, 
									 params.setUseSingleChannel, 
									 enabled=(not rerun_analysis)));
	controls.append(MyControlDefinition("Metadata source: ", 
									 MyControlDefinition.RadioButtonGroup, 
									 params.metadata_source, 
									 params.setMetadataSource, 
									 choices=["Image metadata", "Acquisition metadata"]));
	controls.append(MyControlDefinition("Constrain anchors close to manual selections?", 
									 MyControlDefinition.Checkbox, 
									 params.constrain_anchors, 
									 params.toggleConstrainAnchors, 
									 enabled=(not rerun_analysis)));
	controls.append(MyControlDefinition("Filter out negative curvatures", 
									 MyControlDefinition.Checkbox, 
									 params.filter_negative_curvatures, 
									 params.setFilterNegativeCurvatures));
	controls.append(MyControlDefinition("Account for photobleaching?", 
									 MyControlDefinition.Checkbox, 
									 params.photobleaching_correction, 
									 params.togglePhotobleachingCorrection));
	controls.append(MyControlDefinition("Perform quality control of membrane edges?", 
									 MyControlDefinition.Checkbox, 
									 params.perform_user_qc, 
									 params.togglePerformUserQC, 
									 enabled=(not rerun_analysis)));
	controls.append(MyControlDefinition("Perform spatial cropping?", 
									 MyControlDefinition.Checkbox, 
									 params.perform_spatial_crop, 
									 params.toggleSpatialCrop, 
									 enabled=(not rerun_analysis)));
	controls.append(MyControlDefinition("Perform time cropping?",
									 MyControlDefinition.Checkbox, 
									 params.perform_time_crop, 
									 params.toggleTimeCrop, 
									 enabled=(not rerun_analysis)));
	controls.append(MyControlDefinition("Close images on completion?", 
									 MyControlDefinition.Checkbox, 
									 params.close_on_completion, 
									 params.toggleCloseOnCompletion));
	controls.append(MyControlDefinition("Compare inner and outer curvature regions?", 
									 MyControlDefinition.Checkbox, 
									 params.inner_outer_comparison, 
									 params.setDoInnerOuterComparison, 
									 enabled=(not rerun_analysis)));
	for control in controls:
		control.addControl(dialog);
	dialog.showDialog();
	if dialog.wasCanceled():
		raise KeyboardInterrupt("Run canceled");

	numeric_controls = [c for c in controls if c.control_type==MyControlDefinition.Numeric];
	for nc, nf in zip(numeric_controls, dialog.getNumericFields()):
		nc.setter(float(nf.getText()));
	string_controls = [c for c in controls if c.control_type==MyControlDefinition.String];
	for sc, sf in zip(string_controls, dialog.getStringFields()):
		sc.setter(sf.getText());
	choice_controls = [c for c in controls if c.control_type==MyControlDefinition.Choice];
	for cc, cf in zip(choice_controls, dialog.getChoices()):
		cc.setter(cf.getSelectedItem());
	checkbox_controls = [c for c in controls if c.control_type==MyControlDefinition.Checkbox];
	for cbc, cbf in zip(checkbox_controls, dialog.getCheckboxes()):
		cbc.setter(cbf.getState());
	radiobuttongroup_controls = [c for c in controls if c.control_type==MyControlDefinition.RadioButtonGroup];
	for rbc in radiobuttongroup_controls:
		rbc.setter(rbc.checkboxes[[cb.getState() for cb in rbc.checkboxes].index(True)].getLabel());

	params.persistParameters();
	return params;

def choose_series(filepath, params):
	"""if input file contains more than one image series (xy position), prompt user to choose which one to use"""
	# todo: if necessary (e.g. if lots of series), can improve thumbnail visuals based loosely on https://github.com/ome/bio-formats-imagej/blob/master/src/main/java/loci/plugins/in/SeriesDialog.java
	import_opts = ImporterOptions();
	import_opts.setId(filepath);
	
	reader = ImageReader();
	ome_meta = MetadataTools.createOMEXMLMetadata();
	reader.setMetadataStore(ome_meta);
	reader.setId(filepath);
	no_series = reader.getSeriesCount();
	if no_series == 1:
		return import_opts, params;
	else:
		series_names = [ome_meta.getImageName(idx) for idx in range(no_series)];
		dialog = GenericDialog("Select series to load...");
		dialog.addMessage("There are multiple series in this file! \n" + 
						"This is probably because there are multiple XY stage positions. \n " + 
						"Please choose which series to load: ");
		thumbreader = BufferedImageReader(reader);
		cbg = CheckboxGroup();
		for idx in range(no_series):
			p = Panel();
			p.add(Box.createRigidArea(Dimension(thumbreader.getThumbSizeX(), thumbreader.getThumbSizeY())));
			ThumbLoader.loadThumb(thumbreader, idx, p, True);
			dialog.addPanel(p);
			cb = Checkbox(series_names[idx], cbg, idx==0);
			p.add(cb);

		dialog.showDialog();
		if dialog.wasCanceled():
			raise KeyboardInterrupt("Run canceled");
		if dialog.wasOKed():
			selected_item = cbg.getSelectedCheckbox().getLabel();
			selected_index = series_names.index(selected_item);
			params.setSelectedSeriesIndex(selected_index);
			for idx in range(0, no_series):
				import_opts.setSeriesOn(idx, True) if (idx==selected_index) else import_opts.setSeriesOn(idx, False);
	reader.close();
	return import_opts, params

def crop_review():
	"""handle UI for reviewing cropping"""
	print("doing crop review...");
	dialog = NonBlockingGenericDialog("Review cropping")
	dialog.enableYesNoCancel("Keep this crop", "Revert to uncropped image");
	dialog.setCancelLabel("Cancel analysis");
	dialog.addMessage("Please check whether this cropping is as expected, \n" + 
					"and choose whether to press on or revert to using the \n" + 
					"full, uncropped image. ");
	dialog.showDialog();
	if dialog.wasCanceled():
		raise KeyboardInterrupt("Run canceled");
	elif dialog.wasOKed():
		keep_cropping = True;
	else:
		keep_cropping = False;
	return keep_cropping;
