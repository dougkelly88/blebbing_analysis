# functions for handling user interface for membrane blebbing analysis
#
# D. J. Kelly, 2018-10-15, douglas.kelly@riken.jp

# imports
import math
from ij import IJ
from ij.gui import WaitForUserDialog, GenericDialog, NonBlockingGenericDialog, Roi, PolygonRoi
from ij.plugin import Duplicator
from ij.process import AutoThresholder
from java.awt import GraphicsEnvironment, Panel, Dimension, Checkbox, CheckboxGroup
from javax.swing import Box
from loci.formats import ImageReader, MetadataTools
from loci.formats.gui import BufferedImageReader
from loci.plugins.in import ImporterOptions, ThumbLoader

from Parameters import Parameters
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
	IJ.setTool("rect");
	MyWaitForUser("Crop", "If desired, select an area ROI to crop to...");
	roi = imp.getRoi();
	crop_params = None;
	original_imp = imp.clone();
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
		autoset_zoom(imp);
		imp.killRoi();
		params.setSpatialCrop(crop_params)
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

def time_crop(imp):
	"""trim a time series based on interactively-defined start and end points"""
	MyWaitForUser("First T...", "Choose first time frame and click OK...");
	start_frame = imp.getT();
	MyWaitForUser("Last T...", "Now choose last time frame and click OK...");
	end_frame = imp.getT();
	slices = imp.getNSlices();
	channels = imp.getNChannels();
	dupimp = Duplicator().run(imp, 1, channels, 1, slices, start_frame, end_frame);
	imp.changes = False;
	imp.close();
	dupimp.show()
	autoset_zoom(dupimp);
	return dupimp, (start_frame, end_frame);

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

def perform_user_qc(imp, edges, alt_edges, fixed_anchors_list, params):
	"""allow the user to intervene to fix erroneously identified membrane edges"""
	output_folder = params.output_path;
	current_edges = edges;
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
	qcd_edges[imp.getT() - 1] = last_roi;
	imp.removeImageListener(listener);
	for fridx in range(0, imp.getNFrames()):
		if qcd_edges[fridx].getType() == Roi.FREELINE:
			if (fridx == 0):
				anchors = params.manual_anchor_positions;
			else:
				anchors = fixed_anchors_list[fridx - 1];
			qcd_edges[fridx] = mb.flip_edge(qcd_edges[fridx], anchors);
			fixed_anchors = mb.fix_anchors_to_membrane(anchors, qcd_edges[fridx], params);
			fixed_anchors_list[fridx] = fixed_anchors;
			poly =  qcd_edges[fridx].getInterpolatedPolygon(0.25, False);
			polypoints = [(x,y) for x,y in zip(poly.xpoints, poly.ypoints)];
			idx = [polypoints.index(fixed_anchors[0]), polypoints.index(fixed_anchors[1])];
			idx.sort();
			polypoints = polypoints[idx[0]:idx[1]];
			newedge = PolygonRoi([x for (x,y) in polypoints], 
									[y for (x,y) in polypoints], 
									Roi.POLYLINE);
			imp.setPosition(fridx + 1);
			imp.setRoi(newedge);
			IJ.run(imp, "Interpolate", "interval=1.0 smooth adjust");
			IJ.run(imp, "Fit Spline", "");
			qcd_edges[fridx] = imp.getRoi();
	mbio.save_qcd_edges(qcd_edges, output_folder);
	imp.close();
	return qcd_edges, fixed_anchors_list;

def analysis_parameters_gui():
	"""GUI for setting analysis parameters at the start of a run. TODO: more effectively separate model and view"""
	params = Parameters(load_last_params = True);
	dialog = GenericDialog("Analysis parameters");
	dialog.addNumericField("Curvature length parameter (um):", 
							round(params.curvature_length_um, 2), 
							2);
	dialog.addNumericField("Width of region for intensity analysis (um):", 
							round(params.intensity_profile_width_um, 2), 
							2);
	dialog.addChoice("Threshold method: ", 
						params.listThresholdMethods(),
						params.threshold_method);
	dialog.addChoice("Curvature overlay LUT: ", 
						IJ.getLuts(), 
						params.curvature_overlay_lut_string);
	dialog.addChoice("Curvature kymograph LUT: ", 
						IJ.getLuts(), 
						params.curvature_kymograph_lut_string);
	dialog.addChoice("Labelled species kymograph LUT: ", 
						IJ.getLuts(), 
						params.actin_kymograph_lut_string);
	dialog.addStringField("Labelled species for intensity analysis: ", 
							params.labeled_species);
	dialog.addCheckbox("Use intensity channel for segmentation too?", 
						params.use_single_channel);
	dialog.addRadioButtonGroup("Metadata source: ", 
								["Image metadata", "Acquisition metadata"], 
								1, 2, 
								params.metadata_source);
	dialog.addCheckbox("Filter out negative curvatures", 
						params.filter_negative_curvatures);
	dialog.addCheckbox("Account for photobleaching?", 
						params.photobleaching_correction);
	dialog.addCheckbox("Perform quality control of membrane edges?", 
						params.perform_user_qc);
	dialog.addCheckbox("Perform spatial cropping?", 
						params.perform_spatial_crop);
	dialog.addCheckbox("Perform time cropping?", 
						params.perform_time_crop);
	dialog.addCheckbox("Close images on completion?", 
						params.close_on_completion);
	dialog.addCheckbox("Compare inner and outer curvature regions?", 
						params.inner_outer_comparison);
	dialog.showDialog();
	if dialog.wasCanceled():
		raise KeyboardInterrupt("Run canceled");
	choices =  dialog.getChoices();
	params.setCurvatureLengthUm(dialog.getNextNumber()); # check whether label of numeric field is contained in getNextNumber?
	params.setIntensityProfileWidthUm(dialog.getNextNumber());
	params.setThresholdMethod(choices[0].getSelectedItem());
	params.setCurvatureOverlayLUT(choices[1].getSelectedItem());
	params.setCurvatureKymographLUT(choices[2].getSelectedItem());
	params.setActinKymographLUT(choices[3].getSelectedItem()); # similarly, whether getNextChoice has method to get label - this way, less dependent on order not changing...
	params.setLabeledSpecies(dialog.getNextString());
	params.setUseSingleChannel(dialog.getNextBoolean());
	params.setFilterNegativeCurvatures(dialog.getNextBoolean());
	params.togglePhotobleachingCorrection(dialog.getNextBoolean());
	params.togglePerformUserQC(dialog.getNextBoolean());
	params.toggleSpatialCrop(dialog.getNextBoolean());
	params.toggleTimeCrop(dialog.getNextBoolean());
	params.toggleCloseOnCompletion(dialog.getNextBoolean());
	params.setMetadataSource(dialog.getNextRadioButton());
	params.setDoInnerOuterComparison(dialog.getNextBoolean());
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