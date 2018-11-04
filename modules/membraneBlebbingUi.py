# functions for handling user interface for membrane blebbing analysis
#
# D. J. Kelly, 2018-10-15, douglas.kelly@riken.jp

# imports
import math
from ij import IJ
from ij.gui import WaitForUserDialog, GenericDialog, NonBlockingGenericDialog
from ij.plugin import Duplicator
from ij.process import AutoThresholder

from Parameters import Parameters
import membraneBlebbingEngine as mb

# set the zoom of the current imageplus to give a reasonable window size, 
# based on reasonable guess at screen resolution
def autoset_zoom(imp):
	h = imp.getHeight();
	w = imp.getWidth();
	zy = 100* math.floor(1080 / (2 * h));
	zx = 100* math.floor(1920 / (2 * w));
	zoom_factor = min([zx, zy]);
	IJ.run("Set... ", "zoom=" + str(zoom_factor) + 
						" x=" + str(math.floor(w/2)) + 
						" y=" + str(math.floor(h/2)));
	IJ.run("Scale to Fit", "");

# prompt user to select ROI for subsequent analysis
def crop_to_ROI(imp, params):
	IJ.setTool("rect");
	MyWaitForUser("Crop", "If desired, select an area ROI to crop to...");
	roi = imp.getRoi();
	if not roi.isArea():
		raise TypeError("selected ROI should be an area");
	crop_params = None;
	original_imp = imp.clone();
	if roi is not None:
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

# prompt the user to provide a number of points on the image
def prompt_for_points(imp, title, message, n_points):
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
	MyWaitForUser("Choose first time frame and click OK...");
	start_frame = imp.getT();
	MyWaitForUser("Now choose last time frame and click OK...");
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

def MyWaitForUser(title, message):
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

def analysis_parameters_gui():
	"""GUI for setting analysis parameters at the start of a run. TODO: more effectively separate model and view"""
	params = Parameters(load_last_params = True);
	dialog = GenericDialog("Analysis parameters");
	dialog.addNumericField("Curvature length parameter (pix):", 
							params.curvature_length_pix, 
							1);
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
	dialog.addRadioButtonGroup("Metadata source: ", 
								["Image metadata", "Acquisition metadata"], 
								1, 2, 
								params.metadata_source);
	dialog.addCheckbox("Filter out negative curvatures", 
						params.filter_negative_curvatures);
	dialog.addCheckbox("Account for photobleaching?", 
						params.photobleaching_correction)
	dialog.addCheckbox("Perform spatial cropping?", 
						params.perform_spatial_crop)
	dialog.addCheckbox("Perform time cropping?", 
						params.perform_time_crop);
	dialog.addCheckbox("Close images on completion?", 
						params.close_on_completion)
	dialog.showDialog();
	if dialog.wasCanceled():
		raise KeyboardInterrupt("Run canceled");
	choices =  dialog.getChoices();
	params.setCurvatureLengthPix(dialog.getNextNumber()); # check whether label of numeric field is contained in getNextNumber?
	params.setThresholdMethod(choices[0].getSelectedItem());
	params.setCurvatureOverlayLUT(choices[1].getSelectedItem());
	params.setCurvatureKymographLUT(choices[2].getSelectedItem());
	params.setActinKymographLUT(choices[3].getSelectedItem()); # similarly, whether getNextChoice has method to get label - this way, less dependent on order not changing...
	params.setLabeledSpecies(dialog.getNextString());
	params.setFilterNegativeCurvatures(dialog.getNextBoolean());
	params.togglePhotobleachingCorrection(dialog.getNextBoolean());
	params.toggleSpatialCrop(dialog.getNextBoolean());
	params.toggleTimeCrop(dialog.getNextBoolean());
	params.toggleCloseOnCompletion(dialog.getNextBoolean());
	params.setMetadataSource(dialog.getNextRadioButton());
	params.persistParameters();
	return params;
