# functions for handling user interface for membrane blebbing analysis
#
# D. J. Kelly, 2018-10-15, douglas.kelly@riken.jp

# imports
import math
from ij import IJ
from ij.gui import WaitForUserDialog, GenericDialog
from ij.plugin import Duplicator
from ij.process import AutoThresholder

from Parameters import Parameters

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
def crop_to_ROI(imp):
	IJ.setTool("rect");
	WaitForUserDialog("Crop", "If desired, select a rectangular ROI to crop to...").show();
	roi = imp.getRoi();
	crop_params = None;
	original_imp = imp.clone();
	if roi is not None:
		IJ.run(imp, "Crop", "");
		autoset_zoom(imp);
		crop_params = roi.getBounds();
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
		WaitForUserDialog(title, message).show();
		roi = imp.getRoi();
		if roi is not None:
			selected_points = len(roi.getContainedPoints());
		if ((roi is None) or (selected_points != n_points)):
			WaitForUserDialog("Error!", "Wrong number of points selected! Please try again...").show();
			imp.killRoi();
	return [(p.x, p.y) for p in roi.getContainedPoints()];

def time_crop(imp):
	WaitForUserDialog("Choose first time frame and click OK...").show();
	start_frame = imp.getT();
	WaitForUserDialog("Now choose last time frame and click OK...").show();
	end_frame = imp.getT();
	slices = imp.getNSlices();
	channels = imp.getNChannels();
	cropimp = Duplicator().run(imp, 1, channels, 1, slices, start_frame, end_frame)
	return cropimp, (start_frame, end_frame);

def analysis_parameters_gui():
	"""GUI for setting analysis parameters at the start of a run. TODO: more effectively separate model and view"""
	params = Parameters(load_last_params = True);
	dialog = GenericDialog("Analysis parameters");
	dialog.addNumericField("Curvature length parameter (pix):", 
							params.curvature_length_pix, 
							1);
	dialog.addChoice("Threshold method: ", 
						AutoThresholder.getMethods(),
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
	dialog.addCheckbox("Filter out negative curvatures", 
						params.filter_negative_curvatures);
	dialog.addCheckbox("Perform spatial cropping?", 
						params.perform_spatial_crop)
	#dialog.addToSameRow();
	#dialog.addCheckbox("Perform time cropping?", 
	#					params.perform_time_crop)
	dialog.showDialog();
	if dialog.wasCanceled():
		raise KeyboardInterrupt("Run canceled");
	chc =  dialog.getChoices();
	params.setCurvatureLengthPix(dialog.getNextNumber()); # check whether label of numeric field is contained in getNextNumber?
	params.setThresholdMethod(chc[0].getSelectedItem());
	params.setCurvatureOverlayLUT(chc[1].getSelectedItem());
	params.setCurvatureKymographLUT(chc[2].getSelectedItem());
	params.setActinKymographLUT(chc[3].getSelectedItem()); # similarly, whether getNextChoice has method to get label - this way, less dependent on order not changing...
	params.setLabeledSpecies(dialog.getNextString());
	params.setFilterNegativeCurvatures(dialog.getNextBoolean());
	params.toggleSpatialCrop(dialog.getNextBoolean());
	#params.toggleTimeCrop(dialog.getNextBoolean());
	params.persistParameters();
	return params;

#out = analysis_parameters_gui();
#print(out.__dict__);