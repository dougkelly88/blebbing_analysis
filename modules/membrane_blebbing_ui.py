# functions for handling user interface for membrane blebbing analysis
#
# D. J. Kelly, 2018-10-15, douglas.kelly@riken.jp

# imports
import math;
from ij import IJ;
from ij.gui import WaitForUserDialog

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
