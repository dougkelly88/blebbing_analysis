# Code to load, view and edit edges
# Motivation:	strange results are seen with shortened edges. This should help figure out why. 
# D.J. Kelly, 2019-01-24, douglas.kelly@riken.jp

# imports
import os, sys, math
from ij import IJ, ImageStack, ImagePlus
from ij.io import FileSaver
from ij.gui import WaitForUserDialog, PointRoi, NonBlockingGenericDialog
from ij.plugin import ChannelSplitter, Duplicator, ZProjector
from ij import Prefs
from loci.plugins import BF as bf
from java.awt import Color

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
from CalculatedObjects import CalculatedObjects

def main():
	# ensure consistent preference settings
	Prefs.blackBackground = False;
	
	# get locations for previous and new outputs
	params = Parameters();
	params.loadLastParams();
	output_folder_old, output_folder = mbio.rerun_location_chooser(params.input_image_path);
	params.loadParametersFromJson(os.path.join(output_folder_old, 'parameters used.json'));
	params.setOutputPath(output_folder);

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

	# prompt user to select ROI
	original_imp = Duplicator().run(imp);
	_, crop_params = mbui.crop_to_ROI(imp, params);
	imp.show();

	if crop_params is not None:
		params.perform_spatial_crop = True;
		mbui.autoset_zoom(imp);
		imp.updateAndDraw();
		review_crop = mb.check_cropping(output_folder_old, params);
		keep_crop = not review_crop;
		if review_crop:
			keep_crop = mbui.crop_review();
		if not keep_crop:
			imp.changes = False;
			imp.close();
			imp = original_imp;
		else:
			original_imp.close();
	else:
		original_imp.close();

	# prompt user to do time cropping
	imp, start_end_tuple = mbui.time_crop(imp, params);
	params.setTimeCropStartEnd(start_end_tuple);

	params = mbio.get_metadata(params);
	params.setCurvatureLengthUm(round(params.curvature_length_um / params.pixel_physical_size) * params.pixel_physical_size);
	params.persistParameters();
	IJ.run(imp, "Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
	imp.show();
	if imp.getNChannels() > 1:
		imp.setPosition(params.membrane_channel_number, 1, 1);
	mbui.autoset_zoom(imp);
	IJ.run("Enhance Contrast", "saturated=0.35");

	split_channels = mbfs.split_image_plus(imp, params);
	membrane_test_channel_imp = Duplicator().run(split_channels[0]);
	segmentation_channel_imp = split_channels[-1];

	# import edges
	imp.hide();
	membrane_edges = mbio.load_qcd_edges(os.path.join(output_folder_old, "user_defined_edges.json"));
	dummy_anchors = [params.manual_anchor_positions for _ in membrane_edges]
	membrane_edges, fixed_anchors = mbui.perform_user_qc(membrane_test_channel_imp, membrane_edges, membrane_edges, dummy_anchors, params);
	imp.show();
	rgbstack = ImageStack(imp.getWidth(), imp.getHeight());
	for tidx in range(imp.getNFrames()): 
		imp.setT(tidx+1);
		ip = imp.getProcessor();
		rgbip = ip.convertToRGB();
		rgbstack.addSlice(rgbip);
	rgbimp = ImagePlus(("RGB " + imp.getTitle()), rgbstack);
	imp.close();
	rgbimp.show();
	IJ.run("Colors...", "foreground=yellow background=white selection=yellow");
	for tidx in range(rgbimp.getNSlices()):
		rgbimp.setSlice(tidx+1);
		rgbimp.setRoi(membrane_edges[tidx]);
		IJ.run(rgbimp, "Draw", "slice");
	IJ.run("Colors...", "foreground=red background=white selection=yellow");
	for tidx in range(rgbimp.getNSlices()):
		rgbimp.setSlice(tidx+1);
		for anchor in params.manual_anchor_positions:
			rgbimp.setRoi(PointRoi(anchor[0], anchor[1]));
			IJ.run(rgbimp, "Draw", "slice");
	params.saveParametersToJson(os.path.join(params.output_path, "parameters used.json"));
	params.persistParameters();	
	rgbimp.changes = False;

# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main()

