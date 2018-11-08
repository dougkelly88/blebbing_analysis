# functions for handling output image generation for membrane blebbing analysis
#
# D. J. Kelly, 2018-10-15, douglas.kelly@riken.jp

# TODO?: refactor as a class so that it is trivial to configure whether images should be
# shown globally, e.g. fig = FigClass(show_figs = False); fig.overlay_curvatures(...)
# also possibly easier in this case to avoid boilerplate in image saving...

# imports
import os, sys
script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_path);

from ij import IJ, ImagePlus, ImageStack
from ij.gui import Plot, GenericDialog, TextRoi
from ij.io import FileSaver
from ij.measure import Measurements
from ij.plugin import RGBStackMerge, Duplicator
from ij.plugin.frame import RoiManager
from ij.process import FloatProcessor, ByteProcessor, ColorProcessor
from java.awt import Font, Color
import membraneBlebbingEngine as mb
import membraneBlebbingUi as mbui

def pause_for_debug():
	gd = GenericDialog("Continue?");
	gd.showDialog();
	if gd.wasCanceled():
		raise Exception("Run interupted");

def generate_curvature_overlays(curvature_profile, curvature_stack): 
	"""Generate overlays to display curvature"""
	w = curvature_stack.getWidth();
	ip = FloatProcessor(w, curvature_stack.getHeight());
	pix = ip.getPixels();
	mx = max([c for ((x, y), c) in curvature_profile]);
	for ((x, y), c) in curvature_profile:
		pix[int(round(y)) * w + int(round(x))] = c;
	curvature_stack.addSlice(ip);
	return curvature_stack;

def add_colorbar(imp, limits, fraction=0.05):
	"""destructively overwrite <fraction> at right of image with colorbar"""
	for fridx in range(0, imp.getNSlices()):
		imp.setPosition(fridx + 1);
		pix = imp.getProcessor().getPixels();
		w = imp.getWidth();
		h = imp.getHeight();
		for xidx in range(w - int(w*fraction), w):
			for yidx in range(0, h):
				pix[yidx * w + xidx] = float(limits[1] - limits[0]) * (float((h-1) - yidx)/(h-1)) + limits[0];
	return imp

def generate_limit_labels(imp, limits, cb_fraction, params):
	"""generate text ROIs in correct positions to label colorbar"""
	rois = [];
	w = imp.getWidth();
	h = imp.getHeight();
	txt_h_px = float(h)/15;
	txt_sz_pt = int(round((18.0/24.0) * txt_h_px));
	for limit in limits:
		roi = TextRoi(1, 1, str(round(limit, 3)), Font("SansSerif", Font.ITALIC, txt_sz_pt));
		roi.setJustification(TextRoi.RIGHT);
		roi.setFillColor(Color.BLACK);
		roi.setStrokeColor(Color.WHITE);
		rois.append(roi);
	roim = RoiManager(False);
	roim.reset();
	imp.show();
	mbui.autoset_zoom(imp);
	for fridx in range(1, imp.getNSlices()+1):
		imp.setPosition(fridx);
		roi_uu = rois[1].clone();
		xpos = w - cb_fraction * w - float(w)/100 - rois[1].getFloatWidth();
		roi_uu.setLocation(xpos, 1);
		imp.setRoi(roi_uu);
		roim.addRoi(roi_uu);
		roi_ll = rois[0].clone();
		roi_ll.setLocation(xpos, h - rois[0].getFloatHeight());
		imp.setRoi(roi_ll);
		roim.addRoi(roi_ll);
	roim.runCommand("Show All");
	FileSaver(imp).saveAsTiffStack(os.path.join(params.output_path, "overlaid curvature nudged labels.tif"));
	# nudge positions
	roim.reset();
	imp.killRoi();
	for fridx in range(1, imp.getNSlices()+1):
		imp.setPosition(fridx);
		roi_uu = rois[1].clone();
		xpos = w - cb_fraction * w - float(w)/100;
		roi_uu.setLocation(xpos, 1);
		imp.setRoi(roi_uu);
		roim.addRoi(roi_uu);
		roi_ll = rois[0].clone();
		roi_ll.setLocation(xpos, h - rois[0].getFloatHeight());
		imp.setRoi(roi_ll);
		roim.addRoi(roi_ll);
	roim.runCommand("Show All");
	FileSaver(imp).saveAsTiffStack(os.path.join(params.output_path, "overlaid curvature.tif"));
	return imp;
	
def overlay_curvatures(imp, curvature_stack, curvature_profiles, membrane_channel, params, limits = None, annotate=True):
	"""Overlay curvature pixels on membrane image"""
	overlay_base_imp = imp.clone();
	overlay_imp = ImagePlus("Curvature stack", curvature_stack);
	IJ.run(overlay_imp, params.curvature_overlay_lut_string, "");
	if limits is None:
		flat_list_curv = [c[1] for cs in curvature_profiles for c in cs];
		limits = [min(flat_list_curv), max(flat_list_curv)];
	IJ.setMinAndMax(overlay_imp, limits[0], limits[1]);
	raw_overlay = Duplicator().run(overlay_imp);
	if annotate:
		cb_fraction = 0.05;
		overlay_imp = add_colorbar(overlay_imp, limits, cb_fraction);
	IJ.run(overlay_imp, "RGB Color", "");
	overlaid_stack = ImageStack(overlay_imp.width, overlay_imp.height);
	for fridx in range(1, curvature_stack.getSize()+1):
		raw_idx = overlay_base_imp.getStackIndex(membrane_channel, 1, fridx);
		ip = overlay_base_imp.getStack().getProcessor(raw_idx).convertToRGB();
		pix = overlay_imp.getStack().getProcessor(fridx).getPixels();
		base_pix = ip.getPixels();
		for ((x, y), c) in curvature_profiles[fridx-1]:
			if params.filter_negative_curvatures:
				if (c > 0):
					base_pix[int(round(y)) * imp.width + int(round(x))] = pix[int(round(y)) * imp.width + int(round(x))];
			else: 
				base_pix[int(round(y)) * imp.width + int(round(x))] = pix[int(round(y)) * imp.width + int(round(x))];
		if annotate:
			w = overlay_base_imp.getWidth();
			h = overlay_base_imp.getHeight();
			for x in range(w - int(w * cb_fraction), w):
				for y in range(0, h):
					base_pix[int(round(y)) * imp.width + int(round(x))] = pix[int(round(y)) * imp.width + int(round(x))];
		overlaid_stack.addSlice(ip);
	out_imp = ImagePlus("Overlaid curvatures", overlaid_stack);
	if annotate:
		out_imp = generate_limit_labels(out_imp, limits, cb_fraction, params)
	FileSaver(raw_overlay).saveAsTiffStack(os.path.join(params.output_path, "raw curvature.tif"));
	return out_imp, raw_overlay;

def generate_plain_kymograph(data_to_plot, colormap_string, title_string):
	"""Display unnormalised kymograph"""
	ip = FloatProcessor(len(data_to_plot), max([len(data) for data in data_to_plot]));
	pix = ip.getPixels();
	for idx, data in enumerate(data_to_plot):
		for yidx in range(0,len(data)):
			pix[yidx * len(data_to_plot) + idx] = data[yidx][1];
	imp = ImagePlus(title_string, ip);
	IJ.run(imp, colormap_string, "");
	imp.show();
	return imp;

def generate_kymograph(data_to_plot, colormap_string, title_string):
	"""Display one-channel kymograph with point furthest from the edges along the middle of the kymograph """
	kym_height = 2 * max([len(data) for data in data_to_plot]) + 1;
	ip = FloatProcessor(len(data_to_plot), kym_height);
	# normalise such that point furthest from the anchors is in the middle of the kymograph
	for idx, data in enumerate(data_to_plot):
		dist = [mb.vector_length(data[0][0], p) *
				mb.vector_length(data[-1][0], p)
				  for p in [d[0] for d in data]];
		distal_idx = dist.index(max(dist));
		pix = ip.getPixels();
		for kidx, didx in zip(range(((kym_height - 1)/2 + 1), ((kym_height - 1)/2 + 1) + len(data) - distal_idx), 
						range(distal_idx, len(data))):
			pix[kidx * len(data_to_plot) + idx] = data[didx][1];
		for kidx, didx in zip(range(((kym_height - 1)/2 + 1) - distal_idx, ((kym_height - 1)/2 + 1)), 
						range(0, distal_idx)):
			pix[kidx * len(data_to_plot) + idx] = data[didx][1];
	imp = ImagePlus(title_string, ip);
	IJ.run(imp, colormap_string, "")
	imp.show();
	return imp;

def generate_intensity_weighted_curvature(curvature_overlay, curvature_profiles, intensity_channel_imp, colormap_string):
	"""Generate intensity-weighted curvature image"""
	curv_impRGB = curvature_overlay.clone();
	IJ.run(curv_impRGB, colormap_string, "");
	IJ.run(curv_impRGB, "RGB Color", "");

	int_imp16 = intensity_channel_imp.clone();

	base_impRGB = intensity_channel_imp.clone();
	IJ.run(base_impRGB, "Grays", "");
	IJ.run(base_impRGB, "RGB Color", "");
	
	maxes = [];
	for fridx in range(0, int_imp16.getNFrames()):
		int_imp16.setPositionWithoutUpdate(1, 1, fridx +1);
		maxes.append(int_imp16.getStatistics(Measurements.MIN_MAX).max)
	mx = float(max(maxes));

	for idx, profile in enumerate(curvature_profiles):
		print("Frame = " + str(idx));
		curv_impRGB.setPosition(idx+1);
		curvCP = curv_impRGB.getProcessor();
		base_impRGB.setPosition(idx+1);
		baseCP = base_impRGB.getProcessor();
		int_imp16.setPosition(idx+1);

		for chidx in range(0,3):
			c = ['r', 'g', 'b'];
			print("Image channel = " + c[chidx]);
			baseBP = ByteProcessor(base_impRGB.getWidth(), base_impRGB.getHeight());
			curvBP = ByteProcessor(base_impRGB.getWidth(), base_impRGB.getHeight());
			baseBP = baseCP.getChannel(chidx, baseBP);
			curvBP = curvCP.getChannel(chidx, curvBP);
			
			for ((x,y), c) in profile:
				x = int(round(x));
				y = int(round(y));
				baseBP.putPixelValue(x, y, int(curvBP.getPixel(x,y) * float(int_imp16.getPixel(x,y)[0])/mx));
				#baseBP.putPixelValue(x, y, int(curvBP.getPixel(x,y)));
			baseCP.setChannel(chidx, baseBP);
		base_impRGB.setProcessor(baseCP);
	base_impRGB.setTitle("Merged");
	base_impRGB.show();
	curv_impRGB.show();
	pause_for_debug();

def plot_bleb_evolution(ts, prop_to_plots, title):
	"""Plot evolution of a given bleb property over time"""
	title = title.replace("\u00B5", "u");
	plt = Plot((title + " against time"), "Time", title);
	plt.add("line", ts, prop_to_plots);
	plt.show();
	plot_data = [(t, d) for t, d in zip(ts, prop_to_plots)];
	return plt.getImagePlus(), plot_data;

def merge_kymographs(kym1_imp, kym2_imp, params):
	"""Merge two kymographs"""
	mrg_imp = RGBStackMerge().mergeChannels([kym1_imp, kym2_imp], True);
	mrg_imp.setTitle("Merged " + params.labeled_species + " intensity and curvature kymograph");
	mrg_imp.show();
	return mrg_imp;

#from ij import WindowManager as WM
#import membrane_blebbing_fileio as mbio;
#from ij.plugin import ChannelSplitter

#folder = "D:\\data\\Inverse blebbing\\output\\2018-10-29 17-06-15 output";
#curvature_overlay = IJ.openImage(os.path.join(folder, "raw curvature.tif"));
#curvature_profiles = mbio.load_csv_as_profile(os.path.join(folder, "curvatures.csv"));
#imp = IJ.openImage("D:\\data\\Inverse blebbing\\MAX_2dpf marcksl1b-EGFP inj_TgLifeact-mCh_movie e4_split-bleb1.tif");
#intensity_channel_imp = ChannelSplitter.split(imp)[1];
#intensity_channel_imp.show();
#colormap_string = "physics"

#generate_intensity_weighted_curvature(curvature_overlay, curvature_profiles, intensity_channel_imp, colormap_string);