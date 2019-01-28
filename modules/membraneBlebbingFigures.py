# functions for handling output image generation for membrane blebbing analysis
#
# D. J. Kelly, 2018-10-15, douglas.kelly@riken.jp

# imports
import os
#import os, sys
#script_path = os.path.dirname(os.path.realpath(__file__))
#sys.path.insert(0, script_path);

from ij import IJ, ImagePlus, ImageStack
from ij.gui import Plot, GenericDialog, TextRoi, PointRoi, Roi
from ij.io import FileSaver
from ij.measure import Measurements
from ij.plugin import RGBStackMerge, Duplicator
from ij.plugin.frame import RoiManager
from ij.process import FloatProcessor, ByteProcessor
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
	h = curvature_stack.getHeight();
	ip = FloatProcessor(w, h);
	pix = ip.getPixels();
	curv = [c for ((x, y), c) in curvature_profile];
	start_idx = next((i for i, x in enumerate(curv) if x), None);
	end_idx = next((len(curv) - i for i, x in enumerate(reversed(curv)) if x), None);
	for ((x, y), c) in curvature_profile[start_idx:end_idx]:
		# ensure that no rounding issues cause pixels to fall outside image...
		if x > (w - 1):
			x = w - 1;
		if y > (h - 1):
			y = h - 1;
		if x < 0:
			x = 0;
		if y < 0:
			y = 0;
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
	
def overlay_curvatures(imp, curvature_profiles, params, limits = None, annotate=True):
	"""Overlay curvature pixels on membrane image"""
	membrane_channel = params.membrane_channel_number;
	overlay_base_imp = imp.clone();
	curvature_stack = ImageStack(imp.getWidth(), imp.getHeight());
	for profile in curvature_profiles:
		curvature_stack = generate_curvature_overlays(profile, curvature_stack);
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
	w = overlay_imp.getWidth();
	h = overlay_imp.getHeight();
	overlaid_stack = ImageStack(w, h);
	for fridx in range(1, curvature_stack.getSize()+1):
		raw_idx = overlay_base_imp.getStackIndex(membrane_channel, 1, fridx);
		ip = overlay_base_imp.getStack().getProcessor(raw_idx).convertToRGB();
		pix = overlay_imp.getStack().getProcessor(fridx).getPixels();
		base_pix = ip.getPixels();
		curv = [c for ((x, y), c) in curvature_profiles[fridx-1]];
		start_idx = next((i for i, x in enumerate(curv) if x), None);
		end_idx = next((len(curv) - i for i, x in enumerate(reversed(curv)) if x), None);
		for ((x, y), c) in curvature_profiles[fridx-1][start_idx:end_idx]:
			# ensure that no rounding issues cause pixels to fall outside image...
			if x > (w - 1):
				x = w - 1;
			if y > (h - 1):
				y = h - 1;
			if x < 0:
				x = 0;
			if y < 0:
				y = 0;
			if params.filter_negative_curvatures:
				if (c > 0):
					base_pix[int(round(y)) * w + int(round(x))] = pix[int(round(y)) * w + int(round(x))];
			else: 
				base_pix[int(round(y)) * w + int(round(x))] = pix[int(round(y)) * w + int(round(x))];
		if annotate:
			for x in range(w - int(w * cb_fraction), w):
				for y in range(0, h):
					base_pix[int(round(y)) * w + int(round(x))] = pix[int(round(y)) * w + int(round(x))];
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

def generate_kymograph(data_to_plot, colormap_string, title_string, trim=True):
	"""Display one-channel kymograph with point furthest from the edges along the middle of the kymograph """
	kym_height = 2 * max([len(data) for data in data_to_plot]) + 1;
	ip = FloatProcessor(len(data_to_plot), kym_height);
	# normalise such that point furthest from the anchors is in the middle of the kymograph
	maxy = 0; 
	miny = kym_height;
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
		maxy = ((kym_height - 1)/2 + 1) + len(data) - distal_idx if (((kym_height - 1)/2 + 1) + len(data) - distal_idx) > maxy else maxy;
		miny = ((kym_height - 1)/2 + 1) - distal_idx if (((kym_height - 1)/2 + 1) - distal_idx) < miny else miny;
	imp = ImagePlus(title_string, ip);
	IJ.run(imp, colormap_string, "")
	if trim:
		maxtomiddle = maxy - (kym_height - 1)/2 + 1;
		mintomiddle = (kym_height - 1)/2 + 1 - miny;
		if maxtomiddle > mintomiddle:
			miny = (kym_height - 1)/2 + 1 - maxtomiddle;
		else:
			maxy = (kym_height - 1)/2 + 1 + mintomiddle;
		if (maxy - miny)/2 == round((maxy - miny)/2):
			maxy = maxy - 1;
		imp.setRoi(Roi(0,miny,imp.getWidth(),(maxy-miny)));
		imp = imp.crop();
	imp.show();
	return imp;

def generate_intensity_weighted_curvature(curvature_overlay, curvature_profiles, intensity_channel_imp, colormap_string):
	"""Generate intensity-weighted curvature image"""
	w = intensity_channel_imp.getWidth();
	h = intensity_channel_imp.getHeight();
	
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
				# ensure that no rounding issues cause pixels to fall outside image...
				if x > (w - 1):
					x = w - 1;
				if y > (h - 1):
					y = h - 1;
				if x < 0:
					x = 0;
				if y < 0:
					y = 0;
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

def save_membrane_edge_image(membrane_channel_imp, fixed_anchors_list, membrane_edges, area_rois, params):
	"""Save an image with the membrane channel overlaid with original anchor positions, fixed anchor positions, and membrane edge"""
	IJ.run(membrane_channel_imp, "RGB Color", "");
	anchors = params.manual_anchor_positions;
	midpoint = params.manual_anchor_midpoint[0];
	#w = membrane_channel_imp.getWidth();
	#h = membrane_channel_imp.getHeight();
	for fridx in range(0, membrane_channel_imp.getNFrames()):
		membrane_channel_imp.setPosition(fridx + 1);
		#framep = membrane_channel_imp.getProcessor();
		#rpix = framep.toFloat(0,None).getPixels();
		#bpix = framep.toFloat(2,None).getPixels();
		#float_poly = area_rois[fridx].getContainedFloatPoints();
		#for x, y in zip(float_poly.xpoints, float_poly.ypoints):
		#	xidx = int(round(x));
		#	yidx = int(round(y));
		#	rpix[(yidx + 1) * w + (xidx + 1)] = 0;
		#	bpix[(yidx + 1) * w + (xidx + 1)] = 0;
		#framep.setPixels(0, FloatProcessor(w, h, rpix));
		#framep.setPixels(2, FloatProcessor(w, h, bpix));
		membrane_channel_imp.setRoi(membrane_edges[fridx]);
		IJ.setForegroundColor(0, 255, 255);
		IJ.run(membrane_channel_imp, "Draw", "slice");
		membrane_channel_imp.setRoi(area_rois[fridx]);
		IJ.setForegroundColor(0, 255, 0);
		IJ.run(membrane_channel_imp, "Draw", "slice");
		for p in anchors:
			roi = PointRoi(p[0], p[1]);
			membrane_channel_imp.setRoi(roi);
			IJ.setForegroundColor(255, 0, 0);
			IJ.run(membrane_channel_imp, "Draw", "slice");
		roi = PointRoi(midpoint[0], midpoint[1]);
		membrane_channel_imp.setRoi(roi);
		IJ.setForegroundColor(255, 0, 0);
		IJ.run(membrane_channel_imp, "Draw", "slice");
		for p in fixed_anchors_list[fridx]:
			roi = PointRoi(p[0], p[1]);
			membrane_channel_imp.setRoi(roi);
			IJ.setForegroundColor(255, 255, 0);
			IJ.run(membrane_channel_imp, "Draw", "slice");
		# TODO: if coincidence between anchors and fixed anchors, or indeed with edge, consider 
		# rendering in secondary colours (magenta = (255,0,255), orange = (255,200,0), green = (0,255,0)?
	FileSaver(membrane_channel_imp).saveAsTiffStack(os.path.join(params.output_path, "membrane identification check.tif"));
	IJ.setForegroundColor(255, 255, 255);

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