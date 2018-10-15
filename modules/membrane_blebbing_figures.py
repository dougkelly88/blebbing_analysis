# functions for handling output image generation for membrane blebbing analysis
#
# D. J. Kelly, 2018-10-15, douglas.kelly@riken.jp

# TODO: refactor as a class so that it is trivial to configure whether images should be
# shown globally, e.g. fig = FigClass(show_figs = False); fig.overlay_curvatures(...)
# also possibly easier in this case to avoid boilerplate in image saving...

# imports
from ij import IJ, ImagePlus, ImageStack
from ij.gui import Plot
from ij.plugin import RGBStackMerge
from ij.process import FloatProcessor
import membrane_blebbing_engine as mb;

# generate overlays to display curvature
def generate_curvature_overlays(curvature_profile, curvature_stack): 
	w = curvature_stack.getWidth();
	ip = FloatProcessor(w, curvature_stack.getHeight());
	pix = ip.getPixels();
	mx = max([c for ((x, y), c) in curvature_profile]);
	for ((x, y), c) in curvature_profile:
		pix[int(round(y)) * w + int(round(x))] = c;
	curvature_stack.addSlice(ip);
	return curvature_stack;

# overlay curvature pixels on membrane image
def overlay_curvatures(imp, curvature_stack, curvature_profiles, membrane_channel, limits, colormap_string):
	overlay_base_imp = imp.clone();
	overlay_imp = ImagePlus("Curvature stack", curvature_stack);
	IJ.run(overlay_imp, colormap_string, "");
	IJ.setMinAndMax(overlay_imp, int(limits[0]), int(limits[1]));
	raw_overlay = overlay_imp.clone();
	IJ.run(overlay_imp, "RGB Color", "");
	overlaid_stack = ImageStack(overlay_imp.width, overlay_imp.height);
	for fridx in range(1, curvature_stack.getSize()+1):
		raw_idx = overlay_base_imp.getStackIndex(membrane_channel, 1, fridx);
		ip = overlay_base_imp.getStack().getProcessor(raw_idx).convertToRGB();
		pix = overlay_imp.getStack().getProcessor(fridx).getPixels();
		base_pix = ip.getPixels();
		for ((x, y), c) in curvature_profiles[fridx-1]:
			if (c > 0):
				base_pix[int(round(y)) * imp.width + int(round(x))] = pix[int(round(y)) * imp.width + int(round(x))];
		overlaid_stack.addSlice(ip);
	return ImagePlus("Overlaid curvatures", overlaid_stack), raw_overlay;

# display unnormalised kymograph
def generate_plain_kymograph(data_to_plot, colormap_string, title_string):
	ip = FloatProcessor(len(data_to_plot), max([len(data) for data in data_to_plot]));
	pix = ip.getPixels();
	for idx, data in enumerate(data_to_plot):
		for yidx in range(0,len(data)):
			pix[yidx * len(data_to_plot) + idx] = data[yidx][1];
	imp = ImagePlus(title_string, ip);
	IJ.run(imp, colormap_string, "");
	imp.show();
	return imp;
	
# display one-channel kymograph with point furthest from the edges along the middle of the
# kymograph
def generate_kymograph(data_to_plot, colormap_string, title_string):
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

# generate intensity-weighted curvature image
def generate_intensity_weighted_curvature(curvature_overlay, curvature_profiles, intensity_channel_imp, colormap_string):
	overlay_imp = curvature_overlay.clone();
	base_imp = intensity_channel_imp.clone();
	IJ.run(overlay_imp, colormap_string, "");
	IJ.run(overlay_imp, "RGB Color", "");
	
	for idx, profile in enumerate(curvature_profiles):
		overlay_imp.setPosition(idx + 1);
		base_imp.setPosition(idx + 1);
		overlay_pix = overlay_imp.getProcessor().getPixels();
		base_pix = overlay_imp.getProcessor().getPixels();
		w = overlay_imp.getWidth();
		for ((x,y), c) in profile:
			xyidx = int(round(y)) * w + int(round(x));
			base_pix[xyidx] = base_pix[xyidx] * overlay_pix[xyidx];
	return base_imp;

# plot "bleb length"
def plot_bleb_length(membrane_edges):
	ts = [x for x in range(0, len(membrane_edges))];
	bleb_ls = [membrane_edge.getLength() for membrane_edge in membrane_edges];
	plt = Plot("Crude bleb perimeter length against time", "Time, frames", "Bleb length, um");
	plt.add("line", ts, bleb_ls)
	plt.show();
	plot_data = [(t, l) for t, l in zip(ts, bleb_ls)];
	return plt.getImagePlus(), plot_data

# merge two kymographs
def merge_kymographs(kym1_imp, kym2_imp):
	mrg_imp = RGBStackMerge().mergeChannels([kym1_imp, kym2_imp], True);
	mrg_imp.setTitle("Merged actin intensity and curvature kymograph");
	mrg_imp.show();
	return mrg_imp;