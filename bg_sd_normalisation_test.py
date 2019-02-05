from ij import IJ, ImagePlus
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable, Measurements

import math

def get_background_stdev(imp, mask_imp, midpoint_anchor, dilations=5):
	"""get intensity standard deviation from a background region identified by auto-thresholded image and position of midpoint anchor"""
	ip = mask_imp.getProcessor();
	dilations = 5;
	for d in range(dilations):
		ip.dilate();
	mask_imp.setProcessor(ip);
	IJ.run(mask_imp, "Create Selection", "");
	mask_imp.crop();
	ip = mask_imp.getProcessor();
	ip.invert();
	rt = ResultsTable();
	mxsz = mask_imp.getWidth() * mask_imp.getHeight();
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER | ParticleAnalyzer.SHOW_PROGRESS, ParticleAnalyzer.CENTROID, rt, 0, mxsz);
	roim = RoiManager();
	pa.analyze(mask_imp);
	ds_to_anchor = [math.sqrt((x - midpoint_anchor[0])**2 + (y - midpoint_anchor[1])**2) for x, y in zip(rt.getColumn(rt.getColumnIndex("X")).tolist(), rt.getColumn(rt.getColumnIndex("X")).tolist())];
	roi = roim.getRoi(ds_to_anchor.index(min(ds_to_anchor)));
	imp.setRoi(roi);
	std_dev = imp.getStatistics().stdDev;
	print("original image standard deviation = " + str(std_dev));
	roim.reset();
	roim.close();
	rt.reset();
	return std_dev, roi;


midpoint_anchor = (1, 1);

mask_path = "C:\\Users\\dougk\\Desktop\\output\\2019-01-29 17-10-59 output\\Varying curvature length parameter\\curvature_length_um = 1.0\\binary_membrane_stack.tif";
original_img_path = "D:\\source\\vascular_morphogenesis_ij\\blebbing analysis\\test data\\Simple data.tif";

mask_imp = IJ.openImage(mask_path);
mask_imp.show();
mask_imp.killRoi();
imp = IJ.openImage(original_img_path);
imp.show();

print(get_background_stdev(imp, mask_imp, midpoint_anchor));