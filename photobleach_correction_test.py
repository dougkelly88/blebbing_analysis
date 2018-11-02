from ij import IJ
from ij.process import ImageProcessor

class Parameters :

	def __init__(self, photobleaching_correction = True):
		self.photobleaching_correction = True;

def apply_photobleach_correction_framewise(params, actin_channel_imp, membrane_channel_imp, t0_value=None):
	"""if toggled on, scale the current frame in a stack to have the same mean intensity as the first frame within an ROI"""
	if not params.photobleaching_correction:
		return actin_channel_imp, None;
	else:
		print("T = " + str(actin_channel_imp.getT()));
		if (t0_value is None) or (actin_channel_imp.getT() == 1):
			IJ.run(membrane_channel_imp, "Create Selection", "");
			roi = membrane_channel_imp.getRoi();
			actin_channel_imp.setRoi(roi);
			t0_value = actin_channel_imp.getStatistics().mean;
		else:
			IJ.run(membrane_channel_imp, "Create Selection", "");
			roi = membrane_channel_imp.getRoi();
			actin_channel_imp.setRoi(roi);
			mean_value = actin_channel_imp.getStatistics().mean;
			factor = t0_value/mean_value;
			IJ.run(actin_channel_imp, "Multiply...", "value=" + str(factor) + " slice");
		return actin_channel_imp, t0_value;

nx = 3;
ny = 3;
nz = 1;
nc = 1;
nt = 5;
bitdepth = 8;
seg_imp = IJ.createHyperStack("segmentation image", nx, ny, nc, nz, nt, bitdepth);
intensity_imp = IJ.createHyperStack("intensity image", nx, ny, nc, nz, nt, bitdepth);
for idx in range(0, nt):
	seg_imp.setPosition(idx + 1);
	seg_pix = seg_imp.getProcessor().getPixels();
	for pixidx in range(len(seg_pix)):
		seg_pix[pixidx] = 1;
	intensity_imp.setPosition(idx + 1);
	int_pix = intensity_imp.getProcessor().getPixels();
	for pixidx in range(len(int_pix)):
		int_pix[pixidx] = nt - idx;

ip = seg_imp.getProcessor();
ip.setThreshold(8, 8, ImageProcessor.NO_LUT_UPDATE);
#IJ.run(seg_imp, "Convert to Mask", "")
IJ.run(seg_imp, "Make Binary", "method=Default background=Dark calculate");
#	intensity_imp.show();
#	seg_imp.show();
	
# apply photobleaching correction: 
t0_mean = None;
params = []
params = Parameters();
for idx in range(0, nt):
	intensity_imp.setPosition(idx+1);
	seg_imp.setPosition(idx+1);
	intensity_imp, t0_mean = apply_photobleach_correction_framewise(params, intensity_imp, seg_imp, t0_mean);
intensity_imp.show();