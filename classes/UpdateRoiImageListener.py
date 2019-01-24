# class to support updating ROI from list upon change of frame
#
# D. J. Kelly, 2018-11-06, douglas.kelly@riken.jp

from ij import ImageListener

class UpdateRoiImageListener(ImageListener):
	"""class to support updating ROI from list upon change of frame"""
	def __init__(self, roi_list):
		self.last_frame = 1;
		self.roi_list = roi_list;
		print("UpdateRoiImageListener started");

	def imageUpdated(self, imp):
		roi = imp.getRoi();
		if roi is not None and not roi.isArea():
			self.roi_list[self.last_frame - 1] = roi;
		if imp.getNFrames()>imp.getNSlices():
			frame = imp.getT();
		else:
			frame = imp.getZ();
		self.last_frame = frame;
		imp.setRoi(self.roi_list[frame - 1]);

	def imageOpened(self, imp):
		print("UpdateRoiImageListener: image opened");
			
	def imageClosed(self, imp):
		print("UpdateRoiImageListener: image closed");
		imp.removeImageListener(self);

	def getRoiList(self):
		return self.roi_list;

## testing
#from ij import IJ
#from ij.gui import Roi, PolygonRoi

#imp = IJ.openImage("C:\\Users\\dougk\\Desktop\\C2-Complex data.tif");
#imp.show();
#rois = [];
#for fridx in range(imp.getNFrames()):
#	x0 = 2 * fridx + 3;
#	w = imp.getWidth() - 3 - 4 * fridx;
#	y0 = 2 * fridx + 3;
#	h = imp.getHeight() - 3 - 4 * fridx;
#	roi = PolygonRoi([x0, (x0+w)/2, (x0+w)], [y0, (y0+h)/2, (y0+h)], Roi.POLYLINE);
#	#roi = Roi(x0, y0, w, h) 
#	rois.append(roi);
##listener = UpdateRoiImageListener(rois)
##imp.addImageListener(listener);
#print([zip(poly.xpoints, poly.ypoints) for poly in [roi.getPolygon() for roi in rois]])