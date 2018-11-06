# imports
import math
from ij import IJ, ImagePlus, ImageListener;
from ij.gui import PolygonRoi, Roi
from ij.process import FloatPolygon
from ij.plugin import Straightener, Duplicator, ImageCalculator
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.gui import WaitForUserDialog


class MyImageListener(ImageListener):
	def __init__(self, extra_list):
		IJ.log("Listener started");
		self.last_frame = 1;
	
#	def imageUpdated(self, imp):
#		print("frame = " + str(imp.getT()));

	def imageUpdated(self, imp):
		frame = imp.getT();
		print("last frame = " + str(self.last_frame));
		print("frame = " + str(frame));
		extra_list.append(extra_list[-1] + 1);
		self.last_frame = frame;
		#print("Extra list = " + str(extra_list));

	def imageOpened(self, imp):
		print("Image opened");

	def imageClosed(self, imp):
		print("image closed");
		imp.removeImageListener(self);

imp = IJ.openImage("C:\\Users\\dougk\\Desktop\\C2-Complex data.tif");
imp.show();
extra_string = "hello ";
extra_list = [0];
imp.addImageListener(MyImageListener(extra_list));
