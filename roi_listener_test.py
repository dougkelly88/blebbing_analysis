# imports
import math
from ij import IJ, ImagePlus;
from ij.gui import PolygonRoi, Roi, RoiListener
from ij.process import FloatPolygon
from ij.plugin import Straightener, Duplicator, ImageCalculator
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.gui import WaitForUserDialog


class MyRoiListener(RoiListener):
	def __init__(self):
		IJ.log("Listener started");
	
	def roiModified(self, imp, id):
		roi = imp.getRoi();
		print("roi length = " + str(roi.getLength()));

imp = IJ.openImage("C:\\Users\\dougk\\Desktop\\MAX_C2-Complex data.tif");
imp.show();
roi = Roi(100, 100, 100, 100);
imp.setRoi(roi);
roi.addRoiListener(MyRoiListener());