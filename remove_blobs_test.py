# @ImagePlus imp
# python (jython) imports
import os, sys, math
from datetime import datetime

# java imports - aim to have UI components entirely swing, listeners and layouts awt
from java.awt import Dimension, GridBagLayout, GridBagConstraints, GridLayout
import javax.swing as swing
import javax.swing.table.TableModel

# imagej imports
from ij import IJ, WindowManager, ImagePlus
from ij.gui import Roi, PointRoi, PolygonRoi, GenericDialog, WaitForUserDialog
from ij.io import OpenDialog, DirectoryChooser
from ij.plugin import ChannelSplitter
from ij.process import FloatPolygon
from loci.plugins import BF as bf
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable

rt = ResultsTable();
mxsz = imp.width * imp.height;
#pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER | ParticleAnalyzer.DOES_STACKS, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE, rt, 0, mx);
pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE, rt, 0, mxsz);

roim = RoiManager();
for idx in range(1, imp.getImageStackSize()+1):
	roim.reset();
	rt.reset();
	imp.setPosition(idx);
	pa.analyze(imp);
	#rt_slices = [int(r) for r in rt.getColumn(rt.getColumnIndex("Slice")).tolist()]
	rt_areas = rt.getColumn(rt.getColumnIndex("Area")).tolist();
	#rois = roim.getRoisAsArray().tolist();
	mx_ind = rt_areas.index(max(rt_areas))
	indices_to_remove = [a for a in range(0,len(rt_areas)) if a != mx_ind]
	print(indices_to_remove);
	for rem_idx in indices_to_remove:
		roim.select(imp, rem_idx);
		roim.runCommand(imp, "Fill");
roim.close();

#for idx in range(1, imp.getImageStackSize()+1):##
	#print(idx)
	#print("rt slices =  idx")
	#print(rt_slices[rt_slices == idx])
	#print("areas for that slice")
	#print(rt_areas[rt_slices == idx])
#print(rt_areas);	
#slice_areas = [a for a, b in zip(rt_areas, rt_slices) if b==idx];

