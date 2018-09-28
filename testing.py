# @ImagePlus imp
# python (jython) imports
import os, sys, math
from datetime import datetime

# java imports - aim to have UI components entirely swing, listeners and layouts awt
from java.awt import Dimension, GridBagLayout, GridBagConstraints, GridLayout
import javax.swing as swing
import javax.swing.table.TableModel

# imagej imports
from ij import IJ, WindowManager
from ij.gui import Roi, PointRoi, PolygonRoi, GenericDialog, WaitForUserDialog
from ij.io import OpenDialog, DirectoryChooser
from ij.plugin import ChannelSplitter
from ij.process import FloatPolygon
from loci.plugins import BF as bf

IJ.run(imp, "Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
roi = imp.getRoi();
poly = roi.getFloatPolygon();
p1, p2, p3 = ([] for i in range(3))
l = 5.0; # pixel length scale

for idx,(x,y) in enumerate(zip(poly.xpoints,poly.ypoints)):
	if ((idx > 0) and (idx < poly.npoints - 1)): # ignore first and last points: by definition these won't have anything useful on either side
		# look backwards and calculate pathlength at successive points
		db = 0;
		iidx = idx - 1;
		#print("--------");
		#print("test centre: " + str((x,y)));
		while ((iidx >= 0) and (db < l)):
			dbnew = db + math.sqrt(math.pow((poly.xpoints[iidx] - poly.xpoints[iidx+1]), 2)  + 
								math.pow((poly.ypoints[iidx] - poly.ypoints[iidx+1]), 2));
			#print("test point = " + str((poly.xpoints[iidx], poly.ypoints[iidx])));
			if (dbnew >= l):
				xx = poly.xpoints[iidx+1] + ((l - db)/(dbnew - db))*(poly.xpoints[iidx] - poly.xpoints[iidx+1]);
				yy = poly.ypoints[iidx+1] - ((l - db)/(dbnew - db))*(poly.ypoints[iidx] - poly.ypoints[iidx+1]);
				#print("interpolated point = " + str((xx, yy)));
				dbnew = db + math.sqrt(math.pow(((l - db)/(dbnew - db))*(poly.xpoints[iidx] - poly.xpoints[iidx+1]), 2)  + 
								math.pow(((l - db)/(dbnew - db))*(poly.ypoints[iidx] - poly.ypoints[iidx+1]), 2));
			else:
				iidx = iidx-1;
			db = dbnew;
			#print("distance to backward test point = " + str(db));
			#gd = GenericDialog("Continue?");
			#gd.showDialog();
			#if (gd.wasCanceled()):
			#	raise ValueError('stop!');
		if (db == l):
			pp1 = (xx, yy);
			pp2 = (x, y);
			# then look forwards ONLY IF backwards search was successful...
			iidx = idx + 1;
			df = 0;
			while ((iidx < poly.npoints - 1) and (df < l)):
				dfnew = df + math.sqrt(math.pow((poly.xpoints[iidx] - poly.xpoints[iidx-1]), 2)  + 
								math.pow((poly.ypoints[iidx] - poly.ypoints[iidx-1]), 2));
				#print("test point = " + str((poly.xpoints[iidx], poly.ypoints[iidx])));
				if (dfnew >= l):
					xx = poly.xpoints[iidx-1] + ((l - df)/(dfnew - df))*(poly.xpoints[iidx] - poly.xpoints[iidx-1]);
					yy = poly.ypoints[iidx-1] + ((l - df)/(dfnew - df))*(poly.ypoints[iidx] - poly.ypoints[iidx-1]);
					#print("interpolated point = " + str((xx, yy)));
					dfnew = df + math.sqrt(math.pow(((l - df)/(dfnew - df))*(poly.xpoints[iidx] - poly.xpoints[iidx-1]), 2)  + 
									math.pow(((l - df)/(dfnew - df))*(poly.ypoints[iidx] - poly.ypoints[iidx-1]), 2));
				else:
					iidx = iidx+1;
				df = dfnew;
				#print("distance to forward test point = " + str(df));
				#gd = GenericDialog("Continue?");
				#gd.showDialog();
				#if (gd.wasCanceled()):
				#	raise ValueError('stop!');
			if (df == l):
				p1.append(pp1);
				p2.append(pp2);
				p3.append((xx, yy));
				#print("centre point=");
				#print(pp2);
				#print("backwards point=");
				#print(pp1);
				#print("path length between centre and backwards=");
				#print(db);
				#print("forwards point=");
				#print((xx, yy));
				#print("path length between centre and forwards=");
				#print(df);
print("centre points: " + str(p2[0:2]))
print("backwards points: " + str(p1[0:2]))
print("forwards points: " + str(p3[0:2]))