# @ImagePlus imp
import math
from ij import IJ, ImagePlus;
from ij.gui import PolygonRoi, Roi
from ij.process import FloatPolygon
from ij.plugin import Straightener
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.gui import WaitForUserDialog

membrane_edge = imp.getRoi();
poly = membrane_edge.getInterpolatedPolygon(0.25, False);
ys = [y for y in poly.ypoints];
xs = [x for x in poly.xpoints];
#ys = [int(round(y)) for y in poly.ypoints];
#xs = [int(round(x)) for x in poly.xpoints];
print("len xs");
print(len(xs));
area_roi = PolygonRoi(xs, ys, Roi.POLYGON);
imp2 = IJ.createImage("binary", imp.getWidth(), imp.getHeight(), 1, 8);
imp2.setRoi(area_roi);
mask = imp2.createRoiMask();
mskimp = ImagePlus("mask", mask);
IJ.run(mskimp, "Invert", "");

# get the largest area as the main part of the bleb
ip = mskimp.getProcessor();
mskimp.setProcessor(ip.resize(3 * imp.getWidth()));
IJ.run(mskimp, "Erode", "");
IJ.run(mskimp, "Erode", "");

rt = ResultsTable();
mxsz = mskimp.width * mskimp.height;
pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE, rt, 0, mxsz);

roim = RoiManager();

try:
	roim.reset();
except:
	print(roim);
pa.analyze(mskimp);
rt_areas = rt.getColumn(rt.getColumnIndex("Area")).tolist();
mx_ind = rt_areas.index(max(rt_areas))
indices_to_remove = [a for a in range(0,len(rt_areas)) if a != mx_ind]
for rem_idx in indices_to_remove:
	roim.select(mskimp, rem_idx);
	roim.runCommand(mskimp, "Fill");
roim.reset();
roim.close();
rt.reset();

ip = mskimp.getProcessor();
IJ.run(mskimp, "Dilate", "");
IJ.run(mskimp, "Dilate", "");
mskimp.setProcessor(ip.resize(imp.getWidth()));
rt = ResultsTable();
mxsz = imp.width * imp.height;
roim = RoiManager();
pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE, rt, 0, mxsz);
pa.analyze(mskimp);
nroi = roim.getRoi(0);
mskimp.setRoi(nroi);
roim.reset();
#roim.close();
imp.setRoi(membrane_edge);
print(nroi.getStatistics().area)

anchorline_roi = PolygonRoi([float(poly.xpoints[0]), float(poly.xpoints[-1])], [float(poly.ypoints[0]), float(poly.ypoints[-1])], Roi.POLYLINE);
anchorline_poly = anchorline_roi.getInterpolatedPolygon(0.25, False);
anchorline_roi = PolygonRoi([x for x in anchorline_poly.xpoints], [y for y in anchorline_poly.ypoints], Roi.POLYLINE);
# dilate msk by 1, get nnroi from particle analyzer and check for both contained in membrane edge and contained in nnroi...
mskimp.killRoi();
IJ.run(mskimp, "Dilate", "");
pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE, rt, 0, mxsz);
pa.analyze(mskimp);
nnroi = roim.getRoi(0);
imp.setRoi(nnroi);

#anchorline_xys = [(int(round(x)), int(round(y))) for x, y in zip(anchorline_poly.xpoints, anchorline_poly.ypoints)];
anchorline_xys = [(x, y) for x, y in zip(anchorline_poly.xpoints, anchorline_poly.ypoints)];
#print("anchor line:");
#print(anchorline_xys);
#membrane_xys = [(int(round(x)), int(round(y))) for x, y in zip(poly.xpoints, poly.ypoints)];
membrane_xys = [(x, y) for x, y in zip(poly.xpoints, poly.ypoints)];
#print("membrane: ");
#print(membrane_xys);
nnroi_xys = [(pt.x, pt.y) for pt in nnroi.getContainedPoints().tolist()]
#print("bleb roi:");
#print(nnroi_xys);

anchor_crossing_point_idx = set();
for idx, pt2 in enumerate(membrane_xys):
	for pt in anchorline_xys:
		if (abs(pt[0] - pt2[0]) < 0.5) and (abs(pt[1] - pt2[1]) < 0.5) and ((int(round(pt[0])), int(round(pt[1]))) in nnroi_xys):
			anchor_crossing_point_idx.add(idx);
	#if pt in nnroi_xys and pt in membrane_xys:
	#if pt in membrane_xys:
		#anchor_crossing_point_idx.append(idx);
anchor_crossing_point_idx = list(anchor_crossing_point_idx);
anchor_crossing_point_idx.sort();
#print(anchor_crossing_point_idx);
anchor_crossing_point_idx = [anchor_crossing_point_idx[0], anchor_crossing_point_idx[-1]];
print("Anchor crossing points:");
for idx in anchor_crossing_point_idx:
	print('' + str(idx) + ': ' + str(membrane_xys[idx]))
if anchor_crossing_point_idx[0] == 0 and anchor_crossing_point_idx[1] == len(membrane_xys):
	print("area is just nroi.getStatistics().area");
rotangle = membrane_edge.getAngle(int(xs[0]), int(ys[0]), int(xs[-1]), int(ys[-1])) / 180 * math.pi;
#rotangle = membrane_edge.getAngle(xs[0], ys[0], xs[-1], ys[-1]) / 180 * math.pi;
imp.setRoi(membrane_edge);
roim.close();

# check correct rotation - i.e. such that connecting line is horizontal...
print(rotangle)
#rotRoi = PolygonRoi([(x * math.cos(rotangle) - y * math.sin(rotangle)) for (x, y) in zip(poly.xpoints, poly.ypoints)], 
#					[(x * math.sin(rotangle) + y * math.cos(rotangle)) for (x, y) in zip(poly.xpoints, poly.ypoints)],
#					Roi.POLYGON);
rotRoi = PolygonRoi([(x * math.cos(rotangle) - y * math.sin(rotangle)) for (x, y) in membrane_xys], 
					[(x * math.sin(rotangle) + y * math.cos(rotangle)) for (x, y) in membrane_xys],
					Roi.POLYGON);
WaitForUserDialog("WAIT").show();
imp.setRoi(rotRoi); 
rotX = [(x * math.cos(rotangle) - y * math.sin(rotangle)) for (x, y) in membrane_xys];
rotY = [(x * math.sin(rotangle) + y * math.cos(rotangle)) for (x, y) in membrane_xys];
#print([(x, y) for x, y in zip(rotX, rotY)]);
WaitForUserDialog("WAIT").show();
imp.setRoi(membrane_edge);

# work back along the membrane to find the point of inflection...
#print("len rotY = " + str(len(rotY)));
seg1X = rotX[:anchor_crossing_point_idx[0]];
seg1Y = rotY[:anchor_crossing_point_idx[0]];
seg1_roi = PolygonRoi(seg1X, seg1Y, Roi.POLYLINE);
imp.setRoi(seg1_roi);
#print("seg1Y:");
#print((seg1Y));
WaitForUserDialog("WAIT").show();
#print("seg1:");
#print((seg1));
idx1 = seg1Y.index(max(seg1Y));
#for idx1 in range(anchor_crossing_point_idx[0], 3, -1):
#	dY = float(rotY[idx1] - rotY[idx1 - 3]);
#	dX = float(rotX[idx1] - rotX[idx1 - 3]);
#	if dY/dX < 0:
#		 break;
#print("identified point 1 index: " + str(idx1));
# work forward along the membrane to find the point of inflection...
seg2X = rotX[anchor_crossing_point_idx[1]:]
seg2Y = rotY[anchor_crossing_point_idx[1]:]
seg2_roi = PolygonRoi(seg2X, seg2Y, Roi.POLYLINE);
imp.setRoi(seg2_roi);
#print("seg2Y:");
#print((seg2Y));
WaitForUserDialog("WAIT").show();
idx2 = seg2Y.index(max(seg2Y));
imp.setRoi(membrane_edge);
#for idx2 in range(anchor_crossing_point_idx[1], len(rotY)-3, 3):
#	dY = float(rotY[idx2 + 3] - rotY[idx2]);
#	dY = float(rotX[idx2 + 3] - rotX[idx2]);
#	if dY/dX > 0:
#		 break;
print("identified point 2 index: " + str(idx2));
xs = [pt[0] for pt in membrane_xys];
xs = xs[idx1:idx2+1+anchor_crossing_point_idx[1]];
ys = [pt[1] for pt in membrane_xys];
ys = ys[idx1:idx2+1+anchor_crossing_point_idx[1]];
finalRoi = PolygonRoi(xs, ys, Roi.POLYGON);
WaitForUserDialog("WAIT").show();
imp.setRoi(finalRoi);
WaitForUserDialog("WAIT").show();
imp.setRoi(membrane_edge);
