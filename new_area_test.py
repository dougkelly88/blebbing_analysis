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
#IJ.setTool("point");
#WaitForUserDialog("Choose midpoint anchor...").show();
#mpt_roi = imp.getRoi();
#p = mpt_roi.getContainedPoints();
#midpoint_anchor = (p[0].x, p[0].y);
#print("midpoint anchor = " + str(midpoint_anchor));
#IJ.setTool("freeline");

poly = membrane_edge.getPolygon();
xs = [x for x in poly.xpoints];
ys = [y for y in poly.ypoints];
#anchor_line_roi = PolygonRoi([xs[0], xs[-1]], 
#									[ys[0], ys[-1]], 
#									Roi.POLYLINE);
#imp.setRoi(anchor_line_roi);
#WaitForUserDialog("WAIT").show();
#anchor_line_midpoint = (xs[0] + int(round((xs[-1] - xs[0])/2)), ys[0] + int(round((ys[-1] - ys[0])/2)));
#print("anchor line midpoint = " + str(anchor_line_midpoint));
#orientation_line_roi = PolygonRoi([midpoint_anchor[0], anchor_line_midpoint[0]], 
#									[midpoint_anchor[1], anchor_line_midpoint[1]], 
#									Roi.POLYLINE);									
#imp.setRoi(orientation_line_roi);
#WaitForUserDialog("WAIT").show();
#orientation_angle = orientation_line_roi.getAngle(midpoint_anchor[0], midpoint_anchor[1], 
#													anchor_line_midpoint[0], anchor_line_midpoint[1]);
#orientation_angle = 180 * math.atan2(anchor_line_midpoint[1] - midpoint_anchor[1], 
#								anchor_line_midpoint[0] - midpoint_anchor[0]) / math.pi;
#print("Orientation: " + str(orientation_angle));
rotangle = membrane_edge.getAngle(xs[0], ys[0], xs[-1], ys[-1]) / 180 * math.pi;
meanX = sum(xs)/len(xs);
meanY = sum(ys)/len(ys);
rotX = [((x - meanX) * math.cos(rotangle) - (y - meanY) * math.sin(rotangle)) + meanX for x, y in zip(xs, ys)];
rotY = [((x - meanX) * math.sin(rotangle) + (y - meanY) * math.cos(rotangle)) + meanY for x, y in zip(xs, ys)];
rotY = [(x * math.sin(rotangle) + y * math.cos(rotangle)) for x, y in zip(xs, ys)];
meanRotY = sum(rotY)/len(rotY);
#rotRoi = PolygonRoi(rotX, rotY, Roi.POLYLINE);
#imp.setRoi(rotRoi);
#WaitForUserDialog("Showing rotated ROI?").show();
seg1 = rotY[:int(round(len(rotY)/2))];
seg1.reverse();
seg2 = rotY[int(round(len(rotY)/2)):];
#if orientation_angle > 0:
#print("mean rotY > baseline rotY?");
#print(meanRotY > rotY[0]);
if meanRotY > rotY[0]:
	idx1 = len(seg1) - seg1.index(min(seg1));
	idx2 = int(round(len(rotY)/2)) + seg2.index(min(seg2));
else:
	idx1 = len(seg1) - seg1.index(max(seg1));
	idx2 = int(round(len(rotY)/2)) + seg2.index(max(seg2));
#else:
#	idx1 = len(seg1) - seg1.index(min(seg1));
#	idx2 = int(round(len(rotY)/2)) + seg2.index(min(seg2));
area_poly_xs = xs[idx1:idx2+1];
area_poly_ys = ys[idx1:idx2+1];
area_roi = PolygonRoi(area_poly_xs, area_poly_ys, Roi.POLYGON);
area = area_roi.getStatistics().area;
imp.setRoi(area_roi);
print(area);
WaitForUserDialog("WAIT").show();
imp.setRoi(membrane_edge);
