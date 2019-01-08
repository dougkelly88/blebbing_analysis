# @ImagePlus imp
import json, math
from ij.gui import PolygonRoi, Roi, WaitForUserDialog

def load_qcd_edges(input_file_path):
	"""load edges from JSON"""
	f = open(input_file_path, 'r');
	try:
		edges = json.loads(f.read());
	finally:
		f.close();
	membrane_edges = [];
	for edge in edges:
		print(edge);
		print("");
		xs = [pt[0] for pt in edge];
		print("xs");
		print(xs);
		ys = [pt[1] for pt in edge];
		print("ys");
		print(ys);
		membrane_edges.append(PolygonRoi(xs, ys, Roi.POLYLINE));
		print("polyxs = " + str(membrane_edges[-1].getFloatPolygon().xpoints));
	return membrane_edges;

def bleb_area(membrane_edge):
	"""calculate the area of the drawn bleb, accounting for membrane crossing line joining anchors"""
	if isinstance(membrane_edge, Roi):
		poly = membrane_edge.getFloatPolygon();
		xs = [x for x in poly.xpoints];
		ys = [y for y in poly.ypoints];
		rotangle = membrane_edge.getAngle(int(round(xs[0])), int(round(ys[0])), int(round(xs[-1])), int(round(ys[-1]))) / 180 * math.pi;
	else:
		xs = [pt[0] for pt in membrane_edge];
		ys = [pt[1] for pt in membrane_edge];
		dumRoi = PolygonRoi(xs, xs, Roi.POLYGON);
		rotangle = dumRoi.getAngle(int(round(xs[0])), int(round(ys[0])), int(round(xs[-1])), int(round(ys[-1]))) / 180 * math.pi;
	rotY = [(x * math.sin(rotangle) + y * math.cos(rotangle)) for x, y in zip(xs, ys)];
	meanRotY = sum(rotY)/len(rotY);
	seg1 = rotY[:int(round(len(rotY)/2))];
	seg1.reverse();
	seg2 = rotY[int(round(len(rotY)/2)):];
	if meanRotY > rotY[0]:
		idx1 = len(seg1) - seg1.index(min(seg1));
		idx2 = int(round(len(rotY)/2)) + seg2.index(min(seg2));
	else:
		idx1 = len(seg1) - seg1.index(max(seg1));
		idx2 = int(round(len(rotY)/2)) + seg2.index(max(seg2));
	area_poly_xs = xs[idx1:idx2+1];
	area_poly_ys = ys[idx1:idx2+1];
	len_roi = PolygonRoi(area_poly_xs, area_poly_ys, Roi.POLYLINE);
	length = len_roi.getLength();
	area_roi = PolygonRoi(area_poly_xs, area_poly_ys, Roi.POLYGON);
	area = area_roi.getStatistics().area;
	print(area);
	return length, area, area_roi;

def reflect_point_over_line(line_start, line_end, point_to_reflect):
    if (line_end[0]==line_start[0]): # i.e., vertical line
        return (-point_to_reflect[0], point_to_reflect[1]);
    m = (line_end[1] - line_start[1])/(line_end[0] - line_start[0]);
    c = (line_start[1] - m * line_start[0]);
    d = (point_to_reflect[0] + (point_to_reflect[1] - c)*m)/(1 + math.pow(m,2));
    reflected_point = (2*d - point_to_reflect[0],  2*d*m - point_to_reflect[1] + 2*c);
#    return reflected_point, m, c;
    return reflected_point;
	
def directed_bleb_area(membrane_edge, midpoint_anchor):
	"""calculate the area of the drawn bleb, accounting for membrane crossing line joining anchors"""
	if isinstance(membrane_edge, Roi):
		poly = membrane_edge.getFloatPolygon();
		xs = [x for x in poly.xpoints];
		ys = [y for y in poly.ypoints];
		rotangle = membrane_edge.getAngle(int(round(xs[0])), int(round(ys[0])), int(round(xs[-1])), int(round(ys[-1]))) / 180 * math.pi;
	else:
		xs = [pt[0] for pt in membrane_edge];
		ys = [pt[1] for pt in membrane_edge];
		dumRoi = PolygonRoi(xs, xs, Roi.POLYGON);
		rotangle = dumRoi.getAngle(int(round(xs[0])), int(round(ys[0])), int(round(xs[-1])), int(round(ys[-1]))) / 180 * math.pi;
	rotY = [(x * math.sin(rotangle) + y * math.cos(rotangle)) for x, y in zip(xs, ys)];
	rotYmpa = midpoint_anchor[0] * math.sin(rotangle) + midpoint_anchor[1] * math.cos(rotangle)
	meanRotY = sum(rotY)/len(rotY);
	seg1 = rotY[:int(round(len(rotY)/2))];
	seg1.reverse();
	seg2 = rotY[int(round(len(rotY)/2)):];
#	if meanRotY > rotY[0]:
	if rotYmpa > rotY[0]:
		idx1 = len(seg1) - seg1.index(min(seg1));
		idx2 = int(round(len(rotY)/2)) + seg2.index(min(seg2));
	else:
		idx1 = len(seg1) - seg1.index(max(seg1));
		idx2 = int(round(len(rotY)/2)) + seg2.index(max(seg2));
	area_poly_xs = xs[idx1:idx2+1];
	area_poly_ys = ys[idx1:idx2+1];
	len_roi = PolygonRoi(area_poly_xs, area_poly_ys, Roi.POLYLINE);
	length = len_roi.getLength();
	area_roi = PolygonRoi(area_poly_xs, area_poly_ys, Roi.POLYGON);
	area = area_roi.getStatistics().area;
	return length, area, area_roi;

input_edge_path = "D:\\analysis\\2018-12-28 Images for LKP presentation\\Marcksl1-EGFP plasmid injected into Tg(Fli1 Lifeact-mCherry) embryos 1\\A) Control\\user_defined_edges.json";
#input_edge_path = "D:\\analysis\\2018-12-28 Images for LKP presentation\\Marcksl1-EGFP plasmid injected into Tg(Fli1 Lifeact-mCherry) embryos 1\\C) Marksl1b OE - retracting bleb - inc before\\user_defined_edges.json";
edges = load_qcd_edges(input_edge_path);
for idx, edge in enumerate(edges):
	print(edge)
	l, a, aroi = directed_bleb_area(edge, (5,63));
	print(a);
	imp.setZ(idx+1);
	#imp.setRoi(aroi);
	imp.setRoi(edge);
	WaitForUserDialog("pause").show();

