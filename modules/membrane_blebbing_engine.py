# functions for handling image quantification for membrane blebbing analysis
#
# D. J. Kelly, 2018-10-15, douglas.kelly@riken.jp

# imports
import math
from ij import IJ;
from ij.gui import PolygonRoi, Roi
from ij.process import FloatPolygon
from ij.plugin import Straightener
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable

# convert the membrane identification channel into binary for segmentation
def make_and_clean_binary(imp, threshold_method):
	IJ.run(imp, "Make Binary", "method=" + threshold_method + " background=Dark calculate");
	IJ.run(imp, "Fill Holes", "stack");
	IJ.run(imp, "Open", "stack");
	IJ.run(imp, "Close-", "stack");
	keep_largest_blob(imp);
	return imp;

# move user-defined anchor points onto automatically-segmented membrane
# can do this more efficiently but probably not worth optimising...
# use set to implicitly avoid degeneracy, convert back to list. again, 
# inefficient, but OK for small set
def fix_anchors_to_membrane(anchors_list, membrane_roi):
	outline = membrane_roi.getPolygon();
	fixed_anchors_set = set();
	for anchor_idx, anchor in enumerate(anchors_list):
		fixed_anchor = anchor;
		last_dsq = 100000;
		for (x,y) in zip(outline.xpoints,outline.ypoints):
			d2 = math.pow((x - anchor[0]), 2) + math.pow((y - anchor[1]), 2);
			if d2 < last_dsq:
				last_dsq = d2;
				fixed_anchor = (x, y);
		fixed_anchors_set.add(fixed_anchor);
	if (len(fixed_anchors_set) < (anchor_idx+1)):
		raise ValueError('degeneracy between anchor points!');
	sortlist = list(fixed_anchors_set);
	vec_ls = [vector_length((0, IJ.getImage().getHeight()), v) for v in sortlist];
	return [sortlist[vec_ls.index(min(vec_ls))], sortlist[vec_ls.index(max(vec_ls))]];

# figure out which edge of the roi is the membrane, since IJ might start the roi
# from anywhere along the perimeter w.r.t. the user defined anchors
def get_membrane_edge(roi, fixed_anchors, fixed_midpoint):
	poly = roi.getPolygon();
	started = False;
	e1 = FloatPolygon();
	e2 = FloatPolygon();
	term_index_1 = [(x, y) for x, y in zip(poly.xpoints,poly.ypoints)].index( fixed_anchors[0]);
	term_index_2 = [(x, y) for x, y in zip(poly.xpoints,poly.ypoints)].index( fixed_anchors[1]);
	
	for idx in range(min(term_index_1, term_index_2), max(term_index_1, term_index_2)+1):
		e1.addPoint(poly.xpoints[idx], poly.ypoints[idx]);
	for idx in range(max(term_index_1, term_index_2), poly.npoints + min(term_index_1, term_index_2) + 1):
		e2.addPoint(poly.xpoints[idx % (poly.npoints)], poly.ypoints[idx % (poly.npoints)]);

	anchors_midpoint = (fixed_anchors[1][0] - fixed_anchors[0][0], 
						fixed_anchors[1][1] - fixed_anchors[0][1]);
	e1_mean = (sum(e1.xpoints)/e1.npoints, sum(e1.ypoints)/e1.npoints);
	e2_mean = (sum(e2.xpoints)/e2.npoints, sum(e2.ypoints)/e2.npoints);
	
	theta_e1 = angle_between_vecs(fixed_anchors[0], fixed_anchors[1], fixed_anchors[0], e1_mean);
	theta_e2 = angle_between_vecs(fixed_anchors[0], fixed_anchors[1], fixed_anchors[0], e2_mean);
	sign = lambda x: (1, -1)[x < 0]
	if sign(theta_e1) is not sign(theta_e2):
		theta_midpoint = angle_between_vecs(fixed_anchors[0], fixed_anchors[1], fixed_anchors[0], fixed_midpoint);
		use_edge = (e1, e2)[sign(theta_midpoint) == sign(theta_e2)];
	else:
		use_edge = (e1, e2)[vector_length(anchors_midpoint, e1_mean) < vector_length(anchors_midpoint, e2_mean)]
	return 	PolygonRoi(use_edge, Roi.POLYLINE);

# return angle between two vectors
def angle_between_vecs(u_start, u_end, v_start, v_end):
    u = (u_end[0] - u_start[0], u_end[1] - u_start[1]);
    v = (v_end[0] - v_start[0], v_end[1] - v_start[1]);
    return math.atan2(v[1], v[0]) - math.atan2(u[1], u[0]);

# return vector length
def vector_length(start, end):
	return math.sqrt(math.pow((start[0] - end[0]),2) + math.pow((start[1] - end[1]),2));

# generate arrays of points along the membrane that are separated by path length l - currently in pixels
def generate_l_spaced_points(roi, l): 
	poly = roi.getInterpolatedPolygon();
	p1, p2, p3 = ([] for i in range(3));
	for idx,(x,y) in enumerate(zip(poly.xpoints,poly.ypoints)):
		if ((idx > 0) and (idx < poly.npoints - 1)): # ignore first and last points: by definition these won't have anything useful on either side
			# look backwards and calculate pathlength at successive points
			db = 0;
			iidx = idx - 1;
			while ((iidx >= 0) and (db < l)):
				dbnew = db + math.sqrt(math.pow((poly.xpoints[iidx] - poly.xpoints[iidx+1]), 2)  + 
									math.pow((poly.ypoints[iidx] - poly.ypoints[iidx+1]), 2));
				if (dbnew >= l):
					xx = poly.xpoints[iidx+1] + ((l - db)/(dbnew - db))*(poly.xpoints[iidx] - poly.xpoints[iidx+1]);
					yy = poly.ypoints[iidx+1] - ((l - db)/(dbnew - db))*(poly.ypoints[iidx] - poly.ypoints[iidx+1]);
					dbnew = db + math.sqrt(math.pow(((l - db)/(dbnew - db))*(poly.xpoints[iidx] - poly.xpoints[iidx+1]), 2)  + 
									math.pow(((l - db)/(dbnew - db))*(poly.ypoints[iidx] - poly.ypoints[iidx+1]), 2));
				else:
					iidx = iidx-1;
				db = dbnew;
			if (db == l):
				pp1 = (xx, yy);
				pp2 = (x, y);
				# then look forwards ONLY IF backwards search was successful...
				iidx = idx + 1;
				df = 0;
				while ((iidx < poly.npoints - 1) and (df < l)):
					dfnew = df + math.sqrt(math.pow((poly.xpoints[iidx] - poly.xpoints[iidx-1]), 2)  + 
									math.pow((poly.ypoints[iidx] - poly.ypoints[iidx-1]), 2));
					if (dfnew >= l):
						xx = poly.xpoints[iidx-1] + ((l - df)/(dfnew - df))*(poly.xpoints[iidx] - poly.xpoints[iidx-1]);
						yy = poly.ypoints[iidx-1] + ((l - df)/(dfnew - df))*(poly.ypoints[iidx] - poly.ypoints[iidx-1]);
						dfnew = df + math.sqrt(math.pow(((l - df)/(dfnew - df))*(poly.xpoints[iidx] - poly.xpoints[iidx-1]), 2)  + 
										math.pow(((l - df)/(dfnew - df))*(poly.ypoints[iidx] - poly.ypoints[iidx-1]), 2));
					else:
						iidx = iidx+1;
					df = dfnew;
				if (df == l):
					p1.append(pp1);
					p2.append(pp2);
					p3.append((xx, yy));
	return (p1, p2, p3);

# generate a line profile of local curvatures using three-point method and SSS theorem
# (see http://mathworld.wolfram.com/SSSTheorem.html)
def calculate_curvature_profile(curv_points, roi, remove_negative_curvatures):
	poly = roi.getInterpolatedPolygon();
	curvature_profile = [((x,y),0) for (x,y) in zip(poly.xpoints, poly.ypoints)];
	pos = [p for (p, c) in curvature_profile]
	for (cp, p1, p2) in zip(curv_points[1], curv_points[0], curv_points[2]):
		a = math.sqrt( math.pow((cp[0] - p1[0]), 2) + math.pow((cp[1] - p1[1]), 2) );
		b = math.sqrt( math.pow((cp[0] - p2[0]), 2) + math.pow((cp[1] - p2[1]), 2) );
		c = math.sqrt( math.pow((p2[0] - p1[0]), 2) + math.pow((p2[1] - p1[1]), 2) );
		s = 0.5 * (a + b + c);
		K = math.sqrt(s * (s - a) * (s - b) * (s - c));
		if (K == 0):
			curv = 0;
		else:
			R = (a * b * c)/(4 * K);
			c_1 = tuple(r1-rc for r1, rc in zip(p1, cp));
			c_2 = tuple(r2-rc for r2, rc in zip(p2, cp));
			sign = round((c_1[0]*c_2[1] - c_2[0]*c_1[1]) / (abs(c_1[0]*c_2[1] - c_2[0]*c_1[1]) + 1e-10));
			curv = sign * 1/R;
		if (remove_negative_curvatures and (curv < 0)):
			curv = 0;
		curvature_profile[pos.index(cp)] = (cp, curv);
	return curvature_profile;

# remove all blobs other than the largest by area
def keep_largest_blob(imp):
	rt = ResultsTable();
	mxsz = imp.width * imp.height;
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE, rt, 0, mxsz);

	roim = RoiManager();
	for idx in range(1, imp.getImageStackSize()+1):
		roim.reset();
		rt.reset();
		imp.setPosition(idx);
		pa.analyze(imp);
		rt_areas = rt.getColumn(rt.getColumnIndex("Area")).tolist();
		mx_ind = rt_areas.index(max(rt_areas))
		indices_to_remove = [a for a in range(0,len(rt_areas)) if a != mx_ind]
		for rem_idx in indices_to_remove:
			roim.select(imp, rem_idx);
			roim.runCommand(imp, "Fill");
	roim.reset();
	roim.close();

# return a line profile taking the maximum value over n pixels perpendicular to roi line
def maximum_line_profile(imp, roi, pixel_width):
	imp.setRoi(roi);
	ip = Straightener().straightenLine(imp, pixel_width);
	width = ip.getWidth();
	height = ip.getHeight();
	max_profile = [];
	poly = roi.getInterpolatedPolygon();
	for idx,(x,y) in enumerate(zip(poly.xpoints,poly.ypoints)):
		pix = ip.getLine(idx, 0, idx, height);
		max_profile.append(((x, y), max(pix)));
	return max_profile;