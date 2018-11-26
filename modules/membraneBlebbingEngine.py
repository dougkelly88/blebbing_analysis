# functions for handling image quantification for membrane blebbing analysis
#
# D. J. Kelly, 2018-10-15, douglas.kelly@riken.jp

# imports
import math
from ij import IJ, ImagePlus;
from ij.gui import PolygonRoi, Roi
from ij.process import FloatPolygon
from ij.plugin import Straightener, Duplicator, ImageCalculator
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.gui import WaitForUserDialog

def make_and_clean_binary(imp, threshold_method):
	"""convert the membrane identification channel into binary for segmentation"""
	if "Local: " in threshold_method:
		dup = Duplicator();
		imp1 = dup.run(imp);
		imp2 =  dup.run(imp);
		imp.changes = False;
		imp.close();
		threshold_method = threshold_method.split("Local: ")[1];
		IJ.run(imp1, "8-bit", "");
		IJ.run(imp1, "Auto Local Threshold", "method=" + threshold_method + " radius=15 parameter_1=0 parameter_2=0 white stack");
		IJ.run(imp2, "Make Binary", "method=MinError background=Dark calculate");
		ic = ImageCalculator();
		imp = ic.run("AND create stack", imp1, imp2);
		IJ.run(imp, "Invert", "stack");
		IJ.run(imp, "Make Binary", "method=Default background=Default calculate");
	elif "Edge" in threshold_method:
		IJ.run(imp, "Find Edges", "stack");
		IJ.run(imp, "Make Binary", "method=Mean background=Dark calculate");	
	else:
		IJ.run(imp, "Make Binary", "method=" + threshold_method + " background=Dark calculate"); # "calculate" ensures that threshold is calculated image-wise

	IJ.run(imp, "Open", "stack");
	IJ.run(imp, "Close-", "stack");
	IJ.run(imp, "Close-", "stack");
	IJ.run(imp, "Open", "stack");
	IJ.run(imp, "Fill Holes", "stack");
	IJ.run(imp, "Erode", "stack");
	IJ.run(imp, "Erode", "stack");
	keep_largest_blob(imp);
	IJ.run(imp, "Dilate", "stack");
	IJ.run(imp, "Dilate", "stack");
	IJ.run(imp, "Open", "stack");
	if "Edge" in threshold_method:
		IJ.run(imp, "Erode", "stack");
		IJ.run(imp, "Erode", "stack");
	return imp;

def distances_from_point_along_path(path_roi, point):
	"""For each point along a line ROI, calculate the distance from a given point"""
	poly = path_roi.getPolygon();
	pts = zip(poly.xpoints, poly.ypoints)
	print(pts)
	distances = [0] * len(pts);
	try:
		pidx = pts.index(point);
	except ValueError:
		# point doesn't exist on membrane, find closest point on membrane
		dsq = [(math.pow((x - point[0]),2) + math.pow((y - point[1]),2)) for (x, y) in pts];
		pidx = dsq.index(min(dsq));
	for idx in range(pidx+1, poly.npoints):
		distances[idx] = distances[idx-1] + vector_length(pts[idx], pts[idx-1]);
	for idx in range(pidx-1, -1, -1):
		distances[idx] = distances[idx+1] + vector_length(pts[idx], pts[idx+1]);
	return distances;

#def check_anchor_distance_and_gradient(old_anchors_list, old_membrane_rois, new_anchors, new_roi, params):
#	"""check whether latest set of anchors has significant deviation in space and in gradient from last anchors"""
#	maximum_creep_speed_um_s = 0.2;
#	position_tolerance_microns = params.frame_interval * maximum_creep_speed_um_s;
#	gradient_angle_tolerance = 30.0/180 * math.pi;
#	scale_for_gradient_calc_um = 0.5;
	
#	for anc_idx in range(0,2):
#		dr = vector_length(old_anchors_list[-1][anc_idx], new_anchors[anc_idx]);
#		if (dr * params.pixel_physical_size) > position_tolerance_microns:
#			new_anchors[anc_idx] = old_anchors_list[-1][anc_idx];
#		else:

#	return new_anchors, new_roi;

def fix_anchors_to_membrane(anchors_list, membrane_roi, params):
	"""move user-defined anchor points onto automatically-segmented membrane. TODO: make more efficient?"""
	outline = membrane_roi.getPolygon();
	pts = zip(outline.xpoints, outline.ypoints);
	fixed_anchors_set = set();
	for anchor_idx, anchor in enumerate(anchors_list):
		fixed_anchor = anchor;
		dsq = [(math.pow((x - anchor[0]),2) + math.pow((y - anchor[1]),2)) for (x, y) in pts];
		pidx = dsq.index(min(dsq))
		fixed_anchors_set.add(pts[pidx]);
	if (len(fixed_anchors_set) < (anchor_idx+1)):
		raise ValueError('degeneracy between anchor points!');
	sortlist = list(fixed_anchors_set);
	vec_ls = [vector_length((0, IJ.getImage().getHeight()), v) for v in sortlist];
	return [sortlist[vec_ls.index(min(vec_ls))], sortlist[vec_ls.index(max(vec_ls))]];

def get_membrane_edge(roi, fixed_anchors, fixed_midpoint):
	"""figure out which edge of the roi is the membrane, since IJ might start the roi from anywhere along the perimeter w.r.t. the user defined anchors"""
	poly = roi.getPolygon();
	term_index_1 = [(x, y) for x, y in zip(poly.xpoints,poly.ypoints)].index(fixed_anchors[0]);
	term_index_2 = [(x, y) for x, y in zip(poly.xpoints,poly.ypoints)].index(fixed_anchors[1]);
	start_idx = min(term_index_1, term_index_2);
	end_idx = max(term_index_1, term_index_2);
	xs = [x for x in poly.xpoints];
	ys = [y for y in poly.ypoints];
	e1 = FloatPolygon(xs[start_idx:end_idx+1], ys[start_idx:end_idx+1]);
	e2 = FloatPolygon(list(reversed(xs[end_idx:] + xs[:start_idx+1])), 
						list(reversed(ys[end_idx:] + ys[:start_idx+1])));

	anchors_midpoint = (int(round(0.5 * (fixed_anchors[1][0] + fixed_anchors[0][0]))), 
						int(round(0.5 * (fixed_anchors[1][1] + fixed_anchors[0][1]))));
	e1_mean = (sum(e1.xpoints)/e1.npoints, sum(e1.ypoints)/e1.npoints);
	e2_mean = (sum(e2.xpoints)/e2.npoints, sum(e2.ypoints)/e2.npoints);
	
	theta_e1 = angle_between_vecs(fixed_anchors[0], fixed_anchors[1], fixed_anchors[0], e1_mean);
	theta_e2 = angle_between_vecs(fixed_anchors[0], fixed_anchors[1], fixed_anchors[0], e2_mean);
	sign = lambda x: (1, -1)[x < 0]
	if sign(theta_e1) is not sign(theta_e2):
		#print("using angle to decide on edge ID - vectors linking anchor1 and mean edge position lie on either side of anchor line");
		theta_midpoint = angle_between_vecs(fixed_anchors[0], fixed_anchors[1], fixed_anchors[0], fixed_midpoint);
		use_edge = e2 if (sign(theta_midpoint) == sign(theta_e2)) else e1;
	else:
		#print("using distance to decide on edge ID - vectors linking anchor1 and mean edge position lie on same side of anchor line");
		#print("position anchor midpoint = " + str(anchors_midpoint));
		#print("position e1 mean = " + str(e1_mean));
		#print("length anchor midpoint to e1 mean = " + str(vector_length(anchors_midpoint, e1_mean)));
		#print("position e2 mean = " + str(e2_mean));
		#print("length anchor midpoint to e2 mean = " + str(vector_length(anchors_midpoint, e2_mean)));
		use_edge = e1 if (vector_length(anchors_midpoint, e1_mean) > vector_length(anchors_midpoint, e2_mean)) else e2;
	return 	PolygonRoi(use_edge, Roi.POLYLINE);

def angle_between_vecs(u_start, u_end, v_start, v_end):
	"""return angle between two vectors"""
	u = (u_end[0] - u_start[0], u_end[1] - u_start[1]);
	v = (v_end[0] - v_start[0], v_end[1] - v_start[1]);
	return math.atan2(v[1], v[0]) - math.atan2(u[1], u[0]);

def vector_length(start, end):
	"""return vector length, given start and end points"""
	return math.sqrt(math.pow((start[0] - end[0]),2) + math.pow((start[1] - end[1]),2));

def generate_l_spaced_points(roi, l): 
	"""generate arrays of points along the membrane that are separated by path length l - currently in pixels"""
	poly = roi.getInterpolatedPolygon(1.0, True);
	p1, cp, p2 = ([] for i in range(3));
	for idx,(x,y) in enumerate(zip(poly.xpoints,poly.ypoints)):
		if ((idx > 0) and (idx < poly.npoints)): # ignore first and last points: by definition these won't have anything useful on either side
			# look backwards and calculate pathlength at successive points
			db = 0;
			iidx = idx - 1;
			while ((iidx >= 0) and (db < l)):
				dbnew = db + math.sqrt(math.pow((poly.xpoints[iidx] - poly.xpoints[iidx+1]), 2)  + 
									math.pow((poly.ypoints[iidx] - poly.ypoints[iidx+1]), 2));
				if (dbnew >= l):
					xx = poly.xpoints[iidx+1] + ((l - db)/(dbnew - db))*(poly.xpoints[iidx] - poly.xpoints[iidx+1]);
					yy = poly.ypoints[iidx+1] + ((l - db)/(dbnew - db))*(poly.ypoints[iidx] - poly.ypoints[iidx+1]);
					dbnew = db + math.sqrt(math.pow(((l - db)/(dbnew - db))*(poly.xpoints[iidx] - poly.xpoints[iidx+1]), 2)  + 
									math.pow(((l - db)/(dbnew - db))*(poly.ypoints[iidx] - poly.ypoints[iidx+1]), 2));
				else:
					iidx = iidx-1;
				db = dbnew;
			if (db == l):
				pp1 = (xx, yy);
				pcp = (x, y);
				# then look forwards ONLY IF backwards search was successful...
				iidx = idx + 1;
				df = 0;
				while ((iidx < poly.npoints) and (df < l)):
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
					cp.append(pcp);
					p2.append((xx, yy));
	return (p1, cp, p2);

def calculate_curvature_profile(curv_points, roi, remove_negative_curvatures):
	"""generate a line profile of local curvatures using three-point method and SSS theorem (see http://mathworld.wolfram.com/SSSTheorem.html)"""
	poly = roi.getInterpolatedPolygon(1.0, True);
	curvature_profile = [((x,y),0) for (x,y) in zip(poly.xpoints, poly.ypoints)];
	pos = [p for (p, c) in curvature_profile]
	for (cp, p1, p2) in zip(curv_points[1], curv_points[0], curv_points[2]):
		a = math.sqrt( math.pow((cp[0] - p1[0]), 2) + math.pow((cp[1] - p1[1]), 2) );
		b = math.sqrt( math.pow((cp[0] - p2[0]), 2) + math.pow((cp[1] - p2[1]), 2) );
		c = math.sqrt( math.pow((p2[0] - p1[0]), 2) + math.pow((p2[1] - p1[1]), 2) );
		s = 0.5 * (a + b + c);
		try:
			K = math.sqrt(s * (s - a) * (s - b) * (s - c));
		except ValueError:
			if ((s - a) < 0) | ((s - b) < 0) | ((s - c) < 0):
				K = 0;
			if s < 0:
				raise ValueError('s < 0!');
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

def keep_largest_blob(imp):
	"""remove all blobs other than the largest by area"""
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
		indices_to_remove = [a for a in range(0,len(rt_areas)) if a != mx_ind];
		indices_to_remove.reverse();
		for rem_idx in indices_to_remove:
			roim.select(imp, rem_idx);
			IJ.run(imp, "Set...", "value=0 slice");
	imp.killRoi();
	roim.reset();
	roim.close();

def maximum_line_profile(imp, roi, pixel_width):
	"""return a line profile taking the maximum value over n pixels perpendicular to roi line"""
	imp.setRoi(roi);
	IJ.run(imp, "Interpolate", "interval=1.0 smooth adjust");
	ip = Straightener().straightenLine(imp, pixel_width);
	width = ip.getWidth();
	height = ip.getHeight();
	max_profile = [];
	poly = roi.getInterpolatedPolygon(1.0, True);
	for idx,(x,y) in enumerate(zip(poly.xpoints,poly.ypoints)):
		pix = ip.getLine(idx, 0, idx, height);
		max_profile.append(((x, y), max(pix)));
	return max_profile;

def roi_length(membrane_edge):
	"""calculate the length of the drawn membrane"""
	return membrane_edge.getLength();

def bleb_area(membrane_edge):
	"""calculate the area of the drawn bleb, accounting for membrane crossing line joining anchors"""
	poly = membrane_edge.getFloatPolygon();
	xs = [x for x in poly.xpoints];
	ys = [y for y in poly.ypoints];
	rotangle = membrane_edge.getAngle(int(round(xs[0])), int(round(ys[0])), int(round(xs[-1])), int(round(ys[-1]))) / 180 * math.pi;
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
	#print(area);
	return length, area, area_roi;

def apply_photobleach_correction_stack(params, actin_channel_imp):
	"""if toggled on, scale stack of labeled-species images to have constant mean intensity"""
	if not params.photobleaching_correction:
		return actin_channel_imp;
	# get intial mean labeled-species intensity:
	actin_channel_imp.killRoi();
	actin_channel_imp.setPosition(1);
	mean_t0 = actin_channel_imp.getStatistics().mean;
	# scale subsequent frames
	for idx in range(2, actin_channel_imp.getNFrames()+1):
		actin_channel_imp.setPosition(idx);
		mean_val = actin_channel_imp.getStatistics().mean;
		factor = mean_t0/mean_val;
		IJ.run(actin_channel_imp, "Multiply...", "value=" + str(factor) + " slice");
	return actin_channel_imp;

def apply_photobleach_correction_framewise(params, actin_channel_imp, membrane_channel_imp, t0_value=None):
	"""if toggled on, scale the current frame in a stack to have the same mean intensity as the first frame within an ROI"""
	if not params.photobleaching_correction:
		return actin_channel_imp, None;
	else:
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

def calculate_percentile(imp, roi, percentile):
	""" Use a crude (slow) method to get the pixel values in an roi, then calculate a percentile"""
	ip = imp.getProcessor();
	if (percentile > 100):
		raise ValueError("Percentile should be either in % or as a fraction, not > 100!");
	if (percentile > 1):
		percentile = float(percentile) / 100;
	pts = roi.getContainedPoints();
	vals = [];
	for pt in pts:
		vals.append(ip.getf(pt.x, pt.y));
	vals.sort();
	split_idx = float(len(vals)) * percentile + 0.5 - 1; # -1 as indexed from zero
	if (split_idx == math.floor(split_idx)):
		npc_percentile = vals[int(split_idx)];
	else:
		k = int(math.floor(split_idx));
		#f = split_idx - math.floor(split_idx);
		## calculate percentile according to 
		## https://web.stanford.edu/class/archive/anthsci/anthsci192/anthsci192.1064/handouts/calculating%20percentiles.pdf
		# npc_percentile = (1 - f) * vals[k] + f * vals[k+1];

		# calculate excel-style percentile...
		npc_percentile = vals[k] + (vals[k+1] - vals[k]) * (1 - percentile);
	return npc_percentile;

def flip_edge(roi, anchors):
	"""Check whether the edge is drawn in the expected orientation, and flip if not"""
	manual_roi_start = (roi.getPolygon().xpoints[0], roi.getPolygon().ypoints[0]);
	if vector_length(manual_roi_start, anchors[0]) < vector_length(manual_roi_start, anchors[1]):
		new_anchors = list(reversed(anchors));
	else:
		new_anchors = anchors;
	angle = angle_between_vecs(new_anchors[0], new_anchors[1], 
								(roi.getPolygon().xpoints[0], roi.getPolygon().ypoints[0]), 
								(roi.getPolygon().xpoints[-1], roi.getPolygon().ypoints[-1]))
	if angle < 0:
		xs = [x for x in roi.getPolygon().xpoints];
		ys = [y for y in roi.getPolygon().ypoints];
		xs.reverse();
		ys.reverse();
		roi = PolygonRoi(xs, ys, Roi.POLYLINE);
	return roi;