from ij import IJ, ImageStack
from ij.gui import PolygonRoi, Roi, WaitForUserDialog
import math, os, sys

release = False;

if not release:
	script_path = os.path.dirname(os.path.realpath(__file__));
else: 
	script_path = os.getcwd();
if "Fiji.app" in script_path:
	ss = script_path.split("Fiji.app");
	final_folder = "blebbing analysis";
	script_path = os.path.join(ss[0], "Fiji.app", "plugins", "Scripts", "Plugins", final_folder);
sys.path.insert(0, os.path.join(script_path, 'modules'));
sys.path.insert(0, os.path.join(script_path, 'classes'));

import membraneBlebbingFileio as mbio;
import membraneBlebbingUi as mbui;
import membraneBlebbingEngine as mb;
import membraneBlebbingFigures as mbfig;
from Parameters import Parameters

params = Parameters(output_path="C:\\Users\\dougk\\Desktop");
length_param_pix = 50;
imp = IJ.createImage("test", 500, 500, 1, 8);
imp.show()
R = 200;
x1 = []; y1 = []; c1 = (0, int(float(imp.getHeight())/2));
x2 = []; y2 = []; c2 = (int(float(imp.getWidth())/2),  0);
x3 = []; y3 = []; c3 = (imp.getWidth(), int(float(imp.getHeight())/2));
x4 = []; y4 = []; c4 = (int(float(imp.getWidth())/2), imp.getHeight());
x5 = []; y5 = []; c5 = (0, imp.getHeight());
x6 = []; y6 = []; c6 = (0, 0);
x7 = []; y7 = []; c7 = (imp.getWidth(), 0);
x8 = []; y8 = []; c8 = (imp.getWidth(), imp.getHeight());
c8 = (int(3 * float(imp.getWidth()) / 4), int(3 * float(imp.getHeight()) / 4));

for tidx in range(1000):
	theta = ((float(tidx)/1000) * (2 * math.pi)) - math.pi/4;
	if theta < math.pi/4:
		x1.append(R * math.cos(theta) + c1[0]);
		y1.append(R * math.sin(theta) + c1[1]);
	elif theta < 3 * math.pi/4:
		x2.append(R * math.cos(theta) + c2[0]);
		y2.append(R * math.sin(theta) + c2[1]);
	elif theta < 5 * math.pi/4:
		x3.append(R * math.cos(theta) + c3[0]);
		y3.append(R * math.sin(theta) + c3[1]);
	elif theta < 7 * math.pi/4:
		x4.append(R * math.cos(theta) + c4[0]);
		y4.append(R * math.sin(theta) + c4[1]);

	theta2 = ((float(tidx)/1000) * (2 * math.pi)) - math.pi/2
	if theta2 < 0:
		x5.append(R * math.cos(theta2) + c5[0]);
		y5.append(R * math.sin(theta2) + c5[1]);
	elif theta2 < math.pi/2:
		x6.append(R * math.cos(theta2) + c6[0]);
		y6.append(R * math.sin(theta2) + c6[1]);
	elif theta2 < math.pi:
		x7.append(R * math.cos(theta2) + c7[0]);
		y7.append(R * math.sin(theta2) + c7[1]);
	elif theta2 <  3 * math.pi/2:
		x8.append(R * math.cos(theta2) + c8[0]);
		y8.append(R * math.sin(theta2) + c8[1]);

edges = [PolygonRoi(x1, y1, Roi.POLYLINE), 
			PolygonRoi(x2, y2, Roi.POLYLINE), 
			PolygonRoi(x3, y3, Roi.POLYLINE), 
			PolygonRoi(x4, y4, Roi.POLYLINE), 
			PolygonRoi(x5, y5, Roi.POLYLINE), 
			PolygonRoi(x6, y6, Roi.POLYLINE), 
			PolygonRoi(x7, y7, Roi.POLYLINE), 
			PolygonRoi(x8, y8, Roi.POLYLINE)];
cs = [c1, c2, c3, c4, c5, c6, c7, c8];
anchorses = [[(x1[0], y1[0]), (x1[-1], y1[-1])], 
				[(x2[0], y2[0]), (x2[-1], y2[-1])], 
				[(x3[0], y3[0]), (x3[-1], y3[-1])], 
				[(x4[0], y4[0]), (x4[-1], y4[-1])], 
				[(x5[0], y5[0]), (x5[-1], y5[-1])], 
				[(x6[0], y6[0]), (x6[-1], y6[-1])], 
				[(x7[0], y7[0]), (x7[-1], y7[-1])], 
				[(x8[0], y8[0]), (x8[-1], y8[-1])]]

for eidx, (edge, c, anchors) in enumerate(zip(edges, cs, anchorses)):
	for mp_dir in range(2):
		if not mp_dir:
			midpoint = c;
			target_sign_str = " NEGATIVE";
		else:
			midpoint = (int(float(imp.getWidth())/2), int(float(imp.getHeight())/2));
			if eidx == 7:
				midpoint = (1, 1);
			target_sign_str = " POSITIVE";
		anchors = mb.order_anchors(anchors, [midpoint]);
		print("Ordered anchors = " + str(anchors));
		for idx, pt in enumerate([anchors[0], anchors[1], midpoint]):
			imp.setRoi(Roi(pt[0] - 3, pt[1] - 3, 7, 7));
			IJ.run(imp, "Set...", "value=" + str(idx+1) + " slice");
		imp.setDisplayRange(0,5);
		edge = mb.check_edge_order(anchors, edge);
		curv_pts = mb.generate_l_spaced_points(edge, length_param_pix);
		curv_profile = mb.calculate_curvature_profile(curv_pts,
														edge, 
														False);
		curvs_only = [cv for cp, cv in curv_profile];
		mean_curv = sum(curvs_only)/len(curvs_only);
		print("Average curve profile value = " + str(mean_curv));
		if mean_curv > 0:
			sign_str = " POSITIVE";
		else:
			sign_str = " NEGATIVE";
		imp.setRoi(PolygonRoi(list(edge.getPolygon().xpoints), list(edge.getPolygon().ypoints), Roi.FREELINE));
		curv_stack = ImageStack(imp.getWidth(), imp.getHeight());
		curv_stack = mbfig.generate_curvature_overlays(curv_profile, curv_stack);
		curv_stack = mbfig.generate_curvature_overlays(curv_profile, curv_stack);
		overlaid_curvature_imp, raw_curvature_imp = mbfig.overlay_curvatures(imp, curv_stack, [curv_profile, curv_profile], 1, params, annotate=True)
		overlaid_curvature_imp.show();
		raw_curvature_imp.show();
		WaitForUserDialog(str(eidx) + ", direction " + str(mp_dir) + ": " + sign_str + ", should be " + target_sign_str).show();
		overlaid_curvature_imp.close();
		raw_curvature_imp.close();
		for idx, pt in enumerate([anchors[0], anchors[1], midpoint]):
			imp.setRoi(Roi(pt[0] - 3, pt[1] - 3, 7, 7));
			IJ.run(imp, "Set...", "value=0 slice");

imp.changes = False;
imp.close();