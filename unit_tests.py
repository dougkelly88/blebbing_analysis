# unit tests for blebbing analysis code
#
# D. J. Kelly, 2018-10-16, douglas.kelly@riken.jp

import math, os, sys, unittest

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
from Parameters import Parameters;

from ij import ImagePlus, IJ
from ij.gui import PolygonRoi, Roi
from ij.process import FloatProcessor, ImageProcessor

class TestMbEngine(unittest.TestCase):
	def test_vector_length(self):
		self.assertEqual(1, mb.vector_length((0,0), (1,0)));
		self.assertEqual(1, mb.vector_length((0,0), (-1, 0)));
		self.assertEqual(0, mb.vector_length((0,0), (0, 0)));
	
	def test_angle_between_vecs(self):
		self.assertEqual(0, mb.angle_between_vecs((0,0),(1,0),(0,0),(1,0)));
		self.assertEqual(math.pi/2, mb.angle_between_vecs((0,0),(1,0),(0,0),(0,1)));
		self.assertEqual(math.pi, mb.angle_between_vecs((0,0),(1,0),(0,0),(-1,0)));
		self.assertEqual(-math.pi/2, mb.angle_between_vecs((0,0),(1,0),(0,0),(0,-1)));
		self.assertEqual(math.pi/2, mb.angle_between_vecs((0,0),(1,0),(0,-1),(0,0)));

	def test_fix_anchors_to_membrane(self):
		"""test that anchors are fixed to expected positions, and that degenerate anchors raise an error"""
		ip = FloatProcessor(3, 5);
		imp = ImagePlus("testimp", ip);
		ypts = [y for y in range(1, 6)];
		xpts = [2 for y in range(1, 6)];
		imp.show();
		params = Parameters(inner_outer_comparison=False);
		roi = PolygonRoi(xpts, ypts, Roi.POLYLINE);
		self.assertEqual([(2,4), (2,2)], mb.fix_anchors_to_membrane([(1,2), (3,4)],roi,params))
		self.assertRaises(ValueError, mb.fix_anchors_to_membrane, [(1,2), (1,2)], roi,params)
		imp.close();

	def test_generate_l_spaced_points(self):
		"""test that start and end of l-spaced points range are as expected"""
		ypts = [y for y in range(0, 101)];
		xpts = [2 for y in range(0, 101)];
		roi = PolygonRoi(xpts, ypts, Roi.POLYLINE);
		l_spaced_points = mb.generate_l_spaced_points(roi, 5);
		test_group1 = [l_spaced_points[0][0], l_spaced_points[1][0], l_spaced_points[2][0]];
		test_group2 = [l_spaced_points[0][-1], l_spaced_points[1][-1], l_spaced_points[2][-1]];
		self.assertEqual([(2,0),(2,5),(2,10)], test_group1);
		self.assertEqual([(2,90),(2,95),(2,100)], test_group2);
		
	def test_calculate_curvature_profile(self):
		"""test that flat line returns 0 curvature, r = 1000 curve returns 1/1000 curvature, and that negative curvatures are removed appropriately"""
		ypts = [y for y in range(0, 3)];
		xpts = [1 for y in range(0, 3)];
		roi = PolygonRoi(xpts, ypts, Roi.POLYLINE);
		cp = mb.calculate_curvature_profile(([(1,0)], [(1,1)], [(1,2)]), roi, True);
		ctest1 = [c for ((x,y),c) in cp];
		self.assertEqual([0,0,0], ctest1);

		theta = [t*math.pi/4 for t in range(0,3)]
		ypts = [1000*math.sin(t) for t in theta];
		xpts = [-1000*math.cos(t) for t in theta];
		roi = PolygonRoi(xpts, ypts, Roi.POLYLINE);
		# N.B this length tends to l = pi * r/4 ~=785.4 as precision/length of theta
		# vector increases. Gets silly, so just do by inspection...
		l = 764.5;
		pts = mb.generate_l_spaced_points(roi, l) 
		cp = mb.calculate_curvature_profile(pts, roi, True)
		ctest2 = cp[int(round(len(cp)/2))][1];
		self.assertAlmostEqual(ctest2, 1.0/1000, 5);

		xpts = [1000*math.cos(t) for t in theta];
		roi = PolygonRoi(xpts, ypts, Roi.POLYLINE);
		pts = mb.generate_l_spaced_points(roi, l) 
		cp = mb.calculate_curvature_profile(pts, roi, False)
		ctest3 = cp[int(round(len(cp)/2))][1];
		self.assertAlmostEqual(ctest3, -1.0/1000, 5);

		cp = mb.calculate_curvature_profile(pts, roi, True)
		ctest4 = cp[int(round(len(cp)/2))][1];
		self.assertEqual(ctest4, 0);
		
	def test2_calculate_curvature_profile(self):
		"""test that signs attributed to curvatures are correct for different edge orientations"""
		w = 500; h = 500;
		R = 200;
		length_param_pix = 50;
		x1 = []; y1 = []; c1 = (0, int(float(h)/2));
		x2 = []; y2 = []; c2 = (int(float(w)/2),  0);
		x3 = []; y3 = []; c3 = (w, int(float(h)));
		x4 = []; y4 = []; c4 = (int(float(w)/2), h);
		thetas = [((float(tidx)/1000) * (2 * math.pi)) - math.pi/4 for tidx in range(1000)]
		for theta in thetas:
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
		edges = [PolygonRoi(x1, y1, Roi.POLYLINE), 
			PolygonRoi(x2, y2, Roi.POLYLINE), 
			PolygonRoi(x3, y3, Roi.POLYLINE), 
			PolygonRoi(x4, y4, Roi.POLYLINE)];
		anchorses = [[(x1[0], y1[0]), (x1[-1], y1[-1])], 
				[(x2[0], y2[0]), (x2[-1], y2[-1])], 
				[(x3[0], y3[0]), (x3[-1], y3[-1])], 
				[(x4[0], y4[0]), (x4[-1], y4[-1])]];
		cs = [c1, c2, c3, c4];
		for eidx, (edge, c, anchors) in enumerate(zip(edges, cs, anchorses)):
			for mp_dir in range(2):
				if not mp_dir:
					midpoint = c;
				else:
					midpoint = (int(float(w)/2), int(float(h)/2));
				anchors = mb.order_anchors(anchors, [midpoint]);
				edge = mb.check_edge_order(anchors, edge);
				curv_pts = mb.generate_l_spaced_points(edge, length_param_pix);
				curv_profile = mb.calculate_curvature_profile(curv_pts,
														edge, 
														False);
				curvs_only = [cv for cp, cv in curv_profile];
				mean_curv = sum(curvs_only)/len(curvs_only);
				self.assertEqual(mean_curv > 0, mp_dir);
				

	def test_roi_length(self):
		ypts = [y for y in range(0, 3)];
		xpts = [1 for y in range(0, 3)];
		roi = PolygonRoi(xpts, ypts, Roi.POLYLINE);
		self.assertEqual(mb.roi_length(roi), 2);

	def test_bleb_area(self):	
		"""test behaviour when membrane crosses the line joining the anchor points"""
		xpts = [1, 1, 2, 2, 3, 4, 5, 6, 6, 7, 7];
		ypts = [2, 1, 1, 3, 3, 3, 3, 3, 1, 1, 2];
		roi = PolygonRoi(xpts, ypts, Roi.POLYLINE);
		length, area, area_roi = mb.bleb_area(roi)
		self.assertEqual(area, 4);

	def test_keep_largest_blob(self):
		"""test to ensure that keep largest blob does indeed only retain largest blob"""
		# create binary image
		w = 5;
		imp = IJ.createImage("one", w, w, 1, 8);
		pix = imp.getProcessor().getPixels();

		for idx in range(0, w):
			pix[idx] = 127;
		pix[w * w - w + 2] = 127;
		pix[w * w - w] = 127;
		pix[int(float(w)/2 * float(w+1)/2 + w)] = 127;
		pix[int(float(w)/2 * float(w+1)/2 + w) + 1] = 127
		IJ.setRawThreshold(imp, 125, 255, "null");
		IJ.run(imp, "Convert to Mask", "");
		mb.keep_largest_blob(imp);
		IJ.run(imp, "Create Selection", "");
		roi = imp.getRoi();
		self.assertEqual(roi.getStatistics().area, w);

	def test_maximum_line_profile(self):
		"""test that maximum line profile is generated as expected for straight line"""
		imp = IJ.createImage("test", 3, 4, 1, 8);
		ip = imp.getProcessor();
		pix = ip.getPixels();
		pix[0] = 10; pix[1] = 20; pix[2] = 30;
		pix[3] = 60; pix[4] = 40; pix[5] = 50;
		pix[6] = 80; pix[7] = 90; pix[8] = 70;
		pix[9] = 100; pix[10] = 110; pix[11] = 120;
		roi = PolygonRoi([1, 1], [0, 3], Roi.POLYLINE);
		imp.setRoi(roi);
		mlp = mb.maximum_line_profile(imp, roi, 3);
		vals = [i for ((x, y), i) in mlp];
		self.assertEqual(vals, [30, 60, 90, 120]);
		
	def test_apply_photobleach_correction_framewise(self):
		"""confirm that photobleaching correction equalises mean values across a stack"""
		# create test images: 
		nx = 3;
		ny = 3;
		nz = 1;
		nc = 1;
		nt = 5;
		bitdepth = 8;
		seg_imp = IJ.createHyperStack("segmentation image", nx, ny, nc, nz, nt, bitdepth);
		intensity_imp = IJ.createHyperStack("intensity image", nx, ny, nc, nz, nt, bitdepth);
		for idx in range(0, nt):
			seg_imp.setPosition(idx + 1);
			seg_pix = seg_imp.getProcessor().getPixels();
			for pixidx in range(len(seg_pix)):
				seg_pix[pixidx] = 1;
			intensity_imp.setPosition(idx + 1);
			int_pix = intensity_imp.getProcessor().getPixels();
			for pixidx in range(len(int_pix)):
				int_pix[pixidx] = nt - idx;
		ip = seg_imp.getProcessor();
		ip.setThreshold(8, 8, ImageProcessor.NO_LUT_UPDATE);
		IJ.run(seg_imp, "Make Binary", "method=Default background=Dark calculate");

		# apply photobleaching correction: 
		t0_mean = None;
		params = Parameters(photobleaching_correction = True);
		for idx in range(0, nt):
			intensity_imp.setPosition(idx+1);
			seg_imp.setPosition(idx+1);
			intensity_imp, t0_mean = mb.apply_photobleach_correction_framewise(params, 
																			intensity_imp, 
																			seg_imp, 
																			t0_mean);
		# test result:
		for idx in range(0, nt):
			intensity_imp.setPosition(idx+1);
			self.assertEqual(nt, intensity_imp.getStatistics().mean);


if __name__ in ['__builtin__','__main__']:
	suite = unittest.TestLoader().loadTestsFromTestCase(TestMbEngine)
	unittest.TextTestRunner(verbosity=2).run(suite)

