import math, os, sys, unittest
script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(script_path, 'modules'));
sys.path.insert(0, os.path.join(script_path, 'classes'));

import membrane_blebbing_fileio as mbio;
import membrane_blebbing_ui as mbui;
import membrane_blebbing_engine as mb;
import membrane_blebbing_figures as mbfig;

from ij import ImagePlus
from ij.gui import PolygonRoi, Roi
from ij.process import FloatProcessor

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
		roi = PolygonRoi(xpts, ypts, Roi.POLYLINE);
		self.assertEqual([(2,4), (2,2)], mb.fix_anchors_to_membrane([(1,2), (3,4)],roi))
		self.assertRaises(ValueError, mb.fix_anchors_to_membrane, [(1,2), (1,2)], roi)
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
		
	def test_roi_length(self):
		ypts = [y for y in range(0, 3)];
		xpts = [1 for y in range(0, 3)];
		roi = PolygonRoi(xpts, ypts, Roi.POLYLINE);
		self.assertEqual(mb.roi_length(roi), 2);

	def test_bleb_area(self):
		ip = FloatProcessor(10, 10);
		imp = ImagePlus("testimp", ip);
		
		xpts = [1, 2, 2, 3, 3, 4];
		ypts = [1, 1, 2, 2, 1, 1];
		roi = PolygonRoi(xpts, ypts, Roi.POLYLINE);
		imp.setRoi(roi);
		self.assertEqual(mb.bleb_area(roi, imp), 1);
		imp.killRoi()

		# test behaviour when membrane crosses the line joining the anchor points
		xpts = [1, 1, 2, 2, 3, 4, 5, 6, 6, 7, 7];
		ypts = [2, 1, 1, 3, 3, 3, 3, 3, 1, 1, 2];
		roi = PolygonRoi(xpts, ypts, Roi.POLYLINE);
		imp.setRoi(roi);
		area = mb.bleb_area(roi, imp)
		self.assertEqual(area, 4);
		imp.close();

if __name__ in ['__builtin__','__main__']:
	suite = unittest.TestLoader().loadTestsFromTestCase(TestMbEngine)
	unittest.TextTestRunner(verbosity=2).run(suite)

