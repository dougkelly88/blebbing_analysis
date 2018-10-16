import math, os, sys, unittest
script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(script_path, 'modules'));

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
		ip = FloatProcessor(3, 5);
		imp = ImagePlus("testimp", ip);
		ypts = [y for y in range(1, imp.getHeight() + 1)];
		xpts = [2 for y in range(1, imp.getHeight() + 1)];
		imp.show();
		roi = PolygonRoi(xpts, ypts, Roi.POLYLINE);
		self.assertEqual([(2,4), (2,2)], mb.fix_anchors_to_membrane([(1,2), (3,4)],roi))
		self.assertRaises(ValueError, mb.fix_anchors_to_membrane, [(1,2), (1,2)], roi)
		imp.close();

	def test_generate_l_spaced_points(self):
		ypts = [y for y in range(0, 101)];
		xpts = [2 for y in range(0, 101)];
		roi = PolygonRoi(xpts, ypts, Roi.POLYLINE);
		l_spaced_points = mb.generate_l_spaced_points(roi, 5);
		test_group1 = [l_spaced_points[0][0], l_spaced_points[1][0], l_spaced_points[2][0]];
		test_group2 = [l_spaced_points[0][-1], l_spaced_points[1][-1], l_spaced_points[2][-1]];
		self.assertEqual([(2,0),(2,5),(2,10)], test_group1);
		self.assertEqual([(2,90),(2,95),(2,100)], test_group2);


if __name__ in ['__builtin__','__main__']:
	suite = unittest.TestLoader().loadTestsFromTestCase(TestMbEngine)
	unittest.TextTestRunner(verbosity=2).run(suite)

