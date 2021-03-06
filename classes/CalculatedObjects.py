# class to contain calculated outputs to improve tidiness passing 
# between functions and to facilitate saving intermediate data if necessary later
#
# D. J. Kelly, 2019-01-09, douglas.kelly@riken.jp

from InnerOuterComparisonData import InnerOuterComparisonData

class CalculatedObjects(object):
	"""simple class to contain all intermediate calculations"""
	def __init__(self, 
				 membrane_edges=None, 
				 inner_outer_data=InnerOuterComparisonData(), 
				 fixed_anchors_list=None,
				 curvature_profiles=None,
				 actin_profiles=None, 
				 bleb_perimeter_lengths=None,
				 bleb_areas=None, 
				 timelist=None,
				 background_sd_profile=None, 
				 full_membrane_lengths=None, 
				 euclidean_membrane_lengths=None):
		self.membrane_edges = membrane_edges;
		self.inner_outer_data = inner_outer_data;
		self.fixed_anchors_list = fixed_anchors_list;
		self.curvature_profiles = curvature_profiles;
		self.actin_profiles = actin_profiles;
		self.bleb_perimeter_lengths = bleb_perimeter_lengths;
		self.bleb_areas = bleb_areas;
		self.timelist = timelist
		self.background_sd_profile = background_sd_profile;
		self.full_membrane_lengths = full_membrane_lengths;
		self.euclidean_membrane_lengths = euclidean_membrane_lengths;