# class to contain inner/outer comparison data
#
# D. J. Kelly, 2018-10-26, douglas.kelly@riken.jp

class InnerOuterComparisonData(object):
	"""simple container for data output from inner/outer comparison"""
	def __init__(self, inner_means=[], outer_means=[], inner_sds=[], outer_sds=[]):
		self.inner_means = inner_means;
		self.outer_means = outer_means;
		self.inner_sds = inner_sds;
		self.outer_sds = outer_sds;