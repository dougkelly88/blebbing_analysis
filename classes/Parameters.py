# class to contain analysis parameters
#
# D. J. Kelly, 2018-10-26, douglas.kelly@riken.jp

import json, os, platform
from ij import IJ
from ij.process import AutoThresholder

class Parameters:
	"""Class to hold analysis parameters"""

	_persist_parameters_filename = "IJ_membrane_blebbing_params.json";

	def __init__(self, load_last_params = False,
						input_image_path = None, 
						output_path = None, 
						curvature_length_pix = 10.0, 
						threshold_method = 'Moments', 
						spatial_crop = None, 
						manual_anchor_positions = None, 
						manual_anchor_midpoint = None, 
						curvature_overlay_lut_string = 'physics', 
						curvature_kymograph_lut_string = 'Yellow', 
						actin_kymograph_lut_string = 'Cyan', 
						labeled_species = 'Actin', 
						filter_negative_curvatures = True):
		self.__blebbingparams__ = True;

		success = True;
		if load_last_params:
			success = self.loadLastParams();
		if (not load_last_params) or not success:
			self.input_image_path = input_image_path;
			self.output_path = output_path;

			self.curvature_length_pix = curvature_length_pix;
			self.threshold_method = threshold_method;
			self.spatial_crop = spatial_crop;
			self.manual_anchor_positions = manual_anchor_positions;
			self.manual_anchor_midpoint = manual_anchor_midpoint;

			self.curvature_overlay_lut_string = curvature_overlay_lut_string;
			self.curvature_kymograph_lut_string = curvature_kymograph_lut_string;
			self.actin_kymograph_lut_string = actin_kymograph_lut_string;

			self.labeled_species = labeled_species;
			
			self.filter_negative_curvatures = filter_negative_curvatures;

	def setFilterNegativeCurvatures(self, dofilter):
		self.filter_negative_curvatures = dofilter;

	def setLabeledSpecies(self, labeled_species):
		if not labeled_species:
			labeled_species = "DefaultLabeledSpecies";
		self.labeled_species = labeled_species;

	def setOutputPath(self, path):
		# note that validating paths isn't trivial: https://stackoverflow.com/questions/9532499/check-whether-a-path-is-valid-in-python-without-creating-a-file-at-the-paths-ta
		self.output_path = path;

	def setInputImagePath(self, path):
		self.input_image_path = path;

	def setSpatialCrop(self, spatial_crop):
		self.spatial_crop = spatial_crop;

	def setManualAnchorPositions(self, manual_anchor_positions):
		self.manual_anchor_positions = manual_anchor_positions;

	def setManualAnchorMidpoint(self, manual_anchor_midpoint):
		self.manual_anchor_midpoint = manual_anchor_midpoint;

	def setThresholdMethod(self, method):
		if method in AutoThresholder.getMethods():
			self.threshold_method = method;
		else: 
			raise ValueError('Requested threshold methd is not a valid IJ threshold method')

	def setCurvatureLengthPix(self, length):
		if length > 0:
			self.curvature_length_pix = length;
		else:
			raise ValueError('Curvature length parameter must be positive');

	def setCurvatureOverlayLUT(self, lut_string):
		if lut_string in IJ.getLuts():
			self.curvature_overlay_lut_string = lut_string;
		else:
			raise ValueError('Requested curvature overlay LUT is not a valid IJ LUT');
			
	def setCurvatureKymographLUT(self, lut_string):
		if lut_string in IJ.getLuts():
			self.curvature_kymograph_lut_string = lut_string;
		else:
			raise ValueError('Requested curvature kymograph LUT is not a valid IJ LUT');

	def setActinKymographLUT(self, lut_string):
		if lut_string in IJ.getLuts():
			self.actin_kymograph_lut_string = lut_string;
		else:
			raise ValueError('Requested actin kymograph LUT is not a valid IJ LUT');

	def saveParametersToJson(self, file_path):
		try:
			f = open(file_path, 'w');
			json.dump(self.__dict__, f);
		finally:
			f.close();

	def loadParametersFromJson(self, file_path):
		try:
			f = open(file_path, 'r');
			dct = json.loads(f.read());
			if "__blebbingparams__" in dct:
				self.setInputImagePath(dct["input_image_path"]);
				self.setOutputPath(dct["output_path"]);
				self.setCurvatureLengthPix(dct["curvature_length_pix"]);
				self.setThresholdMethod(dct["threshold_method"])
				self.setSpatialCrop(dct["spatial_crop"]);
				self.setManualAnchorPositions(dct["manual_anchor_positions"]);
				self.setManualAnchorMidpoint(dct["manual_anchor_midpoint"]);
				self.setCurvatureOverlayLUT(dct["curvature_overlay_lut_string"]);
				self.setCurvatureKymographLUT(dct["curvature_kymograph_lut_string"]);
				self.setActinKymographLUT(dct["actin_kymograph_lut_string"]);
				self.setLabeledSpecies(dct["labeled_species"]);
				self.setFilterNegativeCurvatures(dct["filter_negative_curvatures"]);
			else:
				raise ValueError("JSON file doesn't translate to membrane blebbing analysis parameters")
		except IOError:
			print("IOError reading from JSON file");
			return False;
		except: 
			return False;
		finally:
			f.close();
		return True;

	def loadLastParams(self):
		success = True;
		try:
			st = platform.mac_ver()[0];
			if not st:
				# windows
				temp_path = os.path.join(os.getenv('APPDATA'), "IJ_membrane_blebbing");
			else:
				# mac
				temp_path = os.path.join('~/Library/', "IJ_membrane_blebbing");
			if not os.path.isdir(temp_path):
				os.mkdir(temp_path);
			temp_params_path = os.path.join(temp_path, Parameters._persist_parameters_filename);
			if os.path.isfile(temp_params_path):
				success = self.loadParametersFromJson(temp_params_path);
		except:
			raise Warning("Error loading previous settings, reverting to default...");
			return False;
		if not success:
			raise Warning("Error loading previous settings, reverting to default...");
		return success;

	def persistParameters(self):
		st = platform.mac_ver()[0];
		if not st:
			# windows
			temp_path = os.path.join(os.getenv('APPDATA'), "IJ_membrane_blebbing");
		else:
			# mac
			temp_path = os.path.join('~/Library/', "IJ_membrane_blebbing");
		if not os.path.isdir(temp_path):
			os.mkdir(temp_path);
		temp_params_path = os.path.join(temp_path, Parameters._persist_parameters_filename);
		self.saveParametersToJson(temp_params_path);

#params = Parameters();
#print(params.__dict__);
#params.saveParametersToJson("C:\\Users\\dougk\\Desktop\\test.json");
#params.loadParametersFromJson("C:\\Users\\dougk\\Desktop\\test2.json");
#print(params.__dict__);