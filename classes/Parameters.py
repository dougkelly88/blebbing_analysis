# class to contain analysis parameters
#
# D. J. Kelly, 2018-10-26, douglas.kelly@riken.jp

import json, os, platform
from ij import IJ
from ij.process import AutoThresholder

class Parameters:
	"""Class to hold analysis parameters"""

	_persist_parameters_filename = "IJ_membrane_blebbing_params.json";
	_persist_parameters_folder = "IJ_membrane_blebbing";
	_version_string = "1.0.12";

	def __init__(self, load_last_params = False,
						input_image_path = None, 
						output_path = None, 
						curvature_length_um = 1.0, 
						threshold_method = 'Moments', 
						spatial_crop = None, 
						manual_anchor_positions = None, 
						manual_anchor_midpoint = None, 
						curvature_overlay_lut_string = 'physics', 
						curvature_kymograph_lut_string = 'physics', 
						actin_kymograph_lut_string = 'Grays', 
						labeled_species = 'Actin', 
						filter_negative_curvatures = True, 
						perform_spatial_crop = False, 
						perform_time_crop = False, 
						time_crop_start_end = None, 
						close_on_completion = False, 
						metadata_source = "Image metadata", 
						metadata_source_file = None, 
						photobleaching_correction = False, 
						perform_user_qc = False, 
						intensity_profile_width_um = 0.325, 
						membrane_channel_number = 2, 
						use_single_channel = False, 
						inner_outer_comparison = False, 
						selected_series_index = 0):
		self.__blebbingparams__ = True;

		success = True;
		if load_last_params:
			success = self.loadLastParams();
		if (not load_last_params) or not success:
			if input_image_path is None:
				self.input_image_path = os.path.expanduser("~");
			else:
				self.input_image_path = input_image_path;

			self.output_path = output_path;
			self.close_on_completion = close_on_completion;

			self.curvature_length_um = curvature_length_um;
			self.threshold_method = threshold_method;
			self.spatial_crop = spatial_crop;
			self.time_crop_start_end = time_crop_start_end;
			self.manual_anchor_positions = manual_anchor_positions;
			self.manual_anchor_midpoint = manual_anchor_midpoint;

			self.curvature_overlay_lut_string = curvature_overlay_lut_string;
			self.curvature_kymograph_lut_string = curvature_kymograph_lut_string;
			self.actin_kymograph_lut_string = actin_kymograph_lut_string;

			self.labeled_species = labeled_species;
			
			self.perform_user_qc = perform_user_qc;
			self.photobleaching_correction = photobleaching_correction;
			self.filter_negative_curvatures = filter_negative_curvatures;
			self.perform_spatial_crop = perform_spatial_crop;
			self.perform_time_crop = perform_time_crop;
			
			self.metadata_source = metadata_source;
			self.metadata_source_file = metadata_source_file;
			self.pixel_physical_size = 1.0;
			self.pixel_unit = "um";
			self.frame_interval = 1.0;
			self.interval_unit = "s";
			self.intensity_profile_width_um = intensity_profile_width_um;

			self.membrane_channel_number = membrane_channel_number;
			self.use_single_channel = use_single_channel;

			self.inner_outer_comparison = inner_outer_comparison;
			self.selected_series_index = selected_series_index;

		self.software_version = Parameters._version_string;

	def setSelectedSeriesIndex(self, selected_series_index):
		self.selected_series_index =selected_series_index;

	def setIntensityProfileWidthUm(self, width_um):
		self.intensity_profile_width_um = width_um;

	def setDoInnerOuterComparison(self, inner_outer_comparison):
		self.inner_outer_comparison = inner_outer_comparison;

	def setUseSingleChannel(self, use_single_channel):
		self.use_single_channel = use_single_channel;

	def setMembraneChannelNumber(self, membrane_channel_number):
		self.membrane_channel_number = membrane_channel_number;

	def togglePerformUserQC(self, do_qc):
		self.perform_user_qc = do_qc;

	def togglePhotobleachingCorrection(self, photobleaching_correction):
		self.photobleaching_correction = photobleaching_correction;

	def setMetadataSource(self, metadata_source):
		self.metadata_source = metadata_source;
		
	def setMetadataSourceFile(self, metadata_source_file):
		self.metadata_source_file = metadata_source_file;

	def toggleCloseOnCompletion(self, do_close_on_completion):
		self.close_on_completion = do_close_on_completion;

	def setTimeCropStartEnd(self, start_and_end_tuple):
		self.time_crop_start_end = start_and_end_tuple;

	def toggleTimeCrop(self, do_time_crop):
		self.perform_time_crop = do_time_crop;

	def toggleSpatialCrop(self, do_spatial_crop):
		self.perform_spatial_crop = do_spatial_crop;

	def setPixelSizeUnit(self, pixel_unit):
		self.pixel_unit = pixel_unit;

	def setIntervalUnit(self, interval_unit):
		self.interval_unit = interval_unit;

	def setPixelPhysicalSize(self, pixel_physical_size):
		self.pixel_physical_size = pixel_physical_size;

	def setFrameInterval(self, frame_interval):
		self.frame_interval = frame_interval;

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
		if method in self.listThresholdMethods():
			self.threshold_method = method;
		else: 
			raise ValueError('Requested threshold methd is not a valid IJ threshold method')

	def listThresholdMethods(self):
		threshold_methods = AutoThresholder.getMethods();
		local_threshold_methods = ["Bernsen", "Contrast", "Mean", "Median", "MidGrey", "Niblack", "Otsu", "Phansalkar", "Sauvola"]; # from https://github.com/fiji/Auto_Local_Threshold/blob/master/src/main/java/fiji/threshold/Auto_Local_Threshold.java (2018-11-02)
		for meth in local_threshold_methods:
			threshold_methods.append("Local: " + meth);
		threshold_methods.append("Edge + mean");
		return threshold_methods;

	def setCurvatureLengthUm(self, length):
		if length > 0:
			self.curvature_length_um = length;
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
				self.setCurvatureLengthUm(dct["curvature_length_um"]);
				self.setThresholdMethod(dct["threshold_method"])
				self.setSpatialCrop(dct["spatial_crop"]);
				self.setManualAnchorPositions(dct["manual_anchor_positions"]);
				self.setManualAnchorMidpoint(dct["manual_anchor_midpoint"]);
				self.setCurvatureOverlayLUT(dct["curvature_overlay_lut_string"]);
				self.setCurvatureKymographLUT(dct["curvature_kymograph_lut_string"]);
				self.setActinKymographLUT(dct["actin_kymograph_lut_string"]);
				self.setLabeledSpecies(dct["labeled_species"]);
				self.setFilterNegativeCurvatures(dct["filter_negative_curvatures"]);
				self.setTimeCropStartEnd(dct["time_crop_start_end"]);
				self.toggleSpatialCrop(dct["perform_spatial_crop"]);
				self.toggleTimeCrop(dct["perform_time_crop"]);
				self.toggleCloseOnCompletion(dct["close_on_completion"]);
				self.setMetadataSource(dct["metadata_source"]);
				self.setMetadataSourceFile(dct["metadata_source_file"]);
				self.togglePhotobleachingCorrection(dct["photobleaching_correction"]);
				self.togglePerformUserQC(dct["perform_user_qc"]);
				self.setIntensityProfileWidthUm(dct["intensity_profile_width_um"]);
				self.setMembraneChannelNumber(dct["membrane_channel_number"]);
				self.setUseSingleChannel(dct["use_single_channel"]);
				self.setDoInnerOuterComparison(dct["inner_outer_comparison"]);
				self.setSelectedSeriesIndex(dct["selected_series_index"]);
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
			temp_path = self.getPeristenceFileLocation();
			if temp_path:
				temp_params_path = os.path.join(temp_path, Parameters._persist_parameters_filename);
				if os.path.isfile(temp_params_path):
					success = self.loadParametersFromJson(temp_params_path);
				else:
					success = False;
			else:
				success = False;
		except Exception as e:
			print("Warning: Error loading previous settings, reverting to default...");
			raise e;
			return False;
		if not success:
			print("Warning: Error loading previous settings, reverting to default...");
		return success;

	def persistParameters(self):
		temp_path = self.getPeristenceFileLocation();
		if temp_path:
			temp_params_path = os.path.join(temp_path, Parameters._persist_parameters_filename);
			self.saveParametersToJson(temp_params_path);

	def getPeristenceFileLocation(self):
		try:
			st = platform.mac_ver()[0];
			if not st:
				# windows
				temp_path = os.path.join(os.getenv('APPDATA'), Parameters._persist_parameters_folder);
			else:
				# mac
				temp_path = os.path.join(os.path.expanduser("~"), "Library", Parameters._persist_parameters_folder);
			if not os.path.isdir(temp_path):
				os.mkdir(temp_path);
		except Exception as e:
			print("Error: " + e.message);
			return "";
		return temp_path;

#params = Parameters();
#print(params.__dict__);
#params.saveParametersToJson("C:\\Users\\dougk\\Desktop\\test.json");
#params.loadParametersFromJson("C:\\Users\\dougk\\Desktop\\test2.json");
#print(params.__dict__);