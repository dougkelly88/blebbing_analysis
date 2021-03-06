# functions for handling file input/output for membrane blebbing analysis
#
# D. J. Kelly, 2018-10-15, douglas.kelly@riken.jp

# imports 
import csv, json, os
from datetime import datetime
from ij.io import OpenDialog, DirectoryChooser
from ij.gui import PolygonRoi, Roi
from loci.formats import ImageReader, MetadataTools
from ij.plugin.frame import RoiManager
from ome.units import UNITS

import membraneBlebbingUi as mbui

def file_location_chooser(default_directory):
	"""choose folder locations and prepare output folder"""
	# input
	file_path = input_file_location_chooser(default_directory);
	# output
	DirectoryChooser.setDefaultDirectory(os.path.dirname(file_path));
	dc = DirectoryChooser('Select the root folder for saving output');
	output_root = dc.getDirectory();
	if output_root is None:
		raise IOError('no output path chosen');
	timestamp = datetime.strftime(datetime.now(), '%Y-%m-%d %H-%M-%S');
	output_folder = os.path.join(output_root, (timestamp + ' output'));
	os.mkdir(output_folder);
	return file_path, output_folder;

def input_file_location_chooser(default_directory, filt='*.tif', message=None):
	if message==None:
		message = 'Choose original file...'
	od = OpenDialog(message, 
					default_directory, 
					filt);
	file_path = od.getPath();
	if file_path is None:
		raise IOError('no input file chosen');
	return file_path;

def rerun_location_chooser(default_filepath):
	"""choose folder containing a previous analysis run to reanalyse"""
	DirectoryChooser.setDefaultDirectory(os.path.dirname(default_filepath));
	dc = DirectoryChooser('Select the folder containing the previous analysis output...');
	old_output_folder = dc.getDirectory();
	if old_output_folder is None:
		raise IOError('no input path chosen');
	# check that chosen folder contains the right files...
	files_lst = os.listdir(old_output_folder);
	#if not all([f in files_lst for f in ['user_defined_edges.zip', 'parameters used.json']]):
	if not 'parameters used.json' in files_lst or not(('user_defined_edges.json' in files_lst) or ('user_defined_edges.zip' in files_lst)):
		raise IOError('chosen path isn''t a valid membrane blebbing output folder');
	timestamp = datetime.strftime(datetime.now(), '%Y-%m-%d %H-%M-%S');
	new_output_folder = os.path.join(os.path.dirname(os.path.normpath(old_output_folder)), (timestamp + ' output'));
	os.mkdir(new_output_folder);
	return old_output_folder, new_output_folder;

def save_profile_as_csv(profiles, file_path, data_name, xname='x', yname='y', tname='Frame', time_list=[]):
	"""save profiles with 2d independent variable, e.g. curvature profile"""
	f = open(file_path, 'wb');
	writer = csv.writer(f);
	writer.writerow([tname, xname, yname, data_name]);
	if not time_list:
			time_list = [idx for idx in range(0, len(profiles))];
	if len(time_list)>1:
		for idx, profile in enumerate(profiles):
			for subidx, ((x, y), p) in enumerate(profile):
				if len(profiles) == 1:
					writer.writerow([time_list[subidx], x, y, p]);
				else:
					writer.writerow([time_list[idx], x, y, p]);
	else:
		for idx, ((x,y), p) in enumerate(profiles[0]):
			writer.writerow([time_list[0], x, y, p])
	f.close();

def load_csv_as_profile(file_path):
	"""load profiles with 2d independent variable from csv"""
	f = open(file_path, 'r');
	reader = csv.reader(f);
	old_frame = 0;
	profiles = [];
	profile = [];
	for idx, row in enumerate(reader):
		if idx > 0:
			frame = int(row[0]);
			if frame != old_frame:
				profiles.append(profile);
				profile = [];
				old_frame = frame;
			profile.append(((float(row[1]), float(row[2])), float(row[3])));
	profiles.append(profile);
	return profiles;
			 
def save_1d_profile_as_csv(profile, file_path, column_names):
	"""save profiles with 1d independent variable, e.g. length against time"""
	f = open(file_path, 'wb');
	writer = csv.writer(f);
	writer.writerow([column_names[0], column_names[1]]);
	for p in profile:
		writer.writerow([p[0], p[1]]);
	f.close();

def save_parameters(params, file_path):
	"""save parameters used for this analysis"""
	f = open(file_path, 'w');
	json.dump(params, f);
	f.close();

def save_qcd_edges(edges, output_folder):
	"""save JSON describing user-drawn edges"""
	edge_point_list = [zip(poly.xpoints, poly.ypoints) for poly in [edge.getInterpolatedPolygon(1.0, False) for edge in edges]];
	file_path = os.path.join(output_folder, "user_defined_edges.json");
	f = open(file_path, 'w');
	try:
		json.dump(edge_point_list, f);
	finally:
		f.close();

def save_qcd_edges2(edges, output_folder, filename="user_defined_edges.zip"):
	"""save edges as rois to a *.zip file"""
	roim = RoiManager(False)
	for edge in edges:
		if edge is not None:
			roim.addRoi(edge);
	roim.runCommand("Save", os.path.join(output_folder, filename));
	roim.close();

def load_qcd_edges(input_file_path):
	"""load edges from JSON"""
	f = open(input_file_path, 'r');
	try:
		edges = json.loads(f.read());
	finally:
		f.close();
	membrane_edges = [];
	for edge in edges:
		xs = [pt[0] for pt in edge];
		ys = [pt[1] for pt in edge];
		membrane_edges.append(PolygonRoi(xs, ys, Roi.POLYLINE));
	return membrane_edges;

def load_qcd_edges2(input_file_path):
	"""load edges from roi *.zip file"""
	if os.path.isfile(input_file_path):
		roim = RoiManager(False);
		roim.runCommand("Open", input_file_path);
		edges = roim.getRoisAsArray();
		roim.close();
	else:
		edges = load_qcd_edges(os.path.splitext(input_file_path)[0] + '.json');
	return edges;

def get_metadata(params):
	"""get image metadata, either from the image file or from acquisition-time metadata"""
	if params.metadata_source == "Image metadata":
		try:
			reader = ImageReader();
			ome_meta = MetadataTools.createOMEXMLMetadata();
			reader.setMetadataStore(ome_meta);
			reader.setId(params.input_image_path);
			reader.close();
			params.setFrameInterval(ome_meta.getPixelsTimeIncrement(0).value());
			params.setIntervalUnit(ome_meta.getPixelsTimeIncrement(0).unit().getSymbol());
			params.setPixelPhysicalSize(ome_meta.getPixelsPhysicalSizeX(0).value());
			params.setPixelSizeUnit(ome_meta.getPixelsPhysicalSizeX(0).unit().getSymbol());
			params.setMetadataSourceFile(None);
		except Exception as e:
			print(e.message);
			mbui.warning_dialog(["There was a problem getting metadata from the image: ", 
									e.message, 
								"Please consider using acquisition metadata instead (click OK). ", 
								"Or, quit the analysis run and investigate image metadata by hand. "]);
			params.setMetadataSource("Acquisition metadata")
	if params.metadata_source == "Acquisition metadata":
		od = OpenDialog('Choose acquisition metadata file...', 
					os.path.dirname(params.input_image_path), 
					'*.txt');
		file_path = od.getPath();
		if file_path is None:
			raise IOError('no metadata file chosen');
		acq_metadata_dict = import_iq3_metadata(file_path);
		try:
			params.setFrameInterval(acq_metadata_dict['frame_interval']);
		except KeyError:
			params.setFrameInterval(1.0);
		try:
			params.setIntervalUnit(acq_metadata_dict['time_unit']);
		except KeyError:
			params.setIntervalUnit('frames')
		params.setPixelPhysicalSize(acq_metadata_dict['x_physical_size']);
		params.setPixelSizeUnit(acq_metadata_dict['x_unit']);
		params.setMetadataSourceFile(file_path);
	return params;

def import_iq3_metadata(metadata_path):
	"""import basic image metadata based on the metadata saved by iQ3 software at acquisition time"""
	import re
	x_fmt_str = 'x \: (?P<x_pixels>\d+\.?\d*) \* (?P<x_physical_size>\d+\.?\d*) \: (?P<x_unit>\w+)';
	y_fmt_str = 'y \: (?P<y_pixels>\d+\.?\d*) \* (?P<y_physical_size>\d+\.?\d*) \: (?P<y_unit>\w+)';
	z_fmt_str = '\s*Repeat Z \- (?P<z_extent>[+-]?\d+\.?\d*) (?P<z_unit>\w+) in (?P<z_pixels>\d+\.?\d*) planes \(centre\)';
	t_fmt_str = '\s*Repeat T \- (?P<n_frames>\d+\.?\d*) times \((?P<frame_interval>\d+\.?\d*) (?P<time_unit>(min|sec))\)';
	c_fmt_str = r"\s*Repeat \- Channel \((?P<raw_channels_str>\w.*)\)";
	format_strings = [x_fmt_str, y_fmt_str, z_fmt_str, t_fmt_str, c_fmt_str];
	
	meta_dict = {}
	metadata_file = open(metadata_path, 'r')
	try:
		for line in metadata_file.readlines():
			for fmt_str in format_strings:
				m = re.match(fmt_str, line)
				if (bool(m)):
					meta_dict.update(m.groupdict())
		p_num = re.compile('[+-]?\d+\.?\d*')
		try:
			iteritems = meta_dict.iteritems();
		except:
			iteritems = meta_dict.items();
		for key, value in iteritems:
			if p_num.match(value):
				try:
					meta_dict[key] = float(value)
				except:
					#print("conversion error for key " + key);
					continue;
	finally:
		metadata_file.close();
		if 'raw_channels_str' in meta_dict:
			ch_list = meta_dict['raw_channels_str'].split(",")
			meta_dict['n_channels'] = len(ch_list);
			meta_dict['channel_list'] = ch_list;
		return meta_dict