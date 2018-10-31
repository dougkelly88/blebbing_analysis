# functions for handling file input/output for membrane blebbing analysis
#
# D. J. Kelly, 2018-10-15, douglas.kelly@riken.jp

# imports 
import csv, json, os
from datetime import datetime
from ij.io import OpenDialog, DirectoryChooser

# choose folder locations and prepare output folder
def file_location_chooser(default_directory):
	# input
	od = OpenDialog('Choose original file...', 
					default_directory, 
					'*.tif');
	file_path = od.getPath();
	if file_path is None:
		raise IOError('no input file chosen');
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

# save profiles with 2d independent variable, e.g. curvature profile
def save_profile_as_csv(profiles, file_path, data_name):
	f = open(file_path, 'wb');
	writer = csv.writer(f);
	writer.writerow(['Frame', 'x', 'y', data_name]);
	for idx, profile in enumerate(profiles):
		for ((x, y), p) in profile:
			writer.writerow([idx, x, y, p]);
	f.close();

# load profiles with 2d independent variable from csv
def load_csv_as_profile(file_path):
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
			 

# save profiles with 1d independent variable, e.g. length against time
def save_1d_profile_as_csv(profile, file_path, column_names):
	f = open(file_path, 'wb');
	writer = csv.writer(f);
	writer.writerow([column_names[0], column_names[1]]);
	for p in profile:
		writer.writerow([p[0], p[1]]);
	f.close();

# save parameters used for this analysis
def save_parameters(params, file_path):
	f = open(file_path, 'w');
	json.dump(params, f);
	f.close();

#p = load_csv_as_profile("D:\\data\\Inverse blebbing\\output\\2018-10-17 14-56-47 output\\curvatures.csv")
#print(p[-1])