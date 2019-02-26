# Script to generate whole membrane lengths without trying to be smart about where the neck of blebs lie
# D.J. Kelly, 2019-02-26, douglas.kelly@riken.jp

import os, csv, sys
from ij.io import DirectoryChooser

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
import membraneBlebbingFigures as mbfig;

from Parameters import Parameters


dc = DirectoryChooser('Select the root folder from which to search for blebbing code output...');

basefolder = dc.getDirectory();
use_folders = [];

for root, _, files in os.walk(basefolder):
	if "user_defined_edges.zip" in files and "bleb perimeter length.csv" in files and "parameters used.json" in files:
		times = [];
		params = Parameters();
		params.loadParametersFromJson(os.path.join(root, "parameters used.json"));
		f = open(os.path.join(root, "bleb perimeter length.csv"), 'r');
		reader = csv.reader(f);
		for idx, row in enumerate(reader):
			if idx==0:
				headers=[row[0], row[1]]
			else:
				times.append(float(row[0]));
		f.close();
		edges = mbio.load_qcd_edges2(os.path.join(root, "user_defined_edges.zip"));
#		for edge in edges:
#			lengths.append(params.pixel_physical_size * edge.getLength());
		lengths = [params.pixel_physical_size * edge.getLength() for edge in edges];
		print(headers);
		print(times)
		mbio.save_1d_profile_as_csv([(t, l) for t, l in zip(times, lengths)], 
							 os.path.join(root, "full membrane perimeter length.csv"), 
							 [("Time, " + params.interval_unit), "Length, " + params.pixel_unit]);
