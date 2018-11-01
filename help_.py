# display contents of readme file as a help dialog
#
# D. J. Kelly, 2018-10-16, douglas.kelly@riken.jp

import math, os, sys
from ij.gui import NonBlockingGenericDialog

script_path = os.path.dirname(os.path.realpath(__file__));
if "Fiji.app" in script_path:
	ss = script_path.split("Fiji.app");
	final_folder = os.path.basename(script_path);
	script_path = os.path.join(ss[0], "Fiji.app", "plugins", "Scripts", "Plugins", final_folder);
sys.path.insert(0, os.path.join(script_path, 'modules'));
sys.path.insert(0, os.path.join(script_path, 'classes'));

from Parameters import Parameters

readme_fpath = os.path.join(script_path, "README.txt");
params = Parameters();

title = "Membrane Blebbing version " + Parameters._version_string;

try:
	f = open(readme_fpath, "rb");
	text = f.readlines();
except:
	raise IOError("Error reading README.txt");
finally:
	f.close();
	
dialog = NonBlockingGenericDialog(title);
for line in text:
	dialog.addMessage(line);
dialog.showDialog();