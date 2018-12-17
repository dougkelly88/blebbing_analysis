# display contents of readme file as a help dialog
#
# D. J. Kelly, 2018-10-16, douglas.kelly@riken.jp

import math, os, sys
from ij.gui import NonBlockingGenericDialog
from javax.swing import JTextArea, JScrollPane
from java.awt import Panel, Dimension
from java.lang import StringBuilder

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

sb = StringBuilder();
for line in text:
	sb.append(line);

panel = Panel();

txtArea = JTextArea(sb.toString());
txtArea.setEditable(False);
txtArea.setLineWrap(True);
txtArea.setWrapStyleWord(True);
scrollpane = JScrollPane(txtArea);
scrollpane.setPreferredSize(Dimension(500,200));
panel.add(scrollpane)
dialog = NonBlockingGenericDialog(title);
#for line in text:
#	dialog.addMessage(line);
dialog.addPanel(panel);
dialog.showDialog();