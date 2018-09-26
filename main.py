# port program written by Dr. Satoru Okuda for curvature analysis to ImageJ
# Motivation: let ImageJ do heavy lifting on UI, image processing, file IO side, leaving code 
# that does the clever stuff less obscured. Also makes cross-platform deployment easier...
# D.J. Kelly, 2018-09-26, douglas.kelly@riken.jp

# python (jython) imports
import os, sys

# java imports - aim to have UI components entirely swing, listeners and layouts awt
from java.awt import Dimension, GridBagLayout, GridBagConstraints, GridLayout
import javax.swing as swing
import javax.swing.table.TableModel

# imagej imports
from ij.gui import NonBlockingGenericDialog
from ij.io import DirectoryChooser

def main():
	#print (sys.version_info) # debug
	#print(sys.path) # debug

# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main()