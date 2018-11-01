Version 1.0.1

Installation:

Copy the "blebbing analysis" folder into ..\Fiji.app\plugins\Scripts\Plugins

Usage: 

Run software from ImageJ menu Plugins -> blebbing analysis -> membrane blebbing

*Curvature length parameter (pix): the separation along the membrane between the point at which curvature is to be calculated and the two points on either side from which curvature is calculated using the three-point method and SSS theorem (see http://mathworld.wolfram.com/SSSTheorem.html) - defaults to 10. 

*Threshold method: the (in built ImageJ) thresholding method used to separate membrane from background in the intensity (actin) channel - defaults to "Moments". 

*Curvature overlay LUT: the LUT to use to colorize the curvature when overlaid on the intensity (actin) image - defaults to "physics". 

*Curvature kymograph LUT: the LUT to use for the kymograph showing curvature along the length of the bleb over time - defaults to "Yellow". 

*Labelled species kymograph LUT: the LUT to use for the kymograph showing labelled species intensity along the length of the bleb over time - defaults to "Cyan". 

*Labelled species for intensity analysis: descriptor of the labelled species - defaults to "Actin".
 
*Filter out negative curvatures: coerces negative curvatures to zero in the raw data and does not display them on the overlay. 

*Perform spatial cropping?: toggle whether the software prompts for an ROI at the start of the analysis run. 

*Perform time cropping?: toggle whether the software prompts for start and end frames for a period of interest at the start of the analysis run. 

*Close images on completion?: toggle whether all input/output images involved in the analysis should be closed at the end of the run (after automatically saving to output folder). 

