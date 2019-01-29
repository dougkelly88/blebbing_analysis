Installation:
Copy the "blebbing analysis" folder into ..\Fiji.app\plugins\Scripts\Plugins

Usage: 
Run software from ImageJ menu Plugins -> blebbing analysis -> membrane blebbing

*Curvature length parameter (um): the separation along the membrane between the point at which curvature is to be calculated and the two points on either side. From these three points, curvature is calculated using the three-point method and SSS theorem (see http://mathworld.wolfram.com/SSSTheorem.html). Defaults to 1.0 um. 

*Width of region for intensity analysis (um): the width in microns of the region from which the maximum intensity be extracted; i.e. the software looks 0.5x this width to either side of the identified membrane and assigns the intensity at each position to be the maximum value within the region. Defaults to 0.33 um. 

*Threshold method: the (in built ImageJ) thresholding method used to separate membrane from background in the intensity (actin) channel. Note that these methods are global threshold settings unless marked as "Local". Defaults to "Moments". 

*Curvature overlay LUT: the LUT to use to colorize the curvature when overlaid on the intensity (actin) image. Defaults to "physics". 

*Curvature kymograph LUT: the LUT to use for the kymograph showing curvature along the length of the bleb over time. Defaults to "Yellow". 

*Labelled species kymograph LUT: the LUT to use for the kymograph showing labelled species intensity along the length of the bleb over time. Defaults to "Cyan". 

*Labelled species for intensity analysis: descriptor of the labelled species. Defaults to "Actin".

*Use intensity channel for segmentation too?: Performs segmentation and intensity analysis on the same channel. Not generally recommended as segmentation results may skew intensity measurements. 

*Metadata source: choose the source of image metadata (particularly pixels-micron and frames-time conversion factors); if "Image metadata" is selected but can't be extracted, software will prompt for an iQ3 acquisition metadata file instead. 

*Constrain anchors close to manual selections?: Prevent software from evolving the initial position of the anchors (before fixing to the membrane edge) at each frame. Select if there's no significant movement/drift in the unblebbed membrane. 

*Filter out negative curvatures: coerces negative curvatures to zero in the raw data and does not display them on the overlay. 

*Account for photobleaching?: toggle a correction which ensures constant intensity averaged over frames across the imaging sequence, c.f. ratio method in ImageJ. 

*Perform quality control of membrane edges?: toggle whether the software will prompt for user intervention to review membrane edges determined by the automatic segmentation process, and modify membrane edges by drawing freehand lines if necessary. 

*Perform spatial cropping?: toggle whether the software prompts for an ROI at the start of the analysis run. 

*Perform time cropping?: toggle whether the software prompts for start and end frames for a period of interest at the start of the analysis run. 

*Close images on completion?: toggle whether all input/output images involved in the analysis should be closed at the end of the run (after automatically saving to output folder). 

*Compare inner and outer curvature regions?: If toggled on, loops through edge identification steps twice and saves ratios of inner and outer curvature values, as well as information on the variation of ratios. Results can subsequently be plotted using the scripts in ..\blebbing_analysis\plotting\


To re-run analysis on previously analysed data, bypassing the membrane edge definition steps: 
Run software from ImageJ menu Plugins -> blebbing analysis -> re-run blebbing analysis
When prompted, select the folder containing the previous analysis output
If the original image can't be located, software will also prompt for this
Parameters will be set according to the last analysis; change to taste, or simply run again to get an analysis run with the latest software version. 

To loop analysis over a series of parameter values:
Prepare a text file with a single line in this format: {"parameter name": [value1, value2, value3, ...]}
(e.g. {"intensity_profile_width_um": [0.3, 0.5, 1.0, 2.0]}, or see file example_loop_param.txt)
Run software from ImageJ menu Plugins -> blebbing analysis -> looped blebbing analysis
When prompted, select the folder containing the previous analysis output
If the original image can't be located, software will also prompt for this
Software will then loop over the parameters and perform the analysis as required, saving into the folder above the original output in the form ".\<dated folder>\looped analysis\<loop parameter> = <value>"

To re-draw edges for later re-analysis:
Run software from ImageJ menu Plugins -> blebbing analysis -> view and edit edges
When prompted, select the folder containing the previous analysis output
If the original image can't be located, software will also prompt for this
Then proceed as normal with membrane quality control. On completion, an image will be be displayed showing the resulting membrane edges overlaid on the image, and membrane edges will be saved to a new (dated) output folder in the same parent folder as the prior data. 