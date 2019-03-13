# module containing helper functions for plotting blebbing data
#
# D. J. Kelly, 2018-10-26, douglas.kelly@riken.jp

import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from matplotlib.lines import Line2D
from matplotlib.text import Text
from skimage import io
import pandas
import os
import sys
import re
from textwrap import wrap
import math
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.split(os.getcwd())[0], "classes"))
from Parameters import Parameters

# define string constants:
_um = u'\u00B5m';
_degrees = u'\u00B0';
_squared = u'\u00B2';
_totheminusone = u'\u02C9' + u'\u00B9';
_sigma = u'\u03C3'

def generate_row_titles(paramses):
	"""return the titles for each row of the kymograph/time series montage"""
	return [paramses[0].labeled_species + ' intensity kymographs', 
             'Local curvature kymographs', 
             'Intensity-merged curvature kymographs',
             'Average ' + paramses[0].labeled_species + ' intensity', 
             paramses[0].labeled_species + ' intensity standard deviation', 
             'Ratio +/- curvature intensities', 
             'Bleb area, ' + _um + _squared, 
			 'ML/EL ratio'];

def set_font_sizes(title_size=28, subtitle_size=14, axis_label_size=12, tick_label_size=12):
	"""control all font sizes in the notebook"""
	plt.rcParams.update({'axes.titlesize': subtitle_size});
	plt.rcParams.update({'axes.labelsize': axis_label_size});
	plt.rcParams.update({'xtick.labelsize': tick_label_size});
	plt.rcParams.update({'ytick.labelsize': tick_label_size});
	plt.rcParams.update({'figure.titlesize': title_size});

def generate_average_intensity_plots(actin_kym, curvature_kym, curvature_threshold=0):
    """generate curves showing spatially-averaged actin intensity over time, and ratio of intensities above and below threshold curvature"""
    pos_curv_mask = np.ma.masked_greater(curvature_kym, curvature_threshold).mask;
    zero_intensity_mask = np.ma.masked_less_equal(actin_kym, curvature_threshold).mask;
    neg_curv_actin = np.ma.fix_invalid(actin_kym, 
                                       mask=np.bitwise_or(pos_curv_mask, zero_intensity_mask), 
                                       fill_value=None);
    pos_curv_actin = np.ma.fix_invalid(actin_kym, 
                                       mask=np.bitwise_or(np.bitwise_not(pos_curv_mask), zero_intensity_mask), 
                                       fill_value=None);
    all_actin = np.ma.fix_invalid(actin_kym, 
                                 mask=zero_intensity_mask, 
                                 fill_value=None);
    all_actin_I = all_actin.mean(0) - all_actin.min();
    all_actin_std = all_actin.std(0);
    pos_curv_actin_I = pos_curv_actin.mean(0) - all_actin.min();
    neg_curv_actin_I = neg_curv_actin.mean(0) - all_actin.min();
    return all_actin_I, pos_curv_actin_I, neg_curv_actin_I, all_actin_std;

def load_kymograph_data(experiment_folder, condition_names=None, use_subfolder_as_condition=False, normalise_intensity_to_std=False):
	"""load curvature and intensity kymographs, along with image parameters"""
	actin_kyms = []; 
	curvature_kyms = [];
	paramses = [];
	areas = [];
	ml_el_ratios = [];
	subtitles = [];

	subfolders = [x for x in os.listdir(experiment_folder) if os.path.isdir(os.path.join(experiment_folder, x))];
	num_fmt_str = '\d+\.?\d*';
	split_str = re.split(num_fmt_str, subfolders[0], maxsplit=1);
	if len(split_str) > 1:
		fmt_str = split_str[0] + num_fmt_str + split_str[1];
		order_numerically = all([bool(re.match(fmt_str, subfolder)) for subfolder in subfolders]);
		if order_numerically:
			numbers = [float(re.search(r'\d+\.?\d*', s).group(0)) for s in subfolders];
			idxs = [numbers.index(x) for x in sorted(numbers)]
			subfolders = [subfolders[idx] for idx in idxs];

	for idx, subfolder in enumerate(subfolders):
		if os.path.isdir(os.path.join(experiment_folder, subfolder)):
			params = Parameters();
			params.loadParametersFromJson(os.path.join(experiment_folder, subfolder, "parameters used.json"));
			paramses.append(params);

			actin_kym = io.imread(os.path.join(experiment_folder, subfolder, "normalised position {} kymograph.tif".format(params.labeled_species)));
			curvature_kym = io.imread(os.path.join(experiment_folder, subfolder, "normalised position curvature kymograph.tif"));
			
			if 'physical_curvature_unit' not in dir(params):
				curvature_kym = curvature_kym * (1/paramses[idx].pixel_physical_size);
			elif params.physical_curvature_unit=='':
				curvature_kym = curvature_kym * (1/paramses[idx].pixel_physical_size);
			try:
				areas.append(pd.read_csv(os.path.join(experiment_folder, subfolder, "bleb area.csv")));
			except UnicodeDecodeError as e:
				areas.append(pd.read_csv(os.path.join(experiment_folder, subfolder, "bleb area.csv"), encoding='cp1252'));
			try:
				full_membrane_lengths = pd.read_csv(os.path.join(experiment_folder, subfolder, "full membrane length.csv"));
			except UnicodeDecodeError as e:
				full_membrane_lengths = pd.read_csv(os.path.join(experiment_folder, subfolder, "full membrane length.csv"), encoding='cp1252');
			try:
				full_membrane_euclideans = pd.read_csv(os.path.join(experiment_folder, subfolder, "full membrane euclidean length.csv"));
			except UnicodeDecodeError as e:
				full_membrane_euclideans = pd.read_csv(os.path.join(experiment_folder, subfolder, "full membrane euclidean length.csv"), encoding='cp1252');
			ml_el_ratio = full_membrane_euclideans;
			fml_length_key = [s for s in full_membrane_lengths.columns if "Length" in s][0];
			fme_length_key = [s for s in full_membrane_euclideans.columns if "Length" in s][0];
			ml_el_ratio["Accumulated/Euclidean distance ratio"] = full_membrane_lengths[fml_length_key]/full_membrane_euclideans[fme_length_key];
			ml_el_ratio.drop(fme_length_key, axis=1, inplace=True);
			ml_el_ratios.append(ml_el_ratio);
			if 'qc_background_rois' in dir(params) and normalise_intensity_to_std:
				bg_stds = pd.read_csv(os.path.join(experiment_folder, subfolder, params.labeled_species + " channel background standard deviations.csv"), encoding='cp1252');
				bg_std_median = bg_stds[params.labeled_species + " bg std"].median();
				actin_kym = actin_kym/bg_std_median;
			if condition_names is None:
				if use_subfolder_as_condition:
					subtitles.append(subfolder.replace("um", _um))
				else:
					subtitles.append(os.path.splitext(os.path.basename(params.input_image_path))[0].replace("um", _um));
			else:
				subtitles.append(condition_names[idx].replace("um", _um));
			actin_kyms.append(actin_kym);
			curvature_kyms.append(curvature_kym);
	return actin_kyms, curvature_kyms, paramses, subtitles, areas, ml_el_ratios;  

def curvature_with_intensity(curv_im, actin_im, contrast_enhancement=0):
    """taking imshow results as inputs, return an rgb image of curvature weighted by intensity scaled by a contrast enhancement parameter"""
    curv_range = curv_im.get_clim();
    clrm = curv_im.cmap;
    i_range = actin_im.get_clim();
    curvature_kym = curv_im.get_array();
    actin_kym = actin_im.get_array();
    norm_curv = np.clip((np.round(255*(curvature_kym - min(curv_range))/(max(curv_range) - min(curv_range)))).astype(int), 0, 255);
    norm_i = np.clip((actin_kym - min(i_range))/((max(i_range) - min(i_range)) * (1 - contrast_enhancement)), 0, 1);
    vf = np.vectorize(clrm);
    merged_im = np.stack([rgb_curv * norm_i for rgb_curv in vf(norm_curv)[:3]], 2)
    return merged_im;

def add_scalebar(params, ax, fig, color='w', scale_bar_size_um=1, vertical_scale_bar=False, prenormalised=True):
	"""add a scalebar to the provided axes"""
	if not prenormalised:
		scalebar_len_pix = round(float(scale_bar_size_um) / params.pixel_physical_size);
	else:
		scalebar_len_pix = scale_bar_size_um;
	if vertical_scale_bar:
		axis_height_pt = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted()).width*fig.dpi;
		sb_xdata=[0.95 * max(ax.get_xlim()), 0.95 * max(ax.get_xlim())]
		anchor = min(ax.get_ylim()) + 0.025 * (max(ax.get_ylim()) - min(ax.get_ylim()));
		sb_ydata=[anchor + scalebar_len_pix, anchor];
	else:
		axis_height_pt = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted()).height*fig.dpi;
		sb_xdata=[0.9 * max(ax.get_xlim()) - scalebar_len_pix, 0.9 * max(ax.get_xlim())];
		sb_ydata=[0.95 * max(ax.get_ylim()), 0.95 * max(ax.get_ylim())];
	scalebar = Line2D(color=color, 
						xdata=sb_xdata, 
						ydata=sb_ydata,
						linewidth=axis_height_pt/25);
	ax.add_artist(scalebar);

def plot_kymographs(paramses, actin_kyms, curvature_kyms, axs, intensity_lims, curv_lims, contrast_enhancement=0.35):
	"""from loaded kymograph images, show kymographs in the notebook"""
	for idx in range(len(curvature_kyms)):
		t_offs = 0;
		if paramses[idx].time_crop_start_end is not None:
			if paramses[idx].time_crop_start_end[0] is not None:
				t_offs = float(paramses[idx].time_crop_start_end[0] * paramses[idx].frame_interval);
		actin_im=axs[0][idx].imshow(actin_kyms[idx], 
							   plt.cm.gray, 
							   vmin=intensity_lims[0], 
							   vmax=intensity_lims[1],
							   aspect='auto', 
							   interpolation=None, 
							   extent=[t_offs, 
									  paramses[idx].frame_interval*actin_kyms[idx].shape[1]+t_offs, 
									  -actin_kyms[idx].shape[0] * paramses[idx].pixel_physical_size/2, 
									  actin_kyms[idx].shape[0] * paramses[idx].pixel_physical_size/2]);
		curv_im=axs[1][idx].imshow(curvature_kyms[idx], 
                           plt.cm.jet, 
                           vmin=curv_lims[0], 
                           vmax=curv_lims[1], 
                           aspect='auto', 
                           interpolation=None, 
                           extent=[t_offs,
                                   paramses[idx].frame_interval*actin_kyms[idx].shape[1]+t_offs, 
                                   -actin_kyms[idx].shape[0] * paramses[idx].pixel_physical_size/2, 
                                   actin_kyms[idx].shape[0] * paramses[idx].pixel_physical_size/2]);
		merged_im = curvature_with_intensity(curv_im, actin_im, contrast_enhancement=contrast_enhancement);
		axs[2][idx].imshow(merged_im, 
							aspect='auto', 
							interpolation=None, 
							extent=[t_offs, 
									paramses[idx].frame_interval*actin_kyms[idx].shape[1]+t_offs, 
									-actin_kyms[idx].shape[0] * paramses[idx].pixel_physical_size/2, 
									actin_kyms[idx].shape[0] * paramses[idx].pixel_physical_size/2]);
	return actin_im, curv_im;

def plot_timeseries(paramses, axs, actin_kyms, curvature_kyms, areas, ml_el_ratios, curvature_threshold=0):
	"""show 1D plots of different parameters against time"""
	average_nonzero_intensities = [];
	std_nonzero_intensities = [];
	intensity_ratios = [];
	ml_el_time_variances = [];
	for idx in range(len(actin_kyms)):
		t_offs = 0;
		if paramses[idx].time_crop_start_end is not None:
			if paramses[idx].time_crop_start_end[0] is not None:
				t_offs = float(paramses[idx].time_crop_start_end[0] * paramses[idx].frame_interval);
		t = [t_offs + tidx*paramses[idx].frame_interval for tidx in range(actin_kyms[idx].shape[1])];
		average_nonzero_intensity, avg_nzi_pos_curv, avg_nzi_neg_curv, std_nzi = generate_average_intensity_plots(actin_kyms[idx], curvature_kyms[idx], curvature_threshold=curvature_threshold);
		avgIplt = axs[3][idx].plot(t, average_nonzero_intensity, 'b-');
		varIplt = axs[4][idx].plot(t, std_nzi, 'c-');
		ratioIplt = axs[5][idx].plot(t, avg_nzi_pos_curv/avg_nzi_neg_curv, 'r-');
		axs[5][idx].plot(axs[5][idx].get_xlim(), [1, 1], 'k--');
		areaplt = axs[6][idx].plot(t, areas[idx].iloc[:,1], 'g-');
		mlelplt = axs[7][idx].plot(t, ml_el_ratios[idx]['Accumulated/Euclidean distance ratio'], 'm-');
		ml_el_time_variances.append(ml_el_ratios[idx]['Accumulated/Euclidean distance ratio'].var())
		average_nonzero_intensities.append(average_nonzero_intensity);
		std_nonzero_intensities.append(std_nzi);
		intensity_ratios.append(avg_nzi_pos_curv/avg_nzi_neg_curv);
	return average_nonzero_intensities, std_nonzero_intensities, intensity_ratios, ml_el_time_variances;

def get_normalisation_limits(paramses, actin_kyms, curvature_kyms, make_colorscale_symmetrical=True):
	"""generate limits based on range of membrane lengths, analysis durations, intensities and curvatures for normalising plots"""
	ts = [];
	hs = [];
	for idx in range(len(actin_kyms)):
		t_offs = 0;
		if paramses[idx].time_crop_start_end is not None:
			if paramses[idx].time_crop_start_end[0] is not None:
				t_offs = float(paramses[idx].time_crop_start_end[0] * paramses[idx].frame_interval);
		t = [t_offs + tidx*paramses[idx].frame_interval for tidx in range(actin_kyms[idx].shape[1])];
		ts.append(t)
		hs.append(actin_kyms[idx].shape[0] * paramses[idx].pixel_physical_size);
	space_lims = [-max(hs)/2, max(hs)/2];
	time_lims = [min([min(tt) for tt in ts]), max([max(tt) for tt in ts])];
	curv_lims = (min([im.min() for im in curvature_kyms]), max([im.max() for im in curvature_kyms]));
	if make_colorscale_symmetrical:
		curv_lims = (-max([abs(x) for x in curv_lims]), max([abs(x) for x in curv_lims]))
	intensity_lims = (min([im.min() for im in actin_kyms]), max([im.max() for im in actin_kyms]));
	return space_lims, time_lims, intensity_lims, curv_lims;

def adjust_kymograph_display(fig, axs, paramses, actin_im, curv_im, make_colorscale_symmetrical, column_titles, normalise_intensity_to_std=False):
	"""fine-tune the display properties of kymographs in montages"""
	row_titles = generate_row_titles(paramses);
	for idx in range(axs.shape[1]):
		axs[0][idx].set_yticks([]);
		axs[1][idx].set_yticks([]);
		axs[2][idx].set_yticks([]);
		axs[0][idx].set_title("\n".join(wrap(column_titles[idx].replace("_", " "), int(1200.0/(plt.rcParams['axes.labelsize'] * len(column_titles))))));
	cbar_ax1 = fig.add_axes([0.95, axs[0][axs.shape[1]-1].get_position().y0, 
							0.03, (axs[0][axs.shape[1]-1].get_position().y1 - axs[0][axs.shape[1]-1].get_position().y0)])
	cbar_ax2 = fig.add_axes([0.95, axs[1][axs.shape[1]-1].get_position().y0, 
								0.03, (axs[1][axs.shape[1]-1].get_position().y1 - axs[1][axs.shape[1]-1].get_position().y0)])
	fig.colorbar(actin_im, cax=cbar_ax1);
	fig.colorbar(curv_im, cax=cbar_ax2);
	intensity_cbar_label = paramses[0].labeled_species + " intensity, A.U.";
	if all(['qc_background_rois' in dir(params) for params in paramses]) and normalise_intensity_to_std:
		intensity_cbar_label = paramses[0].labeled_species + " intensity normalised to background " + _sigma + ", A.U.";
	cbar_ax1.set_ylabel("\n".join(wrap(intensity_cbar_label, 30)));
	cbar_ax2.set_ylabel("Local curvature (" + paramses[0].pixel_unit + _totheminusone + ')');
	for row_idx in range(3):
		bbox = axs[row_idx][0].get_position();
		ypos = 0.5*(bbox.ymin + bbox.ymax);
		txt = Text(0.05, ypos, "\n".join(wrap(row_titles[row_idx], int(300.0/plt.rcParams['axes.labelsize']))), 
					rotation='vertical', 
					horizontalalignment='center', 
					va='center', 
					fontsize=plt.rcParams['axes.labelsize']);
		fig.add_artist(txt);

def adjust_time_series_plot_display(fig, axs, paramses, average_nonzero_intensities, std_nonzero_intensities, intensity_ratios, areas, ml_el_ratios):
	"""fine-tune the display properties of time series in montages"""
	row_titles = generate_row_titles(paramses);
	avgd_intensity_lims = [10*np.floor(np.min([i.min() for i in average_nonzero_intensities])/10), 
                       10*np.ceil(max([i.max() for i in average_nonzero_intensities])/10)];
	std_intensity_lims = [10*np.floor(np.min([i.min() for i in std_nonzero_intensities])/10), 
						   10*np.ceil(max([i.max() for i in std_nonzero_intensities])/10)];
	ratio_lims = [np.min([r.min() for r in intensity_ratios]), 
				  np.max([r.max() for r in intensity_ratios])];
	area_lims = [min([a.iloc[:,1].min() for a in areas]), max([a.iloc[:,1].max() for a in areas])];
	ml_el_ratio_lims = [min([aer.iloc[:,1].min() for aer in ml_el_ratios]), max([aer.iloc[:,1].max() for aer in ml_el_ratios])];
	n_y_ticks = 4;
	for idx in range(axs.shape[1]):
		axs[3][idx].set_yticks([int(avgd_intensity_lims[0] + yidx*(avgd_intensity_lims[1] - avgd_intensity_lims[0])/(n_y_ticks-1)) 
							for yidx in range(n_y_ticks)]);
		if idx!=0:
			axs[3][idx].set_yticks([]);
			axs[4][idx].set_yticks([]);
			axs[5][idx].set_yticks([]);
			axs[6][idx].set_yticks([]);
		axs[3][idx].set_ylim(avgd_intensity_lims);
		axs[4][idx].set_ylim(std_intensity_lims)
		axs[5][idx].set_ylim(ratio_lims);
		axs[6][idx].set_ylim(area_lims);
		axs[7][idx].set_ylim([1, ml_el_ratio_lims[1]]);
		axs[7][idx].set_xlabel("Time, " + paramses[idx].interval_unit);
	for row_idx in range(5):
		axs[row_idx+3][0].set_ylabel("\n".join(wrap(row_titles[row_idx+3], int(200.0/plt.rcParams['axes.labelsize']))))

def do_space_time_normalisation(axs, paramses, space_lims, time_lims, make_colorscale_symmetrical):
	"""ensure that time and space coordinates are displayed consistently across the dataset"""
	n_x_ticks = 4;
	for idx in range(axs.shape[1]):
		for iidx in range(len(generate_row_titles(paramses))):
			axs[iidx][idx].set_xlim(time_lims)
		axs[6][idx].set_xticks([time_lims[0] + (xidx*time_lims[1]/(n_x_ticks-1)) for xidx in range(n_x_ticks)]);
		axs[5][idx].plot(time_lims, [1, 1], 'k--');

		axs[0][idx].set_ylim(space_lims)
		axs[1][idx].set_ylim(space_lims)
		axs[2][idx].set_ylim(space_lims)

		axs[0][idx].set_facecolor('k');
		axs[2][idx].set_facecolor('k');
		if make_colorscale_symmetrical:
			axs[1][idx].set_facecolor(cm.jet(128));
	return;

def add_arrow_annotations(axs, annotation_arrows):
	"""Add arrows to annotate points on intensity and curvature kymographs"""
	if annotation_arrows is not None:
		for ann in annotation_arrows:
			direction = -1 if ann[1][1] > 0.5 else 1;
			axs[0][ann[0]].annotate('', xycoords='axes fraction', xy=ann[1], 
									xytext=(ann[1][0], ann[1][1] - direction * 0.1), textcoords='axes fraction',
									arrowprops=dict(facecolor='black', edgecolor='white'));
			axs[1][ann[0]].annotate('', xycoords='axes fraction', xy=ann[1], 
									xytext=(ann[1][0], ann[1][1] - direction * 0.1), textcoords='axes fraction',
									arrowprops=dict(facecolor='black', edgecolor='white'));
	return;

def adjust_overlay_images_display(fig, axs, paramses, curv_lims, experiment_title, row_titles):
	"""fine-tune the display properties of overlay images in montages"""
	plt.tight_layout();
	fig.subplots_adjust(top=0.9, right=0.9, hspace=0.15, left=0.1);
	cmap = plt.get_cmap('jet');
	norm = colors.Normalize(vmin=curv_lims[0],vmax=curv_lims[1]);
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm);
	sm.set_array([]);
	cbar_ax1 = fig.add_axes([0.95, 
								axs[axs.shape[0]-1][axs.shape[1]-1].get_position().y0, 
								0.03, 
								(axs[0][axs.shape[1]-1].get_position().y1 - axs[axs.shape[0]-1][axs.shape[1]-1].get_position().y0)])
	fig.colorbar(sm, cax=cbar_ax1);
	cbar_ax1.set_ylabel("Local curvature (" + paramses[0].pixel_unit + _totheminusone + ")")
	plt.rc('figure', titleweight='bold');
	fig.suptitle(experiment_title + " curvature overlays");
	for rowidx in range(axs.shape[0]):
		bbox = axs[rowidx][0].get_position();
		txt = Text(0.05, 0.5*(bbox.ymin + bbox.ymax), 
					"\n".join(wrap(row_titles[rowidx].replace("_", " ").replace("um", _um), int(800.0/(plt.rcParams['axes.labelsize'] * axs.shape[0])))), 
					rotation='vertical', 
					horizontalalignment='center', 
					va='center');
		fig.add_artist(txt);

def plot_overlay_images(paramses, fig, axs, ims, raw_curvatures, curv_lims, base_chosen_timepoints, scale_bar_color='w', scale_bar_size_um=1):
	"""overlay colormapped curvatures on cell membrane-identification images and tile"""
	for imidx in range(len(ims)):
		t_offs = 0;
		if paramses[imidx].time_crop_start_end is not None:
			if paramses[imidx].time_crop_start_end[0] is not None:
				t_offs = float(paramses[imidx].time_crop_start_end[0]) * paramses[imidx].frame_interval;
		chosen_timepoints = [x + paramses[imidx].frame_interval for x in base_chosen_timepoints];
		chosen_frames = [round((t - t_offs)/ paramses[imidx].frame_interval) for t in chosen_timepoints];
		closest_timepoints = [(fr - 1) * paramses[imidx].frame_interval + t_offs for fr in chosen_frames];
		chosen_frames = [ch if ((ch < ims[imidx].shape[0]) and (ch >= 0)) else float('nan') for ch in chosen_frames];
		for fridx, fr in enumerate(chosen_frames):
			cbw = math.ceil(0.95 * ims[imidx][0].shape[1]);
			if not math.isnan(fr):
				curv_pix = np.nonzero(raw_curvatures[imidx][fr]);
				for pixx, pixy in zip(curv_pix[0], curv_pix[1]):
					curv = raw_curvatures[imidx][fr][pixx][pixy];
					existing_val = ims[imidx][fr][pixx][pixy];
					curv_idx_8bit = int(round(255 * (curv - curv_lims[0]) / (curv_lims[1] - curv_lims[0])));
					r, g, b, _ = cm.jet(curv_idx_8bit);
					ims[imidx][fr][pixx][pixy][0] = int(round(255*r));
					ims[imidx][fr][pixx][pixy][1] = int(round(255*g));
					ims[imidx][fr][pixx][pixy][2] = int(round(255*b));
				axs[imidx][fridx].imshow(ims[imidx][fr][:, :cbw, :]);
				axs[imidx][fridx].set_title("t = " + str(closest_timepoints[fridx]) + " s");
			else:
				axs[imidx][fridx].imshow(ims[imidx][0][:, :cbw, :], alpha=0);
				axs[imidx][fridx].set_frame_on(False);
				axs[imidx][fridx].set_visible(False);
			axs[imidx][fridx].set_xticks([]);
			axs[imidx][fridx].set_yticks([]);
			try: # todo: find out why this throws an error in at least one case and deal with properly rather than just (essentially) ignoring
				if fridx == chosen_frames.index(max([f for f in chosen_frames if not np.isnan(f)])):
					add_scalebar(paramses[imidx], 
								axs[imidx][fridx], 
								fig, 
								color=scale_bar_color, 
								scale_bar_size_um=scale_bar_size_um);
					sb_ydata=[0.95 * ims[imidx][0].shape[0], 0.95 * ims[imidx][0].shape[0]]; 
			except:
				pass;
	return;

# for now, const color scale only. Add ability to toggle consistent length scale later?
def get_overlay_normalisation_limits(paramses, ims, raw_curvatures, make_colorscale_symmetrical=True):
	"""return limits for plotting overlays with constant color/length scales"""
	curv_lims =  (min([raw_cs.min() for raw_cs in raw_curvatures]), max([raw_cs.max() for raw_cs in raw_curvatures]))
	if make_colorscale_symmetrical:
		curv_lims = (-max([abs(x) for x in curv_lims]), max([abs(x) for x in curv_lims]));
	return curv_lims;

def load_data_for_overlay_montage(experiment_folder, 
									condition_names, 
									use_subfolder_as_condition):
	"""load images and raw curvature data for montages of overlays"""
	ims = [];
	raw_curvatures = [];
	paramses = [];
	row_titles = [];

	subfolders = [x for x in os.listdir(experiment_folder) if os.path.isdir(os.path.join(experiment_folder, x))];
	for idx, subfolder in enumerate(subfolders):
		params = Parameters();
		params.loadParametersFromJson(os.path.join(experiment_folder, subfolder, "parameters used.json"));
		paramses.append(params);
		im = io.imread(os.path.join(experiment_folder, subfolder, "overlaid curvature.tif"));
		raw_curvature = io.imread(os.path.join(experiment_folder, subfolder, "raw curvature.tif"));
		if 'physical_curvature_unit' not in dir(params):
				raw_curvature = raw_curvature * (1/paramses[idx].pixel_physical_size);
		elif params.physical_curvature_unit=='':
			raw_curvature = raw_curvature * (1/paramses[idx].pixel_physical_size);
		if condition_names is None:
			if use_subfolder_as_condition:
				row_titles.append(subfolder)
			else:
				row_titles.append(os.path.splitext(os.path.basename(params.input_image_path))[0]);
		else:
			row_titles.append(condition_names[idx]);
		ims.append(im);
		raw_curvatures.append(raw_curvature);
	return ims, raw_curvatures, paramses, row_titles;

def plot_ml_el_ratio_variances(ml_el_time_variances, ax, conditions=None, log_yscale=False):
	"""plot a bar chart showing the time-variances of membrane length:euclidean length ratio"""
	if conditions is None:
		conditions = [x+1 for x in range(len(ml_el_time_variances))];
	bar_container = ax.bar([x+1 for x in range(len(ml_el_time_variances))], ml_el_time_variances, tick_label=conditions);
	if log_yscale:
		ax.set_yscale("log");
	else:
		ax.set_yscale("linear");
	ax.set_title("Membrane length/euclidean length time variance");
	ax.set_ylabel("ML/EL time variance");

	for rect in bar_container:
		height = rect.get_height();
		if height < 0.01:
			fmt_str = "{0:.2e}";
		else:
			fmt_str = "{0:.3f}";
		ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
				fmt_str.format(height),
				ha='center', va='bottom')
	plt.tight_layout()