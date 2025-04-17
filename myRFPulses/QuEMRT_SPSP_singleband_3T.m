%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectral-Spatial RF Pulse Design for MRI and MRSI MATLAB Package
%
% Authors: Adam B. Kerr and Peder E. Z. Larson
%
% (c)2007-2014 Board of Trustees, Leland Stanford Junior University and
%	The Regents of the University of California. 
% All Rights Reserved.
%
% Please see the Copyright_Information and README files included with this
% package.  All works derived from this package must be properly cited.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reset SS package globals
clear
ss_globals;
clc

%% Define name and field strength

% name used for files and plots
filename = 'QuEMRT_SPSP_singleband_3T';

% define field strength
b0_tesla = 3; % B0 field strenth in T


%% Passband and stopband settings

% choose passband settings
pass_f_ppm = 0.0; % passband center freq (ppm)
pass_df_ppm = 0.5; % passband width (ppm)

% choose stopband settings
stop_f_ppm = 10.0; % stopband center freq (ppm)
stop_df_ppm = 4.0; % stopband width (ppm)

%% SPSP settings

% general
SS_OPTs = ss_opt({'Nucleus', 'Carbon', ...
	      'Max Duration', 30e-3, ...
          'Sample Time', 4e-6,...
	      'Num Lobe Iters', 10, ...
	      'Max B1', 0.2, ...
	      'Num Fs Test', 100, ...
	      'Verse Fraction', 0.90, ...
	      'SLR', 0, ...
	      'B1 Verse', 0, ...
	      'Min Order', 1,...
	      'Spect Correct', 1,...
          'Max Grad',4.0,...
          'Max Slew',12.0});

SS_FILTER_TYPE = 'min';         % spectral filter [NOT USED]
SS_G_TYPE = 'Flyback Whole';    % SPSP type
SS_RF_TYPE = 'ex';              % RF type

z_filter_type = 'ls';             % least-squares filter design
z_thickness_cm = 1.25;               % slice thickness (cm)
z_timebandwidth = 3.0;                 % time-bandwidth product of RF pulse
z_pass_ripple = 0.01;       % slice profile pass ripple
z_stop_ripple = 0.01;       % slice profile stop ripple

%% Calculate values and create mets struct

f_B0_Gauss = b0_tesla * 1e4; % convert from Tesla to Gauss
f_delta_Hz = 1e-6 * f_B0_Gauss * SS_GAMMA; % 1 ppm in Hz 

% convert values to Hz and fill mets struct
pass_f_hz = pass_f_ppm * f_delta_Hz;
pass_df_hz = pass_df_ppm * f_delta_Hz;
stop_f_hz = stop_f_ppm * f_delta_Hz;
stop_df_hz = stop_df_ppm * f_delta_Hz;

% metabolite			    frequency (Hz)		        freq bandwidth (Hz)		    flip angle (deg)    allowed ripple
mets(1).name = 'stopband'; 	mets(1).f = -stop_f_hz; 	mets(1).df = stop_df_hz; 	mets(1).ang = 0;    mets(1).d = 0.005;
mets(2).name = 'passband'; 	mets(2).f = pass_f_hz;      mets(2).df = pass_df_hz; 	mets(2).ang = 90;   mets(2).d = 0.01;
mets(3).name = 'stopband'; 	mets(3).f = stop_f_hz; 	    mets(3).df = stop_df_hz;    mets(3).ang = 0;    mets(3).d = 0.005;

% create vectors of angles, ripples, and band edges for input to pulse design
f_center_Hz = 0; % force pulse design to optimize for center of frequency specification
[f_band_edges_Hz, f_band_alpha_Rad, f_band_ripple] = create_freq_specs(mets, f_center_Hz);

clear("pass_f_hz", "pass_df_hz", "stop_f_hz", "stop_df_hz", "pass_f_ppm", "pass_df_ppm", "stop_f_ppm", "stop_df_ppm")

%% DESIGN THE PULSE
% choose 3) Fs:  530.8 B1: 0.127G Power: 1.854e-05 G^2 ms Dur: 26.1ms

[shape_grad_Gausscm, shape_rf_Gauss, f_sampling_Hz] = ...
    ss_design(z_thickness_cm, z_timebandwidth, [z_pass_ripple, z_stop_ripple], ...
    f_band_edges_Hz, f_band_alpha_Rad, f_band_ripple, SS_RF_TYPE, ...
    z_filter_type, [], SS_G_TYPE, f_center_Hz, 0);
set(gcf,'Name', filename);

%% interpolate shapes to 10µs raster time
new_raster_time = 10e-6;

num_samples = numel(shape_rf_Gauss);
total_dur = (num_samples - 1) * SS_TS;
num_samples_new = floor(total_dur / new_raster_time);

time_vector_old = 0:SS_TS:(num_samples - 1) * SS_TS;
time_vector_new = 0:new_raster_time:(num_samples_new - 1) * new_raster_time;

shape_grad_interp = interp1(time_vector_old, shape_grad_Gausscm, time_vector_new, 'linear');
shape_rf_interp = interp1(time_vector_old, shape_rf_Gauss, time_vector_new, 'linear');

clear("num_samples", "total_dur", "num_samples_new", "time_vector_old", "time_vector_new")

% overwrite data with new interpolated ones
SS_TS = new_raster_time;
shape_grad_Gausscm = shape_grad_interp;
shape_rf_Gauss = shape_rf_interp;

clear("new_raster_time", "shape_rf_interp", "shape_grad_interp")

%% plot results for individual metabolites
mets_freq_ppm = [12.6 0 -9.7];
mets_freq_hz = mets_freq_ppm * f_delta_Hz;
mets_names = {'Lactate', 'Pyruvate','Bicarbonate'};

for idx = 1:length(mets_freq_hz)
    rf_shift = ss_shift(shape_grad_Gausscm,shape_rf_Gauss,0,mets_freq_hz(idx));
    ss_plot(shape_grad_Gausscm, rf_shift, SS_TS, SS_RF_TYPE, z_thickness_cm*3, 2.5*[min(mets_freq_hz) max(mets_freq_hz)], SS_GAMMA, mets_freq_hz);
    set(gcf,'Name',[mets_names{idx} ' Excitation'],'NumberTitle','Off')

end

clear("mets_freq_ppm", "mets_freq_hz", "mets_names", "idx", "rf_shift")

%% Create output folder

% Define the name of the main output folder
outputFolderName = 'output';

% Get the current date in YYYY_mm_dd format using datetime
currentDate = string(datetime('now', 'Format', 'yyyy_MM_dd'));

% Construct the name of the subfolder
subfolderName = sprintf('%s_%s', currentDate, filename);

% Construct the full path to the output folder
outputPath = fullfile(pwd, outputFolderName);

% Construct the full path to the subfolder
subfolderPath = fullfile(outputPath, subfolderName);

% Check if the output folder exists, and create it if it doesn't
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
    disp(['Created output folder: ', outputPath]);
end

% Check if the subfolder exists, and create it if it doesn't
if ~exist(subfolderPath, 'dir')
    mkdir(subfolderPath);
    disp(['Created subfolder: ', subfolderPath]);
else
    disp(['Subfolder already exists: ', subfolderPath]);
end

clear("currentDate", "subfolderName", "outputFolderName", "outputPath")

%% Save complete mat file
save(fullfile(subfolderPath, [filename, '.mat']))

%% save as VARIAN RF pulse files for Pulseq X-EPI sequence
ss_save(shape_grad_Gausscm, ...
        shape_rf_Gauss, ...
        max(f_band_alpha_Rad), ...
        z_thickness_cm,...
        [], ...
        'Varian', ...
        f_band_edges_Hz, ...
        max(f_band_alpha_Rad),...
        fullfile(subfolderPath, filename));

%% create RFstruct and save workspace
create_RF_struct
save(fullfile(subfolderPath, [RFstruct.filename '_RF_struct']))

%% export to ParaVision
export_RF_PV(RFstruct, subfolderPath)

%% write JSON file
export_RF_json(RFstruct, subfolderPath)