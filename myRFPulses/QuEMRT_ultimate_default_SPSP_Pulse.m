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

% Reset SS package globals
%
clear
ss_opt([]);
ss_globals;
clc

% GENERAL PULSE PARAMETERS
SS_G_TYPE = 'EP Whole';  % Echo-planar design
SS_RF_TYPE = 'ex';  % excitation pulse
SS_OPTs = ss_opt({...
    'Nucleus', 'Carbon', ...
    'Max Duration', 25e-3, ...
    'Sample Time', 10e-6,...
    'Max Grad', 5 ,...
    'Max Slew', 20});

% force pulse design to optimize for center of frequency specification
f_center_Hz = 0;  

% SPECTRAL PULSE PARAMETERS 
f_B0_Gauss = 3e4; % G
f_delta_Hz = 1e-6 * f_B0_Gauss * SS_GAMMA; % 1 ppm 

% metabolite			    frequency (Hz)		freq bandwidth (Hz)		flip angle (deg)    allowed ripple
mets(1).name = 'passband'; 	mets(1).f = 0;      mets(1).df = 10; 		mets(1).ang = 90;   mets(1).d = 0.01;
mets(2).name = 'stopband'; 	mets(2).f = 250; 	mets(2).df = 175;       mets(2).ang = 0;    mets(2).d = 0.01;
mets(3).name = 'stopband'; 	mets(3).f = -250; 	mets(3).df = 175; 		mets(3).ang = 0;    mets(3).d = 0.01;

% create vectors of angles, ripples, and band edges for input to pulse design
[f_band_edges_Hz, f_band_alpha_Rad, f_band_ripple] = create_freq_specs(mets, f_center_Hz);
SS_FILTER_TYPE = 'min';  % minimimum-phase spectral filter

% SPATIAL PULSE PARAMETERS
z_filter_type = 'ls';  % least-squares filter design
z_thickness_cm = 1;  % thickness (cm)
z_timebandwidth = 5; % time-bandwidth, proportional to profile sharpness
z_pass_ripple = 0.01;  
z_stop_ripple = 0.01;  % slice profile pass and stop-band ripples, respectively

% DESIGN THE PULSE
[shape_grad_Gausscm, shape_rf_Gauss, f_sampling_Hz] = ...
    ss_design(z_thickness_cm, z_timebandwidth, [z_pass_ripple, z_stop_ripple], ...
    f_band_edges_Hz, f_band_alpha_Rad, f_band_ripple, SS_RF_TYPE, ...
    z_filter_type, SS_FILTER_TYPE, SS_G_TYPE, f_center_Hz);

set(gcf,'Name', '13C Ultimate SPSP Pulse');


% Pulse Name
filename = '13C_Ultimate_SPSP_Pulse_QuEMRT';
fprintf('\n\nFilename: %s', filename);

% Pulse description
desc = sprintf('SpSp %s at %dT, ',...
    SS_NUCLEUS,...
    f_B0_Gauss * 1e-4);
for m = 1:numel(mets)
    desc = append(desc, sprintf('%s (%dHz): %ddeg, ',...
        mets(m).name, mets(m).f, mets(m).ang));
end
clear('m')
desc = append(desc, sprintf('slice thickness: %dmm, ', ...
    round(10*z_thickness_cm)));
current_time = round(clock);
desc = append(desc, sprintf('date: %d-%02d-%02d', ...
    current_time(1:3)));
fprintf("\nDescription: %s\n", desc);

creation_time = current_time;
clear('current_time')

% Transform and Normalize Shape Units
shape_rf_uT = 100 * shape_rf_Gauss;
shape_rf_max_uT = max(abs(shape_rf_uT));
shape_rf_normalized = shape_rf_uT / shape_rf_max_uT;

shape_rf_shapeInt = sqrt( sum(real(shape_rf_normalized))^2 + sum(imag(shape_rf_normalized))^2) / numel(shape_rf_normalized);
shape_rf_powerInt = shape_rf_normalized * shape_rf_normalized' / numel(shape_rf_normalized);

shape_numel = numel(shape_rf_normalized);
shape_duration_us =  (shape_numel-1) * SS_TS ; % [us]  % DOUBLE CHECK
shape_timesamples_us = (0:shape_numel-1) * SS_TS ; % [us]

shape_rf_normalized = shape_rf_normalized * 100;
shape_rf_magnitude = abs(shape_rf_normalized);
shape_rf_phase_degree = 180 + angle(shape_rf_normalized) * 180/pi;
shape_grad_mTm = shape_grad_Gausscm * 10;

%% Save fig and mat file
save(filename)
savefig(filename)

%% create RFstruct
create_RF_struct
save([RFstruct.filename '_RF_struct'])

%% export to ParaVision
export_RF_PV

%% write JSON file
export_RF_json
