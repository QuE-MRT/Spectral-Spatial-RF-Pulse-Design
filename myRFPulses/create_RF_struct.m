% MATLAB Script: Organize RF Pulse Variables into a Single RF-Struct Object
%
% Description:
% This script collects multiple RF pulse-related variables from the MATLAB 
% workspace and organizes them into a single structured object (RF-struct). 
% The RF-struct consolidates all relevant parameters, simplifying further 
% processing, storage, or analysis of RF pulse data.
%
% Features:
% - Detects and retrieves prior specified RF pulse-related variables.
% - Organizes variables into a structured format for easy access and manipulation.
% - Ensures consistency and reduces redundancy in data handling.
%
% Author: Christoph A. Müller
% Date: 17.01.2025
% Version: 1.0
%
% Usage:
% 1. Run pulse generation script first
% 2. Run the script to consolidate known variables into an RF-struct object.
% 3. The RF-struct can be saved, analyzed, or exported as needed.

% Create pulse description
desc = filename;

% Save creation time
creation_time = round(clock);

% Determine pulse properties
shape_num_samples = numel(shape_rf_Gauss); % total number of grad / rf samples
shape_duration_s = (shape_num_samples - 1) * SS_TS; % total duration of grad / rf in seconds
shape_duration_us = shape_duration_s * 1e6; % total duration in µs
shape_timepoints_s = (0:shape_num_samples-1) * SS_TS; % time points of grad / rf values
shape_timepoints_us = shape_timepoints_s * 1e6; % time points of grad / rf values
shape_rf_max_Gauss = max(abs(shape_rf_Gauss)); % maximum rf amplitude in Gauss
shape_rf_max_uT = 100 * shape_rf_max_Gauss; % maximum rf amplitude in µT

% Transform and normalize complex valued shapes
shape_rf_complex_Gauss = shape_rf_Gauss;
shape_rf_complex_uT = 100 * shape_rf_complex_Gauss;
shape_rf_complex_norm = shape_rf_complex_uT / shape_rf_max_uT;
shape_rf_complex_percent = shape_rf_complex_norm * 100;

% Calculate shape and power intergral of the pulse
shape_rf_shape_integral = sqrt( sum(real(shape_rf_complex_norm))^2 + sum(imag(shape_rf_complex_norm))^2) / numel(shape_rf_complex_norm);
shape_rf_power_intergral = shape_rf_complex_norm * shape_rf_complex_norm' / numel(shape_rf_complex_norm);

% Split into magnitude and phase values
shape_rf_mag_Gauss = abs(shape_rf_complex_Gauss);
shape_rf_mag_uT = abs(shape_rf_complex_uT);
shape_rf_mag_percent = abs(shape_rf_complex_percent);
shape_rf_phs_deg = 180 + angle(shape_rf_complex_norm) * 180/pi;
shape_rf_phs_rad = deg2rad(shape_rf_phs_deg);

% Convert gradient shape
shape_grad_mTm = shape_grad_Gausscm * 10;

% Initialize the RFstruct
RFstruct = struct();

% Store filename and description
RFstruct.filename = filename;
RFstruct.desc = desc;
RFstruct.creation_time = creation_time;
clear('filename','desc', 'creation_time')

% Add metabolite structure to RFstruct and remove it from the workspace
RFstruct.mets = mets;
clear('mets');

% Define prefixes and corresponding RFstruct fields
prefixes = {'f_', 'SS_', 'z_', 'shape_'};
fields = {'opts', 'opts', 'opts', 'shapes'};

% Loop over each prefix to collect variables and store them in RFstruct
for idx = 1:length(prefixes)
    prefix = prefixes{idx};
    fieldname = fields{idx};
    vars = who([prefix '*']);
    
    % Loop through each variable and add it to the appropriate struct field
    for i = 1:length(vars)
        varname = vars{i};
        modifiedVarname = varname;

        if startsWith(varname, "shape_")
            modifiedVarname = extractAfter(varname, "shape_");
        end
        RFstruct.(fieldname).(modifiedVarname) = eval(varname);
    end
    clear('i')
    
    % Delete the variables from the workspace
    clear(vars{:});
end

clear('idx', "modifiedVarname", "fields","fieldname","prefix","prefixes","vars","varname", "ans")
