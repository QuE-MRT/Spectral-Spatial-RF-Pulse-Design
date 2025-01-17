% MATLAB Script: Export Spatial Spectral Pulse Shapes to JSON Dictionary
%
% Description:
% This script processes spatial spectral pulse shapes and exports them
% into a JSON dictionary file. 
%
% Features:
% - Converts pulse shape data into a structured JSON format.
% - Customizable for specific dictionary key-value mappings.
%
% Author: Christoph A. Müller
% Date: 17.01.2025
% Version: 1.0
%
% Usage:
% 1. Run pulse generate script and create_RF_struct.m first
% 2. Run the script in MATLAB.
% 3. The JSON file will be saved in the specified output directory.

% Create JSON file name
filename_json = [RFstruct.filename, '.json']

% create JSONstruct
JSONstruct.filename = RFstruct.filename;
JSONstruct.desc = RFstruct.desc;
JSONstruct.rf_asym = 0.5;
JSONstruct.rf_dur_fix.dur = RFstruct.shapes.shape_duration_us; % [us]
JSONstruct.rf_dur_fix.fix = true; 
JSONstruct.rf_abs = RFstruct.shapes.shape_rf_magnitude; % [uT]
JSONstruct.rf_phs = RFstruct.shapes.shape_rf_phase_degree; % [rad]
JSONstruct.gradz_v = RFstruct.shapes.shape_grad_mTm; % [mT/m]
JSONstruct.gradz_t = RFstruct.shapes.shape_timesamples_us; % [us]
    
% write string from JSONstruct
string = jsonencode(JSONstruct,"ConvertInfAndNaN",true);

% Replace desc to _desc and filename to _filename 
% (matlab can not create <sth>._name field names)
string = strrep(string, 'desc', '_desc');
string = strrep(string, 'filename', '_filename');

% write file
file_id = fopen(filename_json, 'w');
fprintf(file_id, '%s',string);
fclose(file_id);

% finished to be loaded in gammaSTAR

