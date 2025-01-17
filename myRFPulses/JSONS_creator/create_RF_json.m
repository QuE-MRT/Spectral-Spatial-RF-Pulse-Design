%% convert RF-struct into JSON dictionary 
% for readin in gammaSTAR

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

