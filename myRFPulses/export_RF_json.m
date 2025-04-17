function export_RF_json(RFstruct, path)
%EXPORT_RF_JSON Exports RF pulse data to a JSON file.
%
%   export_RF_json(RFstruct, path) exports the RF pulse data contained in the
%   RFstruct to a JSON file in the specified 'path'.
%
%   Inputs:
%       RFstruct: A struct created using the "create_RF_struct" script
%
%       path: A string specifying the directory where the JSON file will be created.

    % Create JSON file name
    filename_base = RFstruct.filename;
    filename_json = [filename_base, '.json'];

    % create JSONstruct
    JSONstruct.filename = filename_base;
    JSONstruct.desc = RFstruct.desc;
    JSONstruct.rf_asym = 0.5;
    JSONstruct.rf_dur_fix.dur = RFstruct.shapes.duration_us; % [us]
    JSONstruct.rf_dur_fix.fix = true;
    JSONstruct.rf_abs = RFstruct.shapes.rf_mag_percent / 100 * RFstruct.shapes.rf_max_uT; % [uT]
    JSONstruct.rf_phs = RFstruct.shapes.rf_phs_rad; % [rad]
    JSONstruct.gradz_v = RFstruct.shapes.grad_mTm; % [mT/m]
    JSONstruct.gradz_t = RFstruct.shapes.timepoints_us; % [us]

    % write string from JSONstruct
    string = jsonencode(JSONstruct,"ConvertInfAndNaN",true, 'PrettyPrint', true);

    % Construct the full path for the JSON file
    fullFilePath = fullfile(path, filename_json);

    % write file
    file_id = fopen(fullFilePath, 'w');
    fprintf(file_id, '%s',string);
    fclose(file_id);

end