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
        RFstruct.(fieldname).(varname) = eval(varname);
    end
    clear('i')
    
    % Delete the variables from the workspace
    clear(vars{:});
end
clear('idx')
%%
clear("fields","fieldname","prefix","prefixes","vars","varname", "ans")
