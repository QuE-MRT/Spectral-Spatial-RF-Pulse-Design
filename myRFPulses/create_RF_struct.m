% Initialize the RFstruct
RFstruct = struct();

% Store filename and description
RFstruct.filename = filename;
RFstruct.desc = desc;
clear('filename','desc')

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
    
    % Delete the variables from the workspace
    clear(vars{:});
end
%%
clear("fields","fieldname","prefix","prefixes","vars","varname")
