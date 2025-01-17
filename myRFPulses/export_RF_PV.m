% MATLAB Script: Export Spatial Spectral Pulse Shapes to Bruker ParaVision Format
%
% Description:
% This script processes spatial spectral pulse shapes and converts them into
% a format compatible with Bruker ParaVision systems. The resulting files
% can be used for spectroscopic imaging or related applications in ParaVision.
%
% Features:
% - Supports input from custom pulse shape formats.
% - Generates output files in the required ParaVision structure.
% - Includes error checking for input validation.
%
% Author: Christoph A. Müller
% Date: 17.01.2025
% Version: 1.0
%
% Usage:
% 1. Run pulse generate script and create_RF_struct.m first
% 2. Run the script in MATLAB.
% 3. Exported files will be saved in the specified output directories.

%% Define paths
path_waveDir = fullfile(pwd,'wave');
if ~exist(path_waveDir,'dir'), mkdir(path_waveDir); end

path_gpDir = fullfile(pwd,'gp');
if ~exist(path_gpDir,'dir'), mkdir(path_gpDir); end

% define header
jcampHeader = {
    sprintf('##TITLE= %s',fullfile(path_waveDir,RFstruct.filename));
    '##JCAMP-DX= 5.00 Bruker JCAMP library';
    '##DATA TYPE= Shape Data';
    '##ORIGIN= QuE-MRT';
    '##OWNER= NMRSU'
    sprintf('##DATE= %d/%02d/%02d',RFstruct.creation_time(1:3));
    sprintf('##TIME= %02d:%02d:%02d',RFstruct.creation_time(4:6));
    sprintf('##MINX= %.6e',min(RFstruct.shapes.shape_rf_magnitude));
    sprintf('##MAXX= %.6e',max(RFstruct.shapes.shape_rf_magnitude));
    sprintf('##MINY= %.6e',min(RFstruct.shapes.shape_rf_phase_degree));
    sprintf('##MAXY= %.6e',max(RFstruct.shapes.shape_rf_phase_degree));
    '##$SHAPE_EXMODE= Excitation';
    '##$SHAPE_TOTROT= 9.000000e+01';
    '##$SHAPE_BWFAC= 1';
    sprintf('##$SHAPE_INTEGFAC= %.8e',RFstruct.shapes.shape_rf_shapeInt);
    '##$SHAPE_REPHFAC=50';
    '##$SHAPE_TYPE=conventional';
    '##$SHAPE_MODE= 0';
    sprintf('##MAXB1 = %.5e uT',RFstruct.shapes.shape_rf_max_uT);
    sprintf('##NPOINTS= %d',numel(RFstruct.shapes.shape_rf_normalized));
    sprintf('##DURATION= %.5e ms',RFstruct.shapes.shape_duration_us * 1e3);
    sprintf('##NUCLEUS= %s',RFstruct.opts.SS_OPTs{1,2});
    sprintf('##FIELD= %.5e T',RFstruct.opts.f_B0_Gauss*1e-4);
    sprintf('##MAXGRAD= %.5e mT/m',RFstruct.opts.SS_OPTs{2,2}*10);
    sprintf('##MAXSLEW= %.5e T/m/s',RFstruct.opts.SS_OPTs{3,2}*10);
    sprintf('##SLICEWIDTH= %.5e cm',RFstruct.opts.z_thickness_cm);
    sprintf('##SLICEOFFSET= %.5e mm',0);
    '##XYPOINTS= (XY..XY)'};

% open .exc file and write
exc_File = fopen( ...
    fullfile(path_waveDir,sprintf('%s.exc',RFstruct.filename)),'w');
for ii=1:numel(jcampHeader)
    fprintf(exc_File,'%s\n',jcampHeader{ii});
end
for ii=1:numel(RFstruct.shapes.shape_rf_magnitude)
    fprintf(exc_File,'%.6e, %.6e\n', ...
        RFstruct.shapes.shape_rf_magnitude(ii), ...
        RFstruct.shapes.shape_rf_phase_degree(ii) ...
        );
end
fprintf(exc_File,'##END= \n');
fclose(exc_File);

% open .meta file and write
meta_File = fopen( ...
    fullfile(path_waveDir,sprintf('%s.meta',RFstruct.filename)),'w');
fprintf(meta_File,'%.6e %.6e %d', ...
    RFstruct.shapes.shape_duration_us*1e3, ...
    RFstruct.shapes.shape_rf_shapeInt, ...
    RFstruct.shapes.shape_numel);
fclose(meta_File);

% open .gp file and write
gp_File = fopen( ...
    fullfile(path_gpDir,sprintf('%s.gp',RFstruct.filename)),'w');
for ii=1:numel(RFstruct.shapes.shape_grad_mTm)
    fprintf(gp_File,'%.6e\n',RFstruct.shapes.shape_grad_mTm(ii));
end
fclose(gp_File);
fclose('all');

%% clean workspace
clear("jcampHeader","exc_File","gp_File","meta_File", ...
    "i","idx","ii","path*")

