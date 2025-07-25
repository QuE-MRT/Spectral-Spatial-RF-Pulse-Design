function export_RF_PV(RFstruct, path)
%EXPORT_RF_PV Exports RF pulse data to .exc, .meta, and .gp files.
%
%   export_RF_PV(RFstruct, path) exports the RF pulse data contained in the
%   RFstruct to specific file formats within the 'wave' and 'gp'
%   subfolders of the specified 'path'.
%
%   Inputs:
%       RFstruct: A struct created using the "create_RF_struct" script
%
%       path: A string specifying the base directory where the 'wave' and
%             'gp' subfolders will be created.

    % Define subfolder names
    waveFolderName = 'wave';
    gpFolderName = 'gp';

    % Construct full paths for wave and gp directories
    path_waveDir = fullfile(path, waveFolderName);
    path_gpDir = fullfile(path, gpFolderName);

    % Create wave directory if it doesn't exist
    if ~exist(path_waveDir, 'dir')
        mkdir(path_waveDir);
    end

    % Create gp directory if it doesn't exist
    if ~exist(path_gpDir, 'dir')
        mkdir(path_gpDir);
    end

    % define header
    jcampHeader = {
        sprintf('##TITLE= %s',RFstruct.filename);
        '##JCAMP-DX= 5.00 Bruker JCAMP library';
        '##DATA TYPE= Shape Data';
        '##ORIGIN= QuE-MRT';
        '##OWNER= NMRSU';
        sprintf('##DATE= %d/%02d/%02d',RFstruct.creation_time(1:3));
        sprintf('##TIME= %02d:%02d:%02d',RFstruct.creation_time(4:6));
        sprintf('##MINX= %.6e',min(RFstruct.shapes.rf_mag_percent));
        sprintf('##MAXX= %.6e',max(RFstruct.shapes.rf_mag_percent));
        sprintf('##MINY= %.6e',min(RFstruct.shapes.rf_phs_deg));
        sprintf('##MAXY= %.6e',max(RFstruct.shapes.rf_phs_deg));
        '##$SHAPE_EXMODE= Excitation';
        '##$SHAPE_TOTROT= 9.000000e+01';
        '##$SHAPE_BWFAC= 1';
        sprintf('##$SHAPE_INTEGFAC= %.8e',RFstruct.shapes.rf_shape_integral);
        '##$SHAPE_REPHFAC=50';
        '##$SHAPE_TYPE=conventional';
        '##$SHAPE_MODE= 0';
        sprintf('##MAXB1 = %.5e uT',RFstruct.shapes.rf_max_uT);
        sprintf('##NPOINTS= %d',RFstruct.shapes.num_samples);
        sprintf('##DURATION= %.5e ms',RFstruct.shapes.duration_us * 1e-3);
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
    for ii=1:numel(RFstruct.shapes.rf_mag_percent)
        fprintf(exc_File,'%.6e, %.6e\n', ...
            RFstruct.shapes.rf_mag_percent(ii), ...
            RFstruct.shapes.rf_phs_deg(ii) ...
            );
    end
    fprintf(exc_File,'##END= \n');
    fclose(exc_File);
    % open .meta file and write
    meta_File = fopen( ...
        fullfile(path_waveDir,sprintf('%s.meta',RFstruct.filename)),'w');
    fprintf(meta_File,'%.6e %.6e %d', ...
        RFstruct.shapes.duration_us*1e-3, ...
        RFstruct.shapes.rf_shape_integral, ...
        RFstruct.shapes.num_samples);
    fclose(meta_File);
    % open .gp file and write
    gp_File = fopen( ...
        fullfile(path_gpDir,sprintf('%s.gp',RFstruct.filename)),'w');
    for ii=1:numel(RFstruct.shapes.grad_mTm)
        fprintf(gp_File,'%.6e\n',RFstruct.shapes.grad_mTm(ii));
    end
    fclose(gp_File);
    fclose('all');
    % clean workspace
    clear("jcampHeader","exc_File","gp_File","meta_File", ...
        "i","idx","ii","path*");

end