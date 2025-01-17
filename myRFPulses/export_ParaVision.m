%% Write rf and gradient waveform to files
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


current_time = clock;
current_time(6) = round(current_time(6));
filename = sprintf('SpSp_%s_%dT_%s%d_%s%d_%s%d_%s%d_%s%d_%dmm_%d%02d%02d',...
    SS_NUCLEUS,...
    f_B0_Gauss * 1e-4,...
    mets(1).name, mets(1).ang, ...
    mets(2).name, mets(2).ang, ...
    mets(3).name, mets(3).ang, ...
    mets(4).name, mets(4).ang, ...
    mets(5).name, mets(5).ang, ...
    round(10*z_thickness_cm),current_time(1:3));

path_waveDir = fullfile(pwd,'wave');
path_gpDir = fullfile(pwd,'gp');

jcampHeader = {
    sprintf('##TITLE= %s',fullfile(path_waveDir,filename));
    '##JCAMP-DX= 5.00 Bruker JCAMP library';
    '##DATA TYPE= Shape Data';
    '##ORIGIN= NVision';
    '##OWNER= <CAMueller>';
    sprintf('##DATE= %d/%02d/%02d',current_time(1:3));
    sprintf('##TIME= %02d:%02d:%02d',current_time(4:6));
    sprintf('##MINX= %.6e',min(shape_rf_magnitude));
    sprintf('##MAXX= %.6e',max(shape_rf_magnitude));
    sprintf('##MINY= %.6e',min(shape_rf_phase_degree));
    sprintf('##MAXY= %.6e',max(shape_rf_phase_degree));
    '##$SHAPE_EXMODE= Excitation';
    '##$SHAPE_TOTROT= 9.000000e+01';
    '##$SHAPE_BWFAC= 1';
    sprintf('##$SHAPE_INTEGFAC= %.8e',shape_rf_shapeInt);
    '##$SHAPE_REPHFAC=50';
    '##$SHAPE_TYPE=conventional';
    '##$SHAPE_MODE= 0';
    sprintf('##MAXB1 = %.5e uT',shape_rf_max_uT);
    sprintf('##NPOINTS= %d',numel(shape_rf_normalized));
    sprintf('##DURATION= %.5e ms',shape_duration_us * 1e3);
    sprintf('##NUCLEUS= %s',SS_OPTs{1,2});
    sprintf('##FIELD= %.5e T',f_B0_Gauss*1e-4);
    sprintf('##MAXGRAD= %.5e mT/m',SS_OPTs{2,2}*10);
    sprintf('##MAXSLEW= %.5e T/m/s',SS_OPTs{3,2}*10);
    sprintf('##SLICEWIDTH= %.5e cm',z_thickness_cm);
    sprintf('##SLICEOFFSET= %.5e mm',0);
    '##XYPOINTS= (XY..XY)'};
%
if ~exist(path_waveDir,'dir'), mkdir(path_waveDir); end
exc_File = fopen(fullfile(path_waveDir,sprintf('%s.exc',filename)),'w');
for ii=1:numel(jcampHeader)
    fprintf(exc_File,'%s\n',jcampHeader{ii});
end
for ii=1:numel(shape_rf_magnitude)
    fprintf(exc_File,'%.6e, %.6e\n',shape_rf_magnitude(ii),shape_rf_phase_degree(ii));
end
fprintf(exc_File,'##END= \n');
fclose(exc_File);

meta_File = fopen(fullfile(path_waveDir,sprintf('%s.meta',filename)),'w');
fprintf(meta_File,'%.6e %.6e %d', ...
    shape_duration_us*1e3, shape_rf_shapeInt, shape_numel);
fclose(meta_File);

if ~exist(path_gpDir,'dir'), mkdir(path_gpDir); end
gp_File = fopen(fullfile(path_gpDir,sprintf('%s.gp',filename)),'w');
for ii=1:numel(shape_grad_mTm)
    fprintf(gp_File,'%.6e\n',shape_grad_mTm(ii));
end
fclose(gp_File);
fclose('all');

save(fullfile(path_waveDir,sprintf('%s.mat',filename)));
%%
clear("jcampHeader","current_time","exc_File","gp_File","meta_File","i","idx","ii","path*")

