% Demonstration of multiband hyperpolarized C-13 pulse designs for MRSI/CSI, including explanation
% of some key pulse characteristics

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

% multiband pulse
fprintf(1, '************************************************************\n')
fprintf(1, 'Here''s a C13 multiband excitation pulse example, for [1-13C]pyr+13C-urea\n');
fprintf(1, 'dynamic MR spectroscopic or chemical shift imaging on a 3T clinical system\n');
fprintf(1, '************************************************************\n')

ss_opt([]);

% GENERAL PULSE PARAMETERS
ss_type = 'EP Whole';  % Echo-planar design
ptype = 'ex';  % excitation pulse
opt = ss_opt({'Nucleus', 'Carbon', ...
	      'Max Duration', 25e-3});

% force pulse design to optimize for center of frequency specification
fctr = 0;  

% SPECTRAL PULSE PARAMETERS 
B0 = 11.7e4; % G
df = 0.5e-6 * B0 * SS_GAMMA; % 0.5 ppm = gamma_C13 * B0 * 0.5e-6
% metabolite			frequency (Hz)		freq bandwidth (Hz)		flip angle (deg)	allowed ripple
% mets(1).name = 'bic'; 	mets(1).f = -1160; 	mets(1).df = 1.5*df; 	    mets(1).ang = 10; 	mets(1).d = .05;
% mets(2).name = 'urea'; 	mets(2).f = -935; 	mets(2).df = 1.5*df; 		mets(2).ang = 0; 	mets(2).d = .05;
% mets(3).name = 'pyr'; 	mets(3).f = 0;  	mets(3).df = 1.5*df; 	    mets(3).ang = 0; 	mets(3).d = .05;
% mets(4).name = 'ala'; 	mets(4).f = 720; 	mets(4).df = 1.5*df;	    mets(4).ang = 0; 	mets(4).d = .05;
% mets(5).name = 'lac'; 	mets(5).f = 1535; 	mets(5).df = 1.5*df;        mets(5).ang = 0; 	mets(5).d = .05;
mets(1).name = 'ure'; 	mets(1).f = -945; 	mets(1).df = 2*df; 	    mets(1).ang = 0; 	mets(1).d = .05;
mets(2).name = 'pyr'; 	mets(2).f = 0; 	    mets(2).df = 2*df; 		mets(2).ang = 10; 	mets(2).d = .05;
mets(3).name = 'ala'; 	mets(3).f = 710;  	mets(3).df = 2*df; 	    mets(3).ang = 0; 	mets(3).d = .05;
mets(4).name = 'la1c'; 	mets(4).f = 1515; 	mets(4).df = 2*df; 		mets(4).ang = 60; 	mets(4).d = .05;
mets(5).name = 'bic';   mets(5).f = -1155;  mets(5).df = 2*df;      mets(5).ang = 0;    mets(5).d = .05;


% create vectors of angles, ripples, and band edges for input to pulse design
[fspec, a_angs, d] = create_freq_specs(mets, fctr);
s_ftype = 'min';  % minimimum-phase spectral filter

% SPATIAL PULSE PARAMETERS
z_ftype = 'ls';  % least-squares filter design
z_d1 = 0.05;  z_d2 = 0.05;  % slice profile pass and stop-band ripples, respectively

% multiband pulse - thicker slice
z_thk = 1.5;  % thickness (cm)
z_tb = 2; % time-bandwidth, proportional to profile sharpness

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);
set(gcf,'Name', '[1-13C]pyr+13C-urea Multiband - slab');

% for saving pulses:
% ss_save(g,rf,max(a_angs),z_thk, [], 'GE', fspec, a_angs, root_fname);


%% Write rf and gradient waveform to files
max_rfuT = 100*max(abs(rf));
bpB1uT   = 1 / (SS_GAMMA * 4e-5);       % B1 for 1 ms, 90 deg bp [uT]
powerFactor = (max_rfuT / bpB1uT)^2;    % RF power is propto B1^2
rfn = 100 * rf / max_rfuT;
sInt = sqrt( sum(real(rfn))^2 + sum(imag(rfn))^2) / numel(rfn);
pInt = rfn * rfn' / numel(rfn);
length_ms = 1e3 * SS_TS * numel(rfn);

rfn = rfn * 100;
mag = abs(rfn);
ph_deg = 180 + angle(rfn) * 180/pi;
grad_mTm = g * 10;

slOffs  = 0; % mm
nSlices = numel(slOffs);

curTime = clock;
curTime(6) = round(curTime(6));
filName = sprintf('SpSp_13C_11p7T_%s%d_%s%d_%s%d_%s%d_%dmm_%d%02d%02d',...
    mets(1).name, mets(1).ang, ...
    mets(2).name, mets(2).ang, ...
    mets(3).name, mets(3).ang, ...
    mets(4).name, mets(4).ang, ...
    round(10*z_thk),curTime(1:3));

waveDir = fullfile(pwd,'wave');
gpDir = fullfile(pwd,'gp');

jcampHeader = {
    sprintf('##TITLE= %s',fullfile(waveDir,filName));
    '##JCAMP-DX= 5.00 Bruker JCAMP library';
    '##DATA TYPE= Shape Data';
    '##ORIGIN= NVision';
    '##OWNER= <CAMueller>';
    sprintf('##DATE= %d/%02d/%02d',curTime(1:3));
    sprintf('##TIME= %02d:%02d:%02d',curTime(4:6));
    sprintf('##MINX= %.6e',min(mag));
    sprintf('##MAXX= %.6e',max(mag));
    sprintf('##MINY= %.6e',min(ph_deg));
    sprintf('##MAXY= %.6e',max(ph_deg));
    '##$SHAPE_EXMODE= Excitation';
    '##$SHAPE_TOTROT= 9.000000e+01';
    '##$SHAPE_BWFAC= 1';
    sprintf('##$SHAPE_INTEGFAC= %.8e',sInt);
    '##$SHAPE_REPHFAC=50';
    '##$SHAPE_TYPE=conventional';
    '##$SHAPE_MODE= 0';
    sprintf('##MAXB1 = %.5e uT',max_rfuT);
    sprintf('##NPOINTS= %d',numel(rf));
    sprintf('##DURATION= %.5e ms',length_ms);
    sprintf('##NUCLEUS= %s',opt{1,2});
    sprintf('##FIELD= %.5e T',B0*1e-4);
    sprintf('##MAXGRAD= %.5e mT/m',opt{2,2}*10);
    sprintf('##MAXSLEW= %.5e T/m/s',opt{3,2}*10);
    sprintf('##SLICEWIDTH= %.5e cm',z_thk);
    sprintf('##SLICEOFFSET= %.5e mm',0);
    '##XYPOINTS= (XY..XY)'};
%
if ~exist(waveDir,'dir'), mkdir(waveDir); end
rfFile = fopen(fullfile(waveDir,sprintf('%s.exc',filName)),'w');
for ii=1:numel(jcampHeader)
    fprintf(rfFile,'%s\n',jcampHeader{ii});
end
for ii=1:numel(rf)
    fprintf(rfFile,'%.6e, %.6e\n',mag(ii),ph_deg(ii));
end
fprintf(rfFile,'##END= \n');
fclose(rfFile);

metaFile = fopen(fullfile(waveDir,sprintf('%s.meta',filName)),'w');
fprintf(metaFile,'%.6e %.6e %d', ...
    length_ms, sInt, numel(rf));
fclose(metaFile);

if ~exist(gpDir,'dir'), mkdir(gpDir); end
gradFile = fopen(fullfile(gpDir,sprintf('%s.gp',filName)),'w');
for ii=1:numel(g)
    fprintf(gradFile,'%.6e\n',grad_mTm(ii));
end
fclose(gradFile);

save(fullfile(waveDir,sprintf('%s.mat',filName)));
return


