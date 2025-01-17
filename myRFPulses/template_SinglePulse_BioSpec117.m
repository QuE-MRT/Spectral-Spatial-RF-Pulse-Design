% Reset SS package globals
%
clear all
clc
ss_opt([]);
ss_globals;


% single band pulse -only pulse
fprintf(1, '\nHere''s a C13 single-band excitation pulse to excite only one metabolite\n');
fprintf(1, 'for a 11.7T preclinical system\n\n');

% GENERAL PULSE PARAMETERS
ss_type = 'EP Whole';
ptype = 'ex';  % excitation pulse
opt = ss_opt({'Nucleus', 'Carbon', ...
	      'Max Duration', 25e-3, ...
          'Spect Correct', 1, ...
          'Sample Time', 8e-6});
% ,...
%     'Max Grad', 76.185 ,...
%     'Max Slew', 634.875});
      
% SPECTRAL PULSE PARAMETERS  - large pass/stop bands chosen for wide
% supression regions
B0 = 3e4; % G
df = 1e-6 * B0 * SS_GAMMA; % 1 ppm = gamma_C13 * B0 * 0.5e-6
% metabolite			frequency (Hz)		freq bandwidth (Hz)		flip angle (deg)    allowed ripple
mets(1).name = 'passband'; 	mets(1).f = 0;      mets(1).df = 1*df; 		mets(1).ang = 90;   mets(1).d = 0.01;
mets(2).name = 'stopband'; 	mets(2).f = 250; 	mets(2).df = 175;       mets(2).ang = 0;    mets(1).d = 0.05;
mets(3).name = 'stopband'; 	mets(3).f = -250; 	mets(3).df = 175; 		mets(3).ang = 0;    mets(1).d = 0.05;


% create vectors of angles, ripples, and band edges for input to pulse design
[fspec, a_angs, d] = create_freq_specs(mets);
fctr = 1;  % force pulse design to optimize for center of frequency specification
s_ftype = 'min';  % minimum-phase spectral filter

% SPATIAL PULSE PARAMETERS

z_thk = 1;  % thickness (cm)
z_tb = 4; % time-bandwidth, proportional to profile sharpness
z_ftype='ls';  % least-squares filter design
z_d1 = 0.01;  z_d2 = 0.01;  % slice profile pass and stop-band ripples, respectively

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);
set(gcf,'Name', '[1-13C]lac only for 3T clinical system');
