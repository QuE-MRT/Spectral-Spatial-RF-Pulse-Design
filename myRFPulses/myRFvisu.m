% MATLAB Script: Visualize Spatial Spectral (SPSP) RF Pulse
%
% Description:
% Alternative visualization to ss_plot()
%
% Author: Christoph A. Müller
% Date: 17.01.2025
% Version: 1.0
%
% Usage:
% 1. Load RFstruct
% 2. Run the script to generate plots and visualize the pulse.


%% create figure
hfig = figure(123);
clf(hfig)
set(hfig, 'color',  'w')

% plot pulse amplitude over gradient 0th moment, i.e. over k-space
% coordinate
subplot(1,3,1)
gradient0thmoment = cumsum(RFstruct.shapes.shape_grad_mTm);

plot3(RFstruct.shapes.shape_timesamples_us, ...
    gradient0thmoment, ...
    real(RFstruct.shapes.shape_rf_uT))
hold on
plot3(RFstruct.shapes.shape_timesamples_us, ...
    gradient0thmoment, ...
    imag(RFstruct.shapes.shape_rf_uT))

set(gca, 'PlotBoxAspectRatio',[1,1,1], ...
    'FontName', 'MontSerrat',...
    'FontSize',9)
title('3D representation of Pulse envelope on kspace-time trajectory')
xlabel('Time (s)')
ylabel('Exc. Gradient 0th Moment')
zlabel('RF Envelope (Gauss)')

subplot(1,3,2)
plot(RFstruct.shapes.shape_timesamples_us, ...
    real(RFstruct.shapes.shape_rf_uT))
title('Temporal RF envelope representation')
xlabel('Time (s)')
ylabel('RF Envelope (Gauss)')
set(gca, 'PlotBoxAspectRatio',[1,1,1], ...
    'FontName', 'MontSerrat',...
    'FontSize',9)

subplot(1,3,3)
plot(gradient0thmoment, ...
    real(RFstruct.shapes.shape_rf_uT))
title('Spatial kSpace RF envelope representation')
xlabel('Exc. Gradient 0th Moment')
ylabel('RF Envelope (Gauss)')
set(gca, 'PlotBoxAspectRatio',[1,1,1], ...
    'FontName', 'MontSerrat',...
    'FontSize',9)


