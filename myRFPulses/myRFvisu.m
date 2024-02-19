%% create figure
hfig = figure(123);
clf(hfig)

set(hfig, 'color',  'w')

% 3D plot g, rf , t


subplot(1,3,1)
surf(f,z,abs(mxy), 'EdgeColor','flat')
xlabel('Frequency (Hz)')
ylabel('Position (cm)')
zlabel('Magnitude M_{xy}')
set(gca, 'PlotBoxAspectRatio',[1,1,1], ...
    'FontName', 'MontSerrat',...
    'FontSize',9)

subplot(1,3,2)
surf(f,z,angle(mxy), 'EdgeColor','flat')
xlabel('Frequency (Hz)')
ylabel('Position (cm)')
zlabel('Phase M_{xy}')

set(gca, 'PlotBoxAspectRatio',[1,1,1], ...
    'FontName', 'MontSerrat',...
    'FontSize',9)

colormap('viridis')

%
subplot(1,3,3)

t = (0:1:numel(g)-1) .* opt{4,2} .* 1e3;
for n = 1:numel(g)
    g0(n) = sum(g(1:n));
end


plot3(t, g0, real(rf), ...
    'Color', NVisionColor('lightMint'))
hold on
plot3(t, g0, imag(rf), ...
    'Color', NVisionColor('darkMint'))
plot3(t, g0, abs(rf), ...
    'Color', NVisionColor('charCoal'))

set(gca, 'PlotBoxAspectRatio',[1,1,1], ...
    'FontName', 'MontSerrat',...
    'FontSize',9)

xlabel('Time (ms)')
ylabel('Exc. Gradient 0th Moment')
zlabel('RF Envelope (Gauss)')
pause(0.5)

