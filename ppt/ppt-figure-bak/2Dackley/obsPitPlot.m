function obsPitPlot(x0, y0, obsPoints, r, val)

cicX = x0-r:0.005:x0+r;
p = length(cicX);
cicY = zeros(1,p);
cicY2 = zeros(1,p);
vls = zeros(1, p);
for i = 1:p
	cicY(i)  =  sqrt(r^2-(cicX(i)-x0)^2)+y0;
	cicY2(i) = -sqrt(r^2-(cicX(i)-x0)^2)+y0;
	vls(i) = val;
end

plot3(cicX, real(cicY), vls, 'k-', 'LineWidth', 2)
plot3(cicX, real(cicY2), vls, 'k-', 'LineWidth',2)


plot3(x0, y0, val, 'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', [0,0,0], ... 
		'MarkerSize', 8, 'Marker', 'o')

%plot(x0, y0, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1], ... 
%        'MarkerSize', 8, 'Marker', 'o')

for i = 1:1:length(obsPoints(:,1))
	plot3(obsPoints(i,1), obsPoints(i,2), val, ...
	'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', [0,0,0], ... 
		'MarkerSize', 7, 'Marker', 'o')
%    plot(obsPoints(i,2), obsPoints(i,1), ...
%    'MarkerFaceColor', [1,1,1], 'MarkerEdgeColor', [0,0,0], ... 
%        'MarkerSize', 6, 'Marker', 'o')
end
