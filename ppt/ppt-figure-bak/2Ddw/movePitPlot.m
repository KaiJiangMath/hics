function movePitPlot(x0, y0, x1, y1, r, val)

cicX = x0-r:0.005:x0+r;
cicY = [];
p = length(cicX);
vls = zeros(p, 1);
for i = 1:p
	cicY(i)  =  sqrt(r^2-(cicX(i)-x0)^2)+y0;
	cicY2(i) = -sqrt(r^2-(cicX(i)-x0)^2)+y0;
	vls(i) = val;
end

plot3(cicX, real(cicY), vls, 'k-', 'LineWidth',  1.5)
plot3(cicX, real(cicY2), vls, 'k-', 'LineWidth', 1.5)

%plot(cicX, cicY,  'k-', 'LineWidth', 1.5)
%plot(cicX, cicY2, 'k-', 'LineWidth', 1.5)


plot3(x0, y0, val, 'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', 'k', ... 
		'MarkerSize', 5, 'Marker', 'o')
plot3(x1, y1, val, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ... 
		'MarkerSize', 6, 'Marker', 'o')

%plot(x0, y0, 'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', [1,1,1], ... 
%        'MarkerSize', 6, 'Marker', 'o')
%plot(x1, y1, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1], ... 
%        'MarkerSize', 8, 'Marker', 'o')

