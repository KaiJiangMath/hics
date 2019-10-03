function objPlot(xRge, yRge, fun, glb_min)

x = linspace(xRge(1), xRge(2), 400);
y = linspace(yRge(1), yRge(2), 400);
for i = 1:1:length(x)
	for j = 1:1:length(y)   
		z(i,j) = fun_value([x(i),y(j)], fun);
	end
end

[X, Y] = meshgrid(x,y);

contour(x,y,z, 50)
%surfc(x,y,z)
surf(x,y,z)
%colorbar
shading interp
alpha(0.4)


xlabel('$x$','Interpreter','latex', 'HorizontalAlignment','right','Rotation',0,'FontName','Times New Roman', 'FontSize',35);
ylabel('$y$','Interpreter','latex', 'HorizontalAlignment','right','Rotation',0,'FontName','Times New Roman', 'FontSize',35);

%colorbar;

plot(glb_min(1), glb_min(2), 'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0], ...
	'MarkerSize', 8, 'Marker','pentagram','Color',[1 0 0]);

