function cdisPlot( p0, p1, rad, val, obsPoints, xRge, yRge, ... 
				   evalu, plotNum, mflag, mystr, fun, glb_min)


%%% mflag: show the better point or observed points
%%% mflag = 1, show reference point and the better point
%%% mflag = 0, show reference point and sample points

%%% rad(1) search radius
%%% rad(2) initail search radius
%%% val: functional value f(p0) 
%%% evalu(1): # of evaluations of each iterations
%%% evalu(2): total # of evaluations 
%%% [xRge, yRge]: the range of objective function
%%% plotNum: the index of plot figures
%%% obsPoints: the sampling set
%%% mystr: information of iterative process

if strcmp(mystr, 'Refinement')
	mystr = 'Refinement';
elseif strcmp(mystr, 'sample')
	mystr = 'Sampling';
elseif strcmp(mystr, 'Optimal')
	mystr = 'Find Better State'; 
elseif strcmp(mystr, 'Initialization')
	mystr = 'Initial State'; 
end

r=rad(1); r0=rad(2);

radstr = '';
%if r==r0
%    radstr = '';
%elseif r<r0
%    radstr = 'Contract Radius r'; 
%elseif r>r0
%    radstr = 'Expand Radius r'; 
%end

fn = sprintf('%d', plotNum);
radius = sprintf('Radius: %.3f', r);

x0 = p0(1); y0 = p0(2);
x1 = p1(1); y1 = p1(2);

pit0 = sprintf('(x_0,y_0)=(%.5f, %.5f)', x0, y0);
fval = sprintf('f = %e', val);
radius = sprintf('Radius: %.4e', r);
evals = sprintf('Evaluations on each time: %d', evalu(1));
tolevals = sprintf('Total Evaluations: %d', evalu(2));


xleft = xRge(1) - 0.1;
xrigt = xRge(2) + 0.1;
yleft = yRge(1) - 0.1;
yrigt = yRge(2) + 0.1;

figure(1)
hold on;
%% box on; 


set(gcf, 'Color', 'white')
set(gca, 'FontName','Times New Roman', 'FontSize', 20);
set(findall(gcf, 'type', 'line'), 'LineWidth', 2)
set(gcf, 'unit', 'normalized', 'position', [0.05, 0.15, 0.9, 0.7]);
% grid on

axis equal
axis on
%axis square

set(gca, 'xlim', [xleft, xrigt])
set(gca, 'ylim', [yleft, yrigt])
%view(43, 44)
%view(43, -26)
view(0, 90)


objPlot(xRge, yRge, fun, glb_min);

if mflag == 0
	disp('hello')
	obsPitPlot(x0, y0, obsPoints, r, val)
else 
	movePitPlot(x0, y0, x1, y1, r, val)
end

% Create textbox
annotation('textbox',...
	[0.1 0.7 0.8 0.346094946401225],...
	'String',{mystr, radstr, radius, pit0, fval, evals, tolevals},...
	'EdgeColor', 'none', ...
	'FontName','Times New Roman', 'FontSize', 20);

%axis normal;
%imshow(strain_image, 'border', 'tight', 'initialmagnification', 'fit');
%set(gcf, 'position', [0, 0, 500, 500]);
%set(gcf, 'unit', 'normalized', 'position', [0.05, 0.05, 0.85, 0.7]);

% set(gcf, 'unit', 'normalized', 'position', [0.05, 0.15, 0.95, 0.7]);
% set(gca, 'LooseInset', get(gca, 'TightInset'))
%set(gcf, 'position', [0, 0, 1300, 500]);



fn = sprintf('_%d', plotNum);
%set(gcf,'PaperPositionMode','auto'); % make sure that the saved pic are the same as the display fig
%print(gcf, '-dpng', '-r280', [strcat(fun, fn), '.png']);

set(gca,'LineWidth',2);
set(gca,'LooseInset',[0,0,0,0]);
set(gcf,'PaperUnits','inches','PaperPosition',[0,0,10,10]);

if exist('ackley/') == 0
    mkdir('ackley/');
end
saveas(gcf,['ackley/ackley_CHC', fn,'.png'],'png');

