
close all;

a=load('nfevl_ackley10d.txt');
dim = 10;

figure(1) 
hold on 

nl = length(a(:,1));
ind = 0;

for i=1:1:nl
	n = a(i,1);
	x = ind:1:ind+n-1;
	ind = ind+n;

	fv = zeros(1,n);
	rv = zeros(1,n);
	fv(:) = a(i,1+dim+1);
	rv(:) = a(i,1+dim+2);

	semilogy(x, fv, 'k-')
end

%
%
%
%
%semilogy(a(:,1+dim), 'k-')
%%plot(a(:,4), 'b-')
%
%xlabel('Number of Evaluations')
%ylabel('Objective Value')
%
%%ylim([10^(-5.6) 10^(1.5)])
%%xlim([0, 21000])

set(gca, 'FontName','Times New Roman', 'FontSize',30);
set(gcf, 'color', 'white', 'unit', 'normalized', 'position', [0.1, 0.05, 0.8, 0.85])
box on; 
axis on; 
axis square;
set(0,'defaultaxeslinewidth',2)
set(0,'defaultlinelinewidth',3)

%%%%%  插入子图
%axes11 = axes('Position', [0.51 0.6 0.24 0.24],'FontName','Times New Roman', 'FontSize', 25);
%semilogy(a(:,6), 'b-')
%axis square
%%ylim([10^(-5.5) 10^(0.1)])
%%xlim([0, 58000])
%
%ylabel('r','HorizontalAlignment','right','Rotation',0)
%set(gca, 'FontName','Times New Roman', 'FontSize', 30);
