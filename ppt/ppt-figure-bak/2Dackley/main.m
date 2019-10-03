close all;
clear all;
clc
format long

%profile on;
set(0,'DefaultFigureVisible','off');

%addpath(genpath(pwd));

funflag = 'ackley';
%funflag = 'gauss';
%funflag = 'dw';
%funflag = 'bukin6';
%funflag = 'drop';
%funflag = 'egg';
%funflag = 'zakharov';
%funflag = 'camelfun';
%funflag = 'woods';
%funflag = 'arwhead';
%funflag = 'chrosen';
%funflag = 'BDQRTIC';
%funflag = 'BRYBND';
%funflag = 'powell';
%funflag = 'spheref';


if strcmp(funflag, 'ackley')
	DIM = 2; 
	L = 5;
	bdrge = zeros(DIM,2);
	for i=1:1:DIM
		bdrge(i,:) = [ -1.0*L, 1.0*L ];
	end
	glb_min = zeros(1, DIM);
	m = 2;   %%% the maximum number of refinement
	r0 = 1;           %%%% inintal radius %%%%
	%%%% inintal values %%%%
	%p0 = [2.12,5.426,-3.45,3.12,-3.0,0.99,1.57,-4.11,5.1,-1.03];
	p0 = [2.12, 3.426];
	%%%% random inintal values %%%%
%    p0 = rand(1, DIM);
%    for i = 1:1:DIM
%    %   p0(i) = p0(i)*( bdrge(i,2)-bdrge(i,1)) + bdrge(i,1); 
%       p0(i) = p0(i)*( (bdrge(i,2)-r0)- (bdrge(i,1)+r0) ) + (bdrge(i,1)+r0); 
%    end
end

if strcmp(funflag, 'woods')
	DIM = 200; 
	L = 40;
	bdrge = zeros(DIM,2);
	for i=1:1:DIM
		bdrge(i,:) = [ -1.0*L, 1.0*L ];
	end
	glb_min = ones(1, DIM);
	m = 7;   %%% the maximum number of refinement
	r0 = 6;           %%%% inintal radius %%%%
	%%%% inintal values %%%%
	%%%% random inintal values %%%%
%    p0 = rand(1, DIM);
%    for i = 1:1:DIM
%    %   p0(i) = p0(i)*( bdrge(i,2)-bdrge(i,1)) + bdrge(i,1); 
%       p0(i) = p0(i)*( (bdrge(i,2)-r0)- (bdrge(i,1)+r0) ) + (bdrge(i,1)+r0); 
%    end
%    pause()
	p0 = zeros(1, DIM);
	for i = 1:1:DIM
		if mod(i,2) == 1
			p0(i) = -3.0;
		end
		if mod(i,2) == 0
			p0(i) = -1.0;
		end
	end
%    p0 = [-3, -1, -3, -1];  %%% hard init
%    p0 = [1.1, 1.2, 1.3, 1.4]; %%% easy init
end

if strcmp(funflag, 'BDQRTIC')
	DIM = 4;
	L = 10;
	bdrge = zeros(DIM,2);
	for i=1:1:DIM
		bdrge(i,:) = [ -1.0*L, 1.0*L ];
	end
	glb_min = zeros(1, DIM); 
	m = 5;   %%% the maximum number of refinement
	r0 = 3;           %%%% inintal radius %%%%
	%%%% inintal values %%%%
	p0 = ones(1, DIM);
	%%%% random inintal values %%%%
	p0 = rand(1, DIM);
%    pause()
	for i = 1:1:DIM
	%   p0(i) = p0(i)*( bdrge(i,2)-bdrge(i,1)) + bdrge(i,1); 
	   p0(i) = p0(i)*( (bdrge(i,2)-r0)- (bdrge(i,1)+r0) ) + (bdrge(i,1)+r0); 
	end
end

if strcmp(funflag, 'powell')
	DIM = 100; 
	L = 10;
	bdrge = zeros(DIM,2);
	for i=1:1:DIM
		bdrge(i,:) = [ -1.0*4, 1.0*5 ];
	end
	glb_min = zeros(1, DIM);
	m = 5;   %%% the maximum number of refinement
	r0 = 3;           %%%% inintal radius %%%%
	%%%% inintal values %%%%
%    p0 = j[-3, -1, -3, -1];  %%% hard init
	%%%% random inintal values %%%%
	p0 = rand(1, DIM);
%    pause()
	for i = 1:1:DIM
	%   p0(i) = p0(i)*( bdrge(i,2)-bdrge(i,1)) + bdrge(i,1); 
	   p0(i) = p0(i)*( (bdrge(i,2)-r0)- (bdrge(i,1)+r0) ) + (bdrge(i,1)+r0); 
	end
end

if strcmp(funflag, 'spheref')
	DIM = 20; 
	L = 5.12;
	bdrge = zeros(DIM,2);
	for i=1:1:DIM
		bdrge(i,:) = [ -1.0*L, 1.0*L ];
	end
	glb_min = zeros(1, DIM);
	m = 5;   %%% the maximum number of refinement
	r0 = 3;           %%%% inintal radius %%%%
	%%%% inintal values %%%%
%    p0 = [-3, -1, -3, -1];  %%% hard init
	%%%% random inintal values %%%%
	p0 = rand(1, DIM);
	for i = 1:1:DIM
	%   p0(i) = p0(i)*( bdrge(i,2)-bdrge(i,1)) + bdrge(i,1); 
	   p0(i) = p0(i)*( (bdrge(i,2)-r0)- (bdrge(i,1)+r0) ) + (bdrge(i,1)+r0); 
	end
%    pause()
end

if strcmp(funflag, 'arwhead')
	DIM = 1000; 
	L = 10;
	bdrge = zeros(DIM,2);
	for i=1:1:DIM
		bdrge(i,:) = [ -1.0*L, 1.0*L ];
	end
	glb_min = ones(1, DIM);
	glb_min(end) = 0;
	m = 5;   %%% the maximum number of refinement
	r0 = 3;           %%%% inintal radius %%%%
	%%%% inintal values %%%%
	p0 = ones(1, DIM);
	%%%% random inintal values %%%%
%    p0 = rand(1, DIM)
%    pause()
%    for i = 1:1:DIM
%    %   p0(i) = p0(i)*( bdrge(i,2)-bdrge(i,1)) + bdrge(i,1); 
%       p0(i) = p0(i)*( (bdrge(i,2)-r0)- (bdrge(i,1)+r0) ) + (bdrge(i,1)+r0); 
%    end
end

if strcmp(funflag, 'chrosen')
	DIM = 120; 
	L = 50;
	bdrge = zeros(DIM,2);
	for i=1:1:DIM
		bdrge(i,:) = [ -1.0*L, 1.0*L ];
	end
	glb_min = ones(1, DIM);
	m = 5;   %%% the maximum number of refinement
	r0 = 5;           %%%% inintal radius %%%%
	%%%% inintal values %%%%
	p0 = -1*ones(1, DIM);
	%%%% random inintal values %%%%
%    p0 = rand(1, DIM);
%    pause()
%    for i = 1:1:DIM
%    %   p0(i) = p0(i)*( bdrge(i,2)-bdrge(i,1)) + bdrge(i,1); 
%       p0(i) = p0(i)*( (bdrge(i,2)-r0)- (bdrge(i,1)+r0) ) + (bdrge(i,1)+r0); 
%    end
end

if strcmp(funflag, 'zakharov')
	DIM = 10; 
	L = 10;
	bdrge = zeros(DIM,2);
	for i=1:1:DIM
		bdrge(i,:) = [ -1.0*L, 1.0*L ];
	end
	glb_min = zeros(1, DIM);
	m = 5;   %%% the maximum number of refinement
	r0 = 3;           %%%% inintal radius %%%%
	%%%% inintal values %%%%
	%p0 = [2.12,5.426,-3.45,3.12,-3.0,0.99,1.57,-4.11,5.1,-1.03];
	%p0 = [6.12,5.426];
	%%%% random inintal values %%%%
	p0 = rand(1, DIM);
	for i = 1:1:DIM
	%   p0(i) = p0(i)*( bdrge(i,2)-bdrge(i,1)) + bdrge(i,1); 
	   p0(i) = p0(i)*( (bdrge(i,2)-r0)- (bdrge(i,1)+r0) ) + (bdrge(i,1)+r0); 
	end
end

if strcmp(funflag, 'camelfun')
	DIM = 2;

	bdrge = [-3, 3;
			 -2, 2];
	glb_min = [-0.0898, -0.7126];
	r0 = 0.7;
	m = 5;
	p0 = [-2, 0.8];
end

if strcmp(funflag, 'bukin6')
	DIM = 2; 
	bdrge = [-15, -0
			  -5,  5];
	glb_min = [-10, 1];
	r0 = 3;           %%%% inintal radius %%%%
end
if strcmp(funflag, 'drop')
	DIM = 2; 
	bdrge = [-5.12, 5.12
			 -5.12, 5.12];
	glb_min = [0, 0];
	r0 = 3;           %%%% inintal radius %%%%
end
if strcmp(funflag, 'egg')
	DIM = 2; 
	bdrge = [-512, 600
			 -512, 600];
	glb_min = [512, 404.2319];
	r0 = 100;           %%%% inintal radius %%%%
end
if strcmp(funflag, 'gauss')
	DIM = 2;
	L = 1;
	bdrge = zeros(DIM,2);
	for i=1:1:DIM
		bdrge(i,:) = [ -1.0*L, 1.0*L ];
	end
	glb_min = zeros(1,DIM);
	r0 = 0.3;           %%%% inintal radius %%%%
	m = 2;
	p0 = [0.57, -0.6];
%    p0 = [6.7, -8, 3];
%    format short
%    p0 = rand(1, DIM)
%    p0 = [ 0.406726915089364   0.666931533207499 0.933725659545930   0.810950032238265   0.484548271834990 0.756749210065512   0.417047453742767   0.971785992989294 0.987974701230832   0.864147529031203 0.388883775912429 0.454741828039112 0.246687197638079   0.784423093024177 0.882837605583052]
%    size(p0)
%    pause()

%    for i = 1:1:DIM
%    %   p0(i) = p0(i)*( bdrge(i,2)-bdrge(i,1)) + bdrge(i,1); 
%       p0(i) = p0(i)*( (bdrge(i,2)-r0)- (bdrge(i,1)+r0) ) + (bdrge(i,1)+r0); 
%    end
end
if strcmp(funflag, 'dw')
	DIM = 2; 
	L = 5;
	bdrge = zeros(DIM,2);
	for i=1:1:DIM
		bdrge(i,:) = [ -1.0*L, 1.0*L ];
	end
	glb_min = [0, 0];
	r0 = 0.5;           %%%% inintal radius %%%%
	m = 2;   %%% the maximum number of refinement
	p0 = [3.2, 1.5];
end
% bdrge

funflag
DIM
% pause()

%%%%% search region %%%%%

epsVal  = 1.0e-4; %%%% convergent criterion of radius
epsGrad = 1.0e-5;
eta = (sqrt(5)-1)/2;    %%% shrinkage factor > 0
gamma = 0.1;  %%% must be > 0

fname = sprintf('./simplice/%dD_%d_simplice.mat', DIM, m);

tic
try 
	%%%% load from file data
	load(fname);
catch 
	simplex = regular_simplex(DIM);
	rotation_simplice = generate_rotation_simplex( DIM, m, simplex );
end
gen_simplices_time = toc

%%%%%%%%  some flags
if DIM == 2 
	plotflag = 1;
else
	plotflag = 0;
end
flag = [plotflag];

parameter = [r0, epsVal, epsGrad, eta, gamma, m]
fname = sprintf('log_SHC_%s_%dD.txt', funflag, DIM);
fid = fopen(fname, 'w');

tic

[r, p0, pdist] = ASHCgauss(funflag, bdrge, parameter, p0, rotation_simplice, glb_min, fid, flag);
%[r, pdist]
%larg = sqrt(DIM);
%k = 1;
%% pause()
%while pdist > r
%    parameter = [pdist, epsVal, epsGrad, eta, gamma, m];
%    [r, p0, pdist] = ASHC(funflag, bdrge, parameter, p0, rotation_simplice, glb_min, fid, flag);
%    [r, pdist]
%
%    k = k+1
%    fprintf(fid, '\n\n ===== iteration : %d =====\n\n', k)
%
%end

fclose(fid);
toc



%SHC(funflag, bdrge, parameter, p0, rotation_simplice, glb_min, flag);
%ASHC_partsampling(funflag, bdrge, parameter, p0, rotation_simplice, glb_min, flag);
%profile off;
%profile viewer;
