function [r,p0,pdist] = ASHCgauss(funflag, bdrge, parameter, p0, rotation_simplice, glb_min, fid, flag) 

DIM = length(p0); 

fname = sprintf('diary_%s_%dD_simplice.out', funflag, DIM);
diary(fname)

fprintf('===> START Stick Hill-Climbing for %d-D %s function \n\n\n', DIM, funflag)

%%%%%% generate regular simplex 
simplex = regular_simplex(DIM);

r0 = parameter(1);
epsVal = parameter(2);
epsGrad = parameter(3);
eta = parameter(4);
gamma = parameter(5);
m = parameter(6);
porig = p0;

%%%% maximum sampling points in each iteration
Nmax = m*DIM*(DIM+1);
gamma = min(epsVal, gamma);
r = r0;
iter= 0;
plotNum = 0;
grad = 1.0;
plotflag = flag(1);

%fname = sprintf('log_SHC_%s_%dD.txt', funflag, DIM);
%fid = fopen(fname, 'w');
valDiff = inf;
valMin = inf;
minPit = inf;
pdist = inf;

tic;

obsOrig = fun_value(p0, funflag);
funEvalu = 1;
obsPoints=[];
if plotflag == 1
	cdisPlot(p0, p0, [r0,r0], obsOrig, obsPoints, bdrge(1,:), ...
	bdrge(2,:), [0, funEvalu], plotNum, 1, 'Initialization', ...
	funflag, glb_min);
	pause(1);   close all
end

minimizer = fun_value(glb_min, funflag);

%    bdrge
%    pause()


while (r > epsVal)
%while (abs(valDiff) > epsVal)
%while (abs(minimizer-obsOrig) > epsVal)
%while (pdist > epsVal)
%    
	iter = iter+1

%    if(minbd == 0 )
%        r = r0/2;
%        p0 = porig;
%    else
%        r = min(r, min(r0, minbd))
%    end
%    r

	%%%% determin search radius 'r' by considering the boundary of search domain
	for i = 1:1:DIM
		mindist(i) = min(abs(p0(i)-bdrge(i,:)));
	end
	minbd = min(mindist(:));
	r = min(r, min(r0, minbd))
	

	%%%% gain sampling points
	obsPoints = simplex_sampling( p0, r, simplex );
	obsValues = obs_values(obsPoints, funflag);
	n = length(obsPoints(:,1));   %%%%%%%%%%

	%%%% count the number of function evaluation
	funEvalu = funEvalu + n;
	fprintf(fid, '%e %d %d %.10e %.15e  \n', r, n/(DIM+1), n, norm(p0-glb_min,2), obsOrig);
%    fprintf(fid, '%e %d %d %.10e %.10e  \n', r, n/(DIM+1), n, norm(p0-glb_min,2)/sqrt(DIM), obsOrig);

%    if norm(p0-glb_min,2) < r
%        r = r*eta;
%    end
%    r
	
	flagIndex = 0;   %%%  jump out to refinement
	
	[minPit minLoc] = min(obsValues);
	rho = (obsOrig-minPit)/abs(obsOrig);  %%% the ratio of function change

	%%%%%%%%% plot observed points %%%%%%%%%
	if plotflag == 1
		plotNum = plotNum + 1;
		cdisPlot(p0, p0, [r,r0], obsOrig, obsPoints, bdrge(1,:), ...
		bdrge(2,:), [n, funEvalu], plotNum, 0, 'sample',  ...
		funflag, glb_min);
		pause(0.5);   close all
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	valDiff = (minPit-obsOrig);
	pdist = norm(p0-glb_min,2);
%    pdist = norm(p0-glb_min,2)/sqrt(DIM);
	
	%    if ((minPit < obsOrig) & (rho > gamma))
	if (minPit < obsOrig)  %%  find the less value point
		obsOrig = minPit;
	
		p1 = p0;
		p0 = obsPoints(minLoc, :); 

		if plotflag == 1
			plotNum = plotNum +1
			cdisPlot(p1, p0, [r,r0], obsOrig, obsPoints, ...
			bdrge(1,:), bdrge(2,:), [0,funEvalu], plotNum, 1, ...
			'Optimal', funflag, glb_min);
			pause(0.5);  close all
		end

		%%%%  find the minimizer in the iteration
		fprintf('On the start sample points, finding the new point, val=    %e\n\n', obsOrig);
	else   %%%%%% refine
		obsPitsOrig = obsPoints;
		DofMin = 1;   %%% requir more sampling pits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		while DofMin == 1
			if n >= Nmax
				r = r*eta;
				n = 0;
				fprintf('Does not find a better point, shrinkage radius r : %e\n', r)
				DofMin = 0;
				flagIndex = 1;
				fprintf(fid, '\n');
			else
				sss = 0;
				
				for mk = 1:1:m*DIM

					sss = sss+1

				%%%%%%  predetermined rotation matrice or rotation simplex
					obsPoints = rotation_simplice{mk};
					
					obsPoints = simplex_sampling( p0, r, obsPoints );
					obsValues = obs_values(obsPoints, funflag);
					n = n + length(obsPoints(:,1));

					%%%% count the number of function evaluation
					funEvalu = funEvalu + length(obsPoints(:,1));
		%%%%%  to be improved, we do not need compare such many comparsion.
					[minPit, minLoc] = min(obsValues);
					minPit;
					valDiff = (minPit-obsOrig);
%                        pause()

					%%%%%%%%%%%%%%%%%%%%%%%%%
					if plotflag == 1
					plotNum = plotNum + 1;
					cdisPlot(p0, p0, [r,r0], obsOrig, ...
					obsPoints, bdrge(1,:), bdrge(2,:), ...
					[length(obsPoints(:,1)),funEvalu], ...
					plotNum, 0, 'Refinement', funflag, glb_min);
					pause(0.5);  close all
					end
					%%%%%%%%%%%%%%%%%%%%%%%%%

				% if ((minPit < obsOrig) & (rho > gamma))
					if (minPit < obsOrig)  %% find the less value point
						p1 = p0;
						p0 = obsPoints(minLoc, :);
						pit_value = [obsOrig,  minPit];

						obsOrig = minPit;

%                            fprintf(fid, '%e & %d & %d & %d & %.10e \\\\ \n', r, sss, (DIM+1)*sss, funEvalu, obsOrig);
						fprintf(fid, '%e %d %d %.10e %.15e \n', r, sss, (DIM+1)*sss, norm(p0-glb_min,2), obsOrig);
%                            fprintf(fid, '%e %d %d %.10e %.10e \n', r, sss, (DIM+1)*sss, norm(p0-glb_min,2)/sqrt(DIM), obsOrig);
						if plotflag == 1
						plotNum = plotNum + 1;
						cdisPlot(p1, p0, [r,r0], obsOrig, ... 
						obsPoints, bdrge(1,:), bdrge(2,:), ... 
						[0,funEvalu], plotNum, 1, 'optimal', ...
						funflag, glb_min);
						pause(0.5);   close all 
						end

						fprintf('Replace observed point: the new value : %e\n\n', minPit);
						flagIndex = 1;
						break;
					end
					if flagIndex == 1
						break;
					end
 
				end   %%%%% END FOR mk
				
				if flagIndex == 1
					break;
				end

			end  %% END if-else
			if(n>=(DIM+1))
				fprintf(fid, '%e %d %d %.10e %.15e \n', r, n/(DIM+1)-1, n-(DIM+1), norm(p0-glb_min,2), obsOrig);
%                fprintf(fid, '%e %d %d %.10e %.10e \n', r, n/(DIM+1)-1, n-(DIM+1), norm(p0-glb_min,2)/sqrt(DIM), obsOrig);
%                fprintf(fid, '%e & %d & %d & %d & %.15e \\\\ \n', r, n/(DIM+1)-1, n-(DIM+1), funEvalu, obsOrig);
			end

		end %  END WHILE
	end
end

runtime = toc


fprintf(fid, '\n distant between convergent pit and glb_min : %.15e', norm(p0-glb_min,2));
%fprintf(fid, '\n distant between convergent pit and glb_min : %.15e', norm(p0-glb_min,2)/sqrt(DIM));

format short; 
p0;
pdist = norm(p0-glb_min,2)

fprintf( fid, '\n CPU time : %.10f sec', runtime );
fprintf( fid, '\n funEvalu = %d', funEvalu );

obsOrig
funEvalu
%fclose(fid);

diary off

end
