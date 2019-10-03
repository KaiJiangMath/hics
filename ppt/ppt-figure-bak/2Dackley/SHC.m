function SHC(funflag, bdrge, parameter, p0, rotation_simplice, glb_min, flag) 

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
Nmax = (2^m)*(DIM+1)
gamma = min(epsVal, gamma);
r = r0;
iter= 0;
plotNum = 0;
grad = 1.0;
plotflag = flag(1);

fname = sprintf('log_SHC_%s_%dD.txt', funflag, DIM);
fid = fopen(fname, 'w');
valDiff = inf;
valMin = inf;
minPit = inf;

tic;

obsOrig = fun_value(p0, funflag)
funEvalu = 1;
obsPoints=[];
if plotflag == 1
	cdisPlot(p0, p0, [r0,r0], obsOrig, obsPoints, bdrge(1,:), ...
	bdrge(2,:), [0, funEvalu], plotNum, 1, 'Initialization', ...
	funflag, glb_min);
	pause(1);   close all
end


for i = 1:1:DIM

	d1 = p0(i) - bdrge(i,1);
	d2 = p0(i) - bdrge(i,2);

	if (d1 < 0) & (d2 > 0)
		break;
	end
end

flag = 1;

while (flag)
	%while (valDiff > eps)
    
	iter = iter+1

	%%%% gain sampling points
	obsPoints = simplex_sampling( p0, r, simplex );
	obsValues = obs_values(obsPoints, funflag);
	n = length(obsPoints(:,1));   %%%%%%%%%%

	%%%% count the number of function evaluation
	funEvalu = funEvalu + n;
%    fprintf(fid, '%d %e & %d & %d & %.10e & %.10e \n', iter, r, n/(DIM+1), n, norm(p0-glb_min,2), obsOrig);
	
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
	else    %%%%%  refine

		fprintf('Refinement \n\n');

		obsPitsOrig = obsPoints;
		DofMin = 1;   %%% requir more sampling pits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		while DofMin == 1
			if n >= Nmax
				flag = 0;
				DofMin = 0;
			else
				sss = 0;
				for mk = 1:1:m
					for nk = 1:2:2^mk
						sss = sss+1;

					%%%%%%  predetermined rotation matrice or rotation simplex
%                        rmat = rotation_matrice{mk}{nk};
%                        obsPoints = rotation_simplex( simplex, rmat );

						obsPoints = rotation_simplice{mk}{nk};
						
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
						disp('refine');
%                        fprintf(fid, '%e & %d & %d & %.10e & %.10e \\\\ \n', r, n/(DIM+1)-1, n-(DIM+1), norm(p0-glb_min,2), obsOrig);

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
%                            fprintf(fid, '%e & %d & %d & %.10e & %.10e \\\\ \n', r, sss, (DIM+1)*sss, norm(p0-glb_min,2), obsOrig);
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
					end  %%%%% END for

					if flagIndex == 1
						break;
					end
 
				end   %%%%% END FOR
				
				if flagIndex == 1
					break;
				end

			end  %% END if-else
			if(n>=(DIM+1))
				fprintf(fid, '%e  %d  %d  %.10e %.10e  \n', r, n/(DIM+1)-1, n-(DIM+1), norm(p0-glb_min,2), obsOrig);
			end

		end %  END WHILE

	end

end

runtime = toc

fprintf(fid, '\n distant between convergent pit and glb_min : %.15e', norm(p0-glb_min,2));
fprintf( fid, '\n CPU time : %.10f sec', runtime );

obsOrig
funEvalu
fclose(fid);

%                                cdisPlot(p0, p0, [r,r0], obsOrig, ... 
%                                obsPoints, bdrge(1,:), bdrge(2,:), ... 
%                                [0,funEvalu], plotNum, 1, 'optimal', ...
%                                funflag, glb_min);
%%                                pause(0.5);   close all 

diary off

end
