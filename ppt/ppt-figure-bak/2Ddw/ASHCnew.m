function [r,p0,pdist] = ASHCnew(funflag, bdrge, parameter, p0, rotation_simplice, glb_min, fid, flag) 


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

minimizer = fun_value(glb_min, funflag);

while (r > epsVal)

	iter = iter+1

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

	y = findpit(p0, glb_min, r);
	fprintf('%f %.15e %.15e  \n', r, norm(p0-glb_min,2), obsOrig);
	fprintf('[%.15e %.15e %.10e] \n\n\n', norm(y(1,:)-glb_min,2), norm(y(2,:)-glb_min,2), norm(y(1,:)-y(2,:),2) );
% 	pause();
	
	flagIndex = 0;   %%%  jump out to refinement
	
	[minPit minLoc] = min(obsValues);
	rho = (obsOrig-minPit)/abs(obsOrig);  %%% the ratio of function change

	valDiff = (minPit-obsOrig);
	pdist = norm(p0-glb_min,2);
%    pdist = norm(p0-glb_min,2)/sqrt(DIM);
	
	%    if ((minPit < obsOrig) & (rho > gamma))
	if (minPit < obsOrig)  %%  find the less value point
		obsOrig = minPit;
	
		p1 = p0;
		p0 = obsPoints(minLoc, :); 

%        pdist = norm(p0-p1,2)
%        pause()

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

						if (minPit < obsOrig)  %% find the less value point
							p1 = p0;
							p0 = obsPoints(minLoc, :);

							pit_value = [obsOrig,  minPit];

							obsOrig = minPit;
							fprintf(fid, '%e %d %d %.10e %.15e \n', r, sss, (DIM+1)*sss, norm(p0-glb_min,2), obsOrig);

							y = findpit(p0, glb_min, r);
							fprintf('%f %.15e %.15e  \n', r, norm(p0-glb_min,2), obsOrig);
							fprintf('[%.15e %.15e] \n\n\n', norm(y(1,:)-glb_min,2), norm(y(2,:)-glb_min,2) );
% 							pause();

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
				fprintf(fid, '%e %d %d %.10e %.15e \n', r, n/(DIM+1)-1, n-(DIM+1), norm(p0-glb_min,2), obsOrig);
				y = findpit(p0, glb_min, r);
				fprintf('%f %.15e %.15e  \n', r, norm(p0-glb_min,2), obsOrig);
				fprintf('[%.15e %.15e] \n\n\n', norm(y(1,:)-glb_min,2), norm(y(2,:)-glb_min,2) );
				%pause();
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
