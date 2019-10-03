function rmat = generate_rotation_matrice( DIM, t0, m )
	
	rmat = cell(1,m);

	for mk = 1:1:m
		rmat{mk} = diag(ones(1, DIM));
	end

	if DIM == 1
		disp(' DIM 1 : there is no rotaion matrix ')
		return;
	end
	if DIM == 2
		rmat{1} = diag(ones(1,DIM));
		rmat{2} = [
			cos(pi/2) -sin(pi/2)
			sin(pi/2)  cos(pi/2)
		];
		rmat{3} = [
			cos(pi) -sin(pi)
			sin(pi)  cos(pi)
		];
		rmat{4} = [
			cos(-pi/2) -sin(-pi/2)
			sin(-pi/2)  cos(-pi/2)
		];
	end
	if DIM > 2
		rmat{1} = diag(ones(1,DIM));
		for i = 2:1:DIM-1
			A = diag(ones(1,DIM));
			A(i-1,i-1) = cos(pi/2);
			A(i-1,i+1) =-sin(pi/2);
			A(i+1,i-1) = sin(pi/2);
			A(i+1,i+1) = cos(pi/2);
			rmat{i} = A;

			A(i-1,i-1) = cos(-pi/2);
			A(i-1,i+1) =-sin(-pi/2);
			A(i+1,i-1) = sin(-pi/2);
			A(i+1,i+1) = cos(-pi/2);
			rmat{DIM+i} = A;
		end
		A = diag(ones(1,DIM));
		A(DIM-1, DIM-1) = cos(pi/2);
		A(DIM-1, DIM  ) =-sin(pi/2);
		A(DIM,   DIM-1) = sin(pi/2);
		A(DIM,   DIM)   = cos(pi/2);
		rmat{DIM} = A;

		A = diag(ones(1,DIM));
		A(1,1) = cos(pi);
		A(1,2) =-sin(pi);
		A(2,1) = sin(pi);
		A(2,2) = cos(pi);
		rmat{DIM+1} = A;

		A = diag(ones(1,DIM));
		A(DIM-1, DIM-1) = cos(-pi/2);
		A(DIM-1, DIM  ) =-sin(-pi/2);
		A(DIM ,  DIM-1) = sin(-pi/2);
		A(DIM,   DIM  ) = cos(-pi/2);
		rmat{m} = A;

		for i = 2:1:DIM-1
			A = diag(ones(1,DIM));
			A(i-1,i-1) = cos(pi/4);
			A(i-1,i+1) =-sin(pi/4);
			A(i+1,i-1) = sin(pi/4);
			A(i+1,i+1) = cos(pi/4);
			rmat{2*DIM+i} = A;

			A(i-1,i-1) = cos(-pi/4);
			A(i-1,i+1) =-sin(-pi/4);
			A(i+1,i-1) = sin(-pi/4);
			A(i+1,i+1) = cos(-pi/4);
			rmat{3*DIM+i} = A;
		end
		
%        for i = 2:1:DIM-1
%            A = diag(ones(1,DIM));
%            A(i-1,i-1) = cos(3*pi/4);
%            A(i-1,i+1) =-sin(3*pi/4);
%            A(i+1,i-1) = sin(3*pi/4);
%            A(i+1,i+1) = cos(3*pi/4);
%            rmat{4*DIM+i} = A;
%
%            A(i-1,i-1) = cos(-3*pi/4);
%            A(i-1,i+1) =-sin(-3*pi/4);
%            A(i+1,i-1) = sin(-3*pi/4);
%            A(i+1,i+1) = cos(-3*pi/4);
%            rmat{5*DIM+i} = A;
%        end
		
	end


%    rmat = cell(1, m);
%    for mk = 1:1:m
%        rmat{mk} = cell(1,2^mk);
%        for nk = 1:2:2^mk
%            rmat{mk}{nk} = zeros(DIM, DIM);
%        end
%    end
%
%    for mk = 1:1:m
%        for nk = 1:2:2^mk
%            t = t0*nk/(2^mk);
%            mat = rotation_matrix( t );
%            rmat{mk}{nk} = mat;
%        end
%    end

end
