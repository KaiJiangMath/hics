function rsimplex = generate_rotation_simplex( DIM, n, simplex )

%    a1 = simplex(1,:);
%    a2 = simplex(2,:);
%    t = acos(a1*a2'/(norm(a1,2)*norm(a2,2)));
%    theta0 = t*ones(1, DIM);
%
%    rsimplex = cell(1, m);
%    
%    for mk = 1:1:m
%        rsimplex{mk} = cell(1,2^mk);
%        for nk = 1:2:2^mk
%            rsimplex{mk}{nk} = zeros(size(simplex));
%        end
%    end
%
%    for mk = 1:1:m
%        for nk = 1:2:2^mk
%            [ mk, nk ];
%            theta = theta0*nk/(2^mk);
%            mat = rotation_matrix( theta );
%            rsimplex{mk}{nk} = rotation_simplex( simplex, mat);
%        end
%    end

	m = n*DIM;
	rmat = generate_rotation_matrice(DIM, 0, m);

	for mk = 1:1:m
		rsimplex{mk} =  zeros(size(simplex));
	end
	for mk = 1:1:m
		mat = rmat{mk};
		rsimplex{mk} = rotation_simplex(simplex, mat);

	end

end
