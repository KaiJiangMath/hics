function rmat = rotation_matrix( t )

	n = length(t);

	if n == 1
		disp(' DIM 1 : there is no rotaion matrix ')
		return;
	end
	if n == 2
		rmat = [
			cos(t(1)) -sin(t(1))
			sin(t(1))  cos(t(1))
		];
	end
	if n > 2

		A = diag(ones(1,n));
		A(2,2) = cos(t(1));
		A(2,3) =-sin(t(1));
		A(3,2) = sin(t(1));
		A(3,3) = cos(t(1));

		rmat = A;
		for i = 2:1:n-1
			A = diag(ones(1,n));

			A(i-1,i-1) = cos(t(i));
			A(i-1,i+1) =-sin(t(i));
			A(i+1, i-1) = sin(t(i));
			A(i+1, i+1) = cos(t(i));

			rmat = rmat*A;
		end

		A = diag(ones(1,n));
		A(n-2, n-2) = cos(t(n));
		A(n-2, n-1) =-sin(t(n));
		A(n-1, n-2) = sin(t(n));
		A(n-1, n-1) = cos(t(n));

		rmat = rmat*A;


%        A = zeros(n, n, n);
%        for i = 1:1:n
%            A(:,:,i) = diag(ones(1,n));
%        end
%
%        A(2,2,1) = cos(t(1));
%        A(2,3,1) =-sin(t(1));
%        A(3,2,1) = sin(t(1));
%        A(3,3,1) = cos(t(1));
%
%        for i = 2:1:n-1
%            A(i-1,i-1,  i) = cos(t(i));
%            A(i-1,i+1,  i) =-sin(t(i));
%            A(i+1, i-1, i) = sin(t(i));
%            A(i+1, i+1, i) = cos(t(i));
%        end
%
%        A(n-2, n-2, n) = cos(t(n));
%        A(n-2, n-1, n) =-sin(t(n));
%        A(n-1, n-2, n) = sin(t(n));
%        A(n-1, n-1, n) = cos(t(n));
%        
%        rmat = A(:,:,1);
%        for i = 2:1:n
%            rmat = rmat*A(:,:,i);
%        end
	end
end
