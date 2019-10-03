function x = regular_simplex(n)

%*****************************************************************************

%  x = simplex_coordinates1 ( n );

  x(1:n,1:n+1) = 0.0;

  for i = 1 : n
%
%  Set X(I,I) so that sum ( X(1:I,I)^2 ) = 1.
%
    x(i,i) = sqrt ( 1.0 - sum ( x(1:i-1,i).^2 ) );
%
%  Set X(I,J) for J = I+1 to N+1 by using the fact that XI dot XJ = - 1 / N
%
    for j = i + 1 : n + 1
      x(i,j) = ( - 1.0 / n - ( x(1:i-1,i)' * x(1:i-1,j) ) ) / x(i,i);
    end

  end

  x = x';

end
