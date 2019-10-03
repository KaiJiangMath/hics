function [y] = gauss(xx, h, s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gauss FUNCTION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = length(xx);

if (nargin < 3)
%    s(1) = 1.0;
%    s(2) = 0.5;
	s(1:d) = 1.0;
end
if (nargin < 2)
	h = 20;
end

sum = 0;
for ii = 1:d
	x2 = xx(ii)*xx(ii);
	sum = sum + x2/s(ii);
end
y = -h*exp(-sum);

end

