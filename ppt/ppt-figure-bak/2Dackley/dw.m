function [y] = dw(xx, c1, c2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Dennis-Woods FUNCTION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = length(xx);

if (nargin < 3)
	c2 = [-1, 1];
end
if (nargin < 2)
	c1 = [1, -1];
end

if d == 2
%    s1 = norm(xx-c1, 2);
%    s2 = norm(xx-c2, 2);
%    y = 0.5*max(s1, s2);

	x = xx(1);
	y = xx(2);
	z1 = (x-c1(1))^2 + (y-c1(2))^2;
	z2 = (x-c2(1))^2 + (y-c2(2))^2;
	y  = 0.5*max(z1, z2);
else
	disp('dimension is wrong');
	return;
end


end
