
%%%  寻找过圆点的直线与圆相交的点

function [y] = findpit(x0, x1, r)

	n = length(x0);
	y = zeros(2,n);

	t = r/norm(x0-x1,2);

	for i=1:1:n
		y(1,i) = x0(i) + t*(x1(i)-x0(i));
		y(2,i) = x0(i) - t*(x1(i)-x0(i));
	end
end
