function [z, gradx, grady, hesse] = objFun(x, y)

c1 = [1, -1];
c2 = -1*c1;

z1 = (x-c1(1))^2 + (y-c1(2))^2;
z2 = (x-c2(1))^2 + (y-c2(2))^2;
z = 0.5*max(z1, z2);

gradx = 1.0;
grady = 1.0;

gradxx = 1.0;
gradyy = 1.0;
gradxy = 1.0;
hesse = [gradxx, gradxy; gradxy, gradyy];
