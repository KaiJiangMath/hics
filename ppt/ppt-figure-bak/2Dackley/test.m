
%clear all

%DIM = 2;
%%%%%% generate regular simplex 
%simplex = regular_simplex(DIM);
%
%theta = pi/3;
%Rmat = [cos(theta) -sin(theta); sin(theta) cos(theta)];
%
%a1 = simplex(1,:);
%a2 = simplex(2,:);
%a3 = simplex(3,:);
%
%v1 = Rmat*a1'
%v2 = Rmat*a2'
%v3 = Rmat*a3'
%
%theta1 = acos(a2*a3'/(norm(a2,2)*norm(a3,2)))/2;
%theta2 = acos(a1*a3'/(norm(a1,2)*norm(a3,2)))/2;
%theta3 = acos(a1*a2'/(norm(a1,2)*norm(a2,2)))/2;

r0 = 1;
p0 = [0,0];
simplex = [
   1.000000000000000                   0
     -0.500000000000000   0.866025403784439
	   -0.500000000000000  -0.866025403784439 ];


r1 = 2.5;
p1 = [6.120000000000000   5.426000000000000];

r = r1;
p = p1
simplex = simplex*r+p
pause()

cicX = p(1)-r:0.005:p(1)+r;
pn = length(cicX);
for i = 1:pn
	cicY(i)  =  real(sqrt(r^2-(cicX(i)-p(1))^2))+p(2);
	cicY2(i) = -real(sqrt(r^2-(cicX(i)-p(1))^2))+p(2);
end

figure(1)
plot(p(1), p(2), 'ro')
hold on
plot(cicX, cicY,   'k-', 'LineWidth', 1.5)
plot(cicX, cicY2,  'k-', 'LineWidth', 1.5)

for i=1:1:3
	plot(simplex(i,1), simplex(i,2), 'bo');
end

axis square

