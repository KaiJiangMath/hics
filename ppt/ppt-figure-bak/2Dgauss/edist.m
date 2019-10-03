function a = edist(p, p0)
	
	dim = length(p);
	distp = p0-p;
	distp = distp.^2;
	a = 0;
	for i = 1:1:dim
		a = a + distp(i);
	end
	a = sqrt(a);

end
