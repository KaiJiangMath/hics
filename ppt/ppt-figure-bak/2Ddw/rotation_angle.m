function theta = rotation_angle( s0 )

	a1 = s0(1,:);
	a2 = s0(2,:);
	
	theta = acos(a1*a2'/(norm(a1,2)*norm(a2,2)));

end
