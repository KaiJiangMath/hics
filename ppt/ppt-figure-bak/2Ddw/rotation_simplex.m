function pits = rotation_simplex( s0, rmat )

%    rmat = rotation_matrix(t);

	for i = 1:1:length(s0(:,1))
		pits(i,:) = rmat*(s0(i,:)');
	end

end
