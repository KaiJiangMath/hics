function pits = simplex_sampling( p0, r, s0)

	pits = r*s0 + repmat(p0, length(s0(:,1)), 1);

end
