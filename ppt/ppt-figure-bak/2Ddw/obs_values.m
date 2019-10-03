function obsValues = obs_values( pits, fun )
	
	obsValues = zeros(length(pits(:,1)),1);

	for i = 1:1:length(pits(:,1))
		obsValues(i) = fun_value(pits(i,:), fun);
	end

end
