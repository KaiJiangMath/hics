function value = fun_value( pit, fun )
	
	if strcmp(fun, 'ackley')
		value = ackley(pit);
	end
	if strcmp(fun, 'bukin6')
		value = bukin6(pit);
	end
	if strcmp(fun, 'drop')
		value = drop(pit);
	end
	if strcmp(fun, 'egg')
		value = egg(pit);
	end
	if strcmp(fun, 'gauss')
		value = gauss(pit);
	end
	if strcmp(fun, 'dw')
		value = dw(pit);
	end
	if strcmp(fun, 'camelfun')
		value = camel(pit);
	end
	if strcmp(fun, 'zakharov')
		value = zakharov(pit);
	end
	if strcmp(fun, 'woods')
		value = woods(pit);
	end
	if strcmp(fun, 'arwhead')
		value = arwhead(pit);
	end
	if strcmp(fun, 'chrosen')
		value = chrosen(pit);
	end
	if strcmp(fun, 'BDQRTIC')
		value = BDQRTIC(pit);
	end
	if strcmp(fun, 'powell')
		value = powell(pit);
	end
	if strcmp(fun, 'spheref')
		value = powell(pit);
	end
	
end
