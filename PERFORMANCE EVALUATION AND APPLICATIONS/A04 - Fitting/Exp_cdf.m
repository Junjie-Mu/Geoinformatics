function F = Exp_cdf(x, p)
	l = p(1);
	
	F = max(0,1 - exp(-l*x));
end