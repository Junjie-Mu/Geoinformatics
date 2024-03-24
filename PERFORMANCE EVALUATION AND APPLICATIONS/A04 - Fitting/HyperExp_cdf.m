function F = HyperExp_cdf(x, p)
	l1 = p(1);
	l2 = p(2);
	p1 = p(3);
	
	F = max(0,1 - p1 * exp(-l1*x) - (1-p1) * exp(-l2*x));
end