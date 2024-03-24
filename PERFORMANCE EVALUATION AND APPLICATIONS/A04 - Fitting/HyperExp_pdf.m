function F = HyperExp_pdf(x, p)
	l1 = p(1);
	l2 = p(2);
	p1 = p(3);
	
	F = (x > 0) .* (p1 * l1 * exp(-l1*x) + (1-p1) * l2 * exp(-l2*x));
end