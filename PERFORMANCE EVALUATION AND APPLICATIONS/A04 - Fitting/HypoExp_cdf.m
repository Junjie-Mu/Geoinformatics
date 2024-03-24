function F = HypoExp_cdf(x, p)
	l1 = p(1);
	l2 = p(2);
	
	F = (x>0) .* min(1,max(0,1 - l2/(l2-l1) * exp(-l1*x) + l1/(l2-l1) * exp(-l2*x)));
end