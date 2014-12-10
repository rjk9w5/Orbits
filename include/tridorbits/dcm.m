function cart_vec = dcm(polar_vec,inclin,RAAN,arg_periap,nu)
	theta = arg_periap+nu;
	ccc = nan(3,3);

	ccc(1,1) = cos(RAAN)*cos(theta) - sin(RAAN)*sin(theta)*cos(inclin);
	ccc(1,2) = -cos(RAAN)*sin(theta) - sin(RAAN)*cos(theta)*cos(inclin);
	ccc(1,3) = sin(RAAN)*sin(inclin);

	ccc(2,1) = sin(RAAN)*cos(theta) + cos(RAAN)*sin(theta)*cos(inclin);
	ccc(2,2) = -sin(RAAN)*sin(theta) + cos(RAAN)*cos(theta)*cos(inclin);
	ccc(2,3) = -cos(RAAN)*sin(inclin);

	ccc(3,1) = sin(theta)*sin(inclin);
	ccc(3,2) = cos(theta)*sin(inclin);
	ccc(3,3) = cos(inclin);

	cart_vec = ccc*polar_vec;
end
