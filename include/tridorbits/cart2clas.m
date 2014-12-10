function [a,e,i,RAAN,arg_periap,M] = cart2clas(cartvec,mu)
	r_vec = cartvec(1:3);
	v_vec = cartvec(4:6);
	r = norm(r_vec,2);
	v = norm(v_vec,2);
	
	a = -mu/2*(v^2/2 - mu/r)^-1;
	h_vec = cross(r_vec,v_vec);
	h = norm(h_vec);

	e_vec = 1/mu*((v^2-mu/r)*r_vec - (dot(r_vec,v_vec))*v_vec);
	e = norm(e_vec,2);

	i = acos(h_vec(3)/h);

	tmp = cross([0 0 1],h_vec/h);
	a_hat = tmp/norm(tmp,2);

	RAAN = acos(dot(a_hat,[1 0 0]));

	arg_periap = acos(dot(a_hat,e_vec)/e);

	E = acos((1 - r/a)/e);
	% need to implement quad check to determine if body is approaching or leaving
	M = E = e*cos(E);
end