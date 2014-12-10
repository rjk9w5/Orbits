% bool quad_check: =0         -> departing
%			   	   =1 or null -> approaching

function [a,e,i,RAAN,arg_periap,E] = cart2clas(cartvec,mu,quad_check)
	quad_check = 1;
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
	if(quad_check)
			E = -E;
	end
end