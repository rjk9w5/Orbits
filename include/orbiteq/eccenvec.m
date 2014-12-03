function [e_vec] = eccenvec(mu,v_vec,r_vec)
v = norm(v_vec,2);
r = norm(r_vec,2);

e_vec = (1/mu) * ((v^2 - mu/r)*r_vec - dot(r_vec,v_vec)*v_vec);
end