function out = A_star(m,d,a,k)
if k==0
    out = -a*sinh(m*d)/(2*d*m*(N_const(m,k))^0.5);
elseif k>0
    out = -a*sin(m*d)/(2*d*m*(N_const(m,k))^0.5);
end
end