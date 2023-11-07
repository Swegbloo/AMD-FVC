function out = A_star(~,d,k,m)
if k==0
    out = (sinh(m)-sinh(m*d))/(m*(N_const(m,k))^0.5);
elseif k>0
    out = (sin(m)-sin(m*d))/(m*(N_const(m,k))^0.5);
end
end