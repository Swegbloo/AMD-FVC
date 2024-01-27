function out = Z_prime(k,m,z)
if k == 0
    out = m*sinh(m*z)/sqrt(N_const(m,k));
else
    out = -m*sin(m*z)/sqrt(N_const(m*k));
end