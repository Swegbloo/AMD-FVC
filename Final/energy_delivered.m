function ed = energy_delivered(h,m,A)
 % A will cancel on ratio P_h/P_d
g = 9.81;
p = 1025;
cg = sqrt(g*h)/2*(tanh(m)+m*(sech(m))^2)/sqrt(m*tanh(m));
ed = 0.5*p*g*(A^2)*cg;
end