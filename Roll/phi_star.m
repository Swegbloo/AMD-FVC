function out = phi_star(a,n,d)
out = 0;
l = n*pi*a/d;
if n == 0
    out = a^4/4;
elseif n>0
    out = (l*a)*(l*a*besseli(0,l*a)-2*besseli(1,l*a))/(dbesseli(1,l*a)*l^3);
end
end