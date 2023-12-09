function out = phi_star(a,n,d)
out = 0;
l = n*pi/d;
if n == 0
    out = 0.25*a^3;
elseif n>0
    out = a*(l*a*besseli(0,l*a) - 2.0...
            *besseli(1,l*a))/(besseli(1,l*a)*l^2);
end
%disp(out)
end