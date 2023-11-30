function out = alpha_star(n,d,a)
out = 0;
l = n*pi/d;
if n == 0
    eps = 1;
    out = -eps/(2*d^2)*(a*d^3/3-a^3*d/4);
elseif n > 0
    eps = 2;
    out = -eps/(2*d^2)*(-a^3/(4*l)*sin(l*d)+a*sin(l*d)*d^2/l+2*d*a*cos(l*d)/l^2+2*a*sin(l*d)/l^3);
end
end