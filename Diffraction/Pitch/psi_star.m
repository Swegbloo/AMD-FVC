function out = psi_star(a,n,d)
out = 0;
l = n*pi/d;
if n == 0
    out = a^3/(p+3);
elseif n>0
    for i=1:50
        out = out + (n*pi*a/d)^(2*i)/(factorial(i)*factorial(p+i)*4^i*(p+2*i+3));
    end
    out = out*((n*pi*a/d)^(p+3))/(2^p);
end
%disp(out)
end
