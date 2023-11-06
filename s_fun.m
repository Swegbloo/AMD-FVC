function out = s_fun(n,a,d)
out = (n*pi/2)*(dbesseli(0,n*pi*a/d)/besseli(0,n*pi*a/d));
end