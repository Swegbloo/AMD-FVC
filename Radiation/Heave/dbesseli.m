function out = dbesseli(n,z,scal)
% this function calcualtes the derivative of besselk(n,z)
if nargin==3
    out = (besseli(n+1,z,scal) + besseli(n-1,z,scal))/2;
else
    out = (besseli(n+1,z) + besseli(n-1,z))/2;
end
