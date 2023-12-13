function out = dbesselk(n,z,scal)
% this function calcualtes the derivative of besselk(n,z)
if nargin==3
    out = -(besselk(n+1,z,scal) + besselk(n-1,z,scal))/2;
else
    out = -(besselk(n+1,z) + besselk(n-1,z))/2;
end