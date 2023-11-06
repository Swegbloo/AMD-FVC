clc;
clear all;

nData = 20;
m = zeros(nData,1);
m(1,1) = 0.1;
global v;
v = m(1,1)*tanh(m(1,1));
for ik = 2:nData
    
    fun = @dispe;
    nIni = 0.51+(ik-2) ;
    x_0 = [nIni (nIni+0.9)]*pi;
    m(ik,1) = fzero(fun,x_0);
    
end
    