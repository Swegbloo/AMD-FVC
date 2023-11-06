function out = alpha(l,n,m,d,a)
out = zeros(n+1,1);
G = zeros(n+1,1);
E = zeros(n+1);
I = eye(n+1);
for x = 1:n+1
    G(x,1) = g_fun(n,x-1,m,d,l,a);
    for y = 1:n+1
            E(x,y) = e_fun(n,x-1,y-1,a,d,m,l);
    end
end
% disp(I-E);
%disp(I);
% disp(G);
out = (I-E)\G; %Given Matrix Equation: A-E*A=G =>A=G*(I-E)^(-1)
% disp(out)
end