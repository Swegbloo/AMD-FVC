clear all;
clc;
a=0.5;
d=0.5;
n = 10;
l1=100;
z = zeros(1,l1);
m = zeros(n+1,l1);
alpha_vector = zeros(n+1,1);
if l1>1
m(1,:) = linspace(0,8,l1)/a; % index 1 implies 1st m_0
end
m(1,1) = 0.001;
x = m(1,:)*a; % x to be plotted
global v;
fun = @dispe;
for ii = 0:l1-1
    v = m(1,ii+1)*tanh(m(1,ii+1));
    % x_0 = [0 0.49*pi];
    % m(2,i+1) = fzero(fun,x_0);
    for j = 0:n-1
        x_0 = [(j+0.51)*pi (j+1.49)*pi]; % interval where the root lie
        m(j+2,ii+1) = fzero(fun,x_0);
    end
end

for l = 0:l1-1  
alpha_vector(:,1) = alpha(l,n,m(:,l+1),d,a);
for ii = 0:n
z(1,1+l) = z(1,1+l) + AR_fun(a,d,ii,m(:,l+1),n,alpha_vector)*A_star(a,d,ii,m(ii+1,l+1)); % A and R are combined
% ARV(i+1,1) = AR_fun(a,d,i,m(i+1,l+1),l);
% ASV(i+1,1) = A_star(a,d,i,m(i+1,l+1));
end

z(1,l+1) = z(1,l+1)*(-1/(a*(1-d)));

end

plot(x,imag(z));
% plot(x,imag(z));
hold on;
