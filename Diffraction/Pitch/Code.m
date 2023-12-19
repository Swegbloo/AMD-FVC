% n = input('no. of roots (for accuracy) =');
% a = input('a =');
% l1 = input('l =');
% d = input('d =');
s = BasicClass;
n=4;
a=5;
l1=100;
d=0.1;
k = n; % so that the solution is unique
m = zeros(k+1,l1);
if l1>1
m(1,:) = linspace(0,4,l1)/a; % index 1 implies 1st m_0
end
m(1,1) = 0.001;
x = m(1,:)*a; % x to be plotted
S = zeros(n+1,1);
A = zeros(n+1,l1);
y1 = zeros(1,l1);
y2 = zeros(1,l1);
global v;
fun = @dispe;
for i = 1:l1
    v = m(1,i)*tanh(m(1,i));
    for j = 1:k
        x_0 = [(j-0.49)*pi j*pi]; % interval where the root lie
        m(j+1,i) = fzero(fun,x_0);
    end
end
for y = 1:n
    S(y,1) = s_fun(y-1,a,d);
end
for l = 1:l1
        A(:,l) = alpha(l,n,m,d,a);
end
%disp(A);
[y1,y2] = y_out(l1,n,d,a,A,S);
hold on;
% plot(x,y1/a);
plot(x,y2/a);
