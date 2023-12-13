function out = cc_fun(n,m,d,k)
out = 0;
% if k == 0
%     out = (2*(-1)^(n)*sinh(m*d)/(m*d))/((1+((n*pi)/(m*d))^2)*sqrt(N_const(m,k)));
% elseif k >= 1
%     out = (2*sin(m*d-n*pi)*m*d)/((m*d-n*pi)*(m*d+n*pi)*sqrt(N_const(m,k)));
% end
lambda = n*pi/d;
if k == 0
   if n == 0
        out = sinh(m*d)/m;
        out = out*N_const(m,k)^-0.5;
   else
        out = m*sinh(m*d);
        out = out*(-1.0)^n;
        out = out*N_const(m,k)^-0.5;
        out = out/(lambda^2 + m^2);
   end
elseif k > 0
    if n == 0
        out = sin(m*d)/m;
        out = out*N_const(m,k)^-0.5;
    else
        out = m*sin(m*d);
        out = out/(m^2 - lambda^2);
        out = out*(-1.0)^n;
        out = out*N_const(m,k)^-0.5;
    end
end
%disp(out)
end