function [y1,y2] = y_out(l,n,d,a,A,S)
y = zeros(1,l);
for i = 1:l
    for j = 2:n+1
        jj = j-1;
        y(1,i) = y(1,i) + ((-1)^jj)*A(j,i)*S(j,1)/(jj*pi)^2;
    end
    y(1,i) = 2*(0.25*(d/a)-(1/16)*(a/d)+(1/a)*(A(1,i)/4+2*(d/a)*y(1,i)));
end
y1 = real(y);
y2 = imag(y);
end
