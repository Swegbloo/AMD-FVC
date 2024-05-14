function eh = energy_harnessed(b_3,x_3,b_5,x_5,b_2,x_2)
B = zeros(6);
X = zeros(6,1);
%disp(b_3);
%disp(x_3); 
for i = 1:6
    for j = 1:6
        if (i == 3) && (j == 3)
           B(i,j) = b_3;
           X(i,1) = x_3;
        elseif (i == 5) && (j == 5)
           B(i,j) = b_5;
           X(i,1) = x_5;
        elseif (i==2)&&(j==2)
           B(i,j) = b_2;
           X(i,1) = x_2;
    end
end
if det(B) == 0
    for i=1:6
        if B(i,i) ~= 0
            eh = (1/8)*transpose(X)*X/B(i,i);
        end
    end
else
disp(B);
Bi = inv(B);
% Bi(1,1) = 0;
% Bi(2,2) = 0;
% Bi(4,4) = 0;
% Bi(6,6) = 0;
% Bi(5,5) = 0;
%disp((Bi));
eh = (1/8)*transpose(X)*Bi\X;
end
end