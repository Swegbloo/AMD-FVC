function eh = energy_harnessed(b_3,x_3,b_5,x_5,b_2,x_2,b_51)
B = zeros(3);
X = zeros(3,1);
%disp(b_3);
%disp(x_3); 
for i = 1:3
    for j = 1:3
        if (i == 1) && (j == 1)
           B(i,j) = b_2;
           X(i,1) = x_2;
        elseif (i == 2) && (j == 2)
           B(i,j) = b_3;
           X(i,1) = x_3;
        elseif (i==3) && (j==3)
           B(i,j) = b_5;
           X(i,1) = x_5;
        elseif (i==1)&&(j==3)
            B(i,j)=b_51;
            B(j,i)=b_51;
        end
    end
end
    eh = 0;
if b_51==0
    for i=1:3
        if B(i,i) ~= 0
            eh = eh + (1/8)*abs(X(i,1))^2/B(i,i);
        end
    end
elseif b_51~=0
Bi = inv(B);

if b_51~=0
    disp(B);
    disp(inv(B));
end
eh = (1/8)*X'*Bi*X;
end
end