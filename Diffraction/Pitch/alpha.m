function out = alpha(p,n,m,d,a)
out = zeros(n+1,1);
G = zeros(n+1,1);
E = zeros(n+1);
I = eye(n+1);

    for x = 1:n+1
            G(x,1) = alpha_star(p,x-1,d,a,m(1));
        for y = 1:n+1
            E(x,y) = e_fun_roll(p,n,x-1,y-1,a,d,m);
        end
    end
    %disp(E)
     % disp(I-E);
    %disp(I);
    
     %disp(G);
    out = inv(I-E)*G; %Given Matrix Equation: A-E*A=G =>A=G*(I-E)^(-1)
    % disp(out)
end
