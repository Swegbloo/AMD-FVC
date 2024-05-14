function[pDiffMoment] = Fn_diffractionTorque(nEqs,radius,clearance,sigma)


eMatrix = zeros(nEqs, nEqs);




alpha_star = zeros(nEqs,1);

mroots = dispersion_free_surface_vMikeM(sigma,nEqs,1);
nRoots = size(mroots,2);
m0 = -1i*mroots(1);

alpha = zeros(nEqs,3);
alpha(:,3) = mroots(2:end).';

vEigNs = zeros(nRoots, 1);

lambda = zeros(nEqs,1);


jacobiSymbols = 2*ones(nEqs,1);
jacobiSymbols(1) = 1;

psiFuns = zeros(nEqs,2);

for ik = 1:nRoots
    if ik == 1
        vEigNs(ik) = 0.5 + sinh(2*m0)/(4*m0);

    else
        vEigNs(ik) = 0.5 + sin(2*mroots(ik))/(4*mroots(ik));

    end 

end

z0_das = m0*sinh(m0)/sqrt(vEigNs(1));

N0 = vEigNs(1);
alpha(:,4) = vEigNs(2:end).';


for ik = 1:nEqs
    
    lambda(ik) = (ik-1)*pi/clearance;
    
    
    cn0 = (clearance/2)*COMPCY(clearance,2,ik-1,m0,N0,alpha);
    alpha_star(ik) = 4*1i*cn0/(radius*dbesselh(1,1,m0*radius)*z0_das);
%     alpha_star(ik) = 8*1i*sqrt(vEigNs(1))*cn0/(radius*m0^2*(besselh(0,1,m0*radius)-besselh(2,1,m0*radius))*sinh(m0));

    if ik == 1
    psiFuns(ik,1) = radius;
    psiFuns(ik,2) = 1.;
    else
        psiFuns(ik,1) = besseli(1,lambda(ik)*radius);
        psiFuns(ik,2) = lambda(ik)*dbesseli(1,lambda(ik)*radius);
    end
    
    
end



for ik = 1:nEqs
    
    for jk = 1:nEqs
        sum = 0.;

        for pk = 1:nEqs
            r_pj = COMPRY(pk,radius,m0,alpha,1,1);
            r_pj_das = COMPRY(pk,radius,m0,alpha,2,1);
            if pk == 1
                r_pj_das = m0*r_pj_das;
            else
                r_pj_das = alpha(pk-1,3)*r_pj_das;
            end
            
            c_nj = (clearance/2)*COMPCY(clearance,pk,ik-1,m0,N0,alpha);
            c_lj= (clearance/2)*COMPCY(clearance,pk,jk-1,m0,N0,alpha);
            sum = sum + jacobiSymbols(jk)*r_pj*psiFuns(jk,2)*c_lj*c_nj/...
                (clearance*r_pj_das*psiFuns(jk,1));
        end
        
        eMatrix(ik,jk) = sum;
    end
    
end






% for ik = 1:nEqs
%     
%     alpha_star(ik) = -1i*2*pi*besselj(0,m0*radius)*sinh(m0*clearance)*(-1)^(ik-1)/...
%         (sinh(m0)*(m0^2 + lambda(ik)^2));
%     
% end

% gVector
% gVector = zeros(nEqs,1);
% 
% for ik = 1:nEqs
%     
%     sum = 0.0;
%     
%     for ij = 1:nEqs
%         cnj = (clearance/2)*COMPCY(clearance,ij,ik-1,m0,N0,alpha);
%         sum = sum + A0_star(ij)*cnj;
%     end
%     
%     gVector(ik) = -2.0*pi*sum + alpha_star(ik);
%     
% end



% xVector = gmres(eye(nEqs)-eMatrix, alpha_star, nEqs,1.0e-06, nEqs);
% alpha_star(2:end) = 0.;
xVector = (eye(nEqs)-eMatrix)\alpha_star;
% eMatrix


% A_kj
A_kj = zeros(nEqs,1);
for ik = 1:nEqs
    sum = 0.;
    for jk = 1:nEqs
        c_nj = (clearance/2)*COMPCY(clearance,jk,ik-1,m0,N0,alpha);
        eps_kn = jacobiSymbols(jk)*2.0/(2*pi*clearance);
        sum = sum + eps_kn*psiFuns(jk,2)*c_nj*xVector(jk,1)/psiFuns(jk,1);
    end
    
    A_kj(ik,1) = sum;
end

% A_kj
% return

sum = 0.;

for ik = 1:nEqs
    if ik == 1
        z_star = (cosh(m0*clearance) + m0*sinh(m0) - cosh(m0) - clearance*m0*sinh(m0*clearance))/(sqrt(N0)*m0^2);
        z_star_0 = z_star;
        r_pj_das = m0*COMPRY(ik,radius,m0,alpha,2,1);
    else
        mk = alpha(ik-1,3);
        Nk = alpha(ik-1,4);
        z_star = (cos(mk) + mk*sin(mk) - clearance*mk*sin(mk*clearance) - cos(mk*clearance))/(sqrt(Nk)*mk^2);
        r_pj_das = alpha(ik-1,3)*COMPRY(ik,radius,m0,alpha,2,1);
    end    
    
    r_pj = COMPRY(ik,radius,m0,alpha,1,1);

    sum = sum + A_kj(ik,1)*r_pj*z_star/r_pj_das;
end

sideTorque = 4.0*sigma*z_star_0/(pi*radius^3*dbesselh(1,1,m0*radius)*z0_das)...
    -1i*sigma*sum/radius^2; % 4.0*...

sum = 0.;
for ik = 1:nEqs
    if ik == 1
        psi_n_star = 0.25*radius^3;
    else
        
        psi_n_star = radius*(lambda(ik)*radius*besseli(0,lambda(ik)*radius) - 2.0...
            *besseli(1,lambda(ik)*radius))/(besseli(1,lambda(ik)*radius)*lambda(ik)^2);
    end
    
    
    sum = sum + xVector(ik,1)*jacobiSymbols(ik)*(-1)^(ik)*psi_n_star;
end



bottomTorque = -1i*sigma*sum/(pi*clearance*radius^3);

% pDiffMoment =  sideTorque - bottomTorque ; %-+ 

pDiffMoment(1) = sideTorque;
pDiffMoment(2) = bottomTorque;
pDiffMoment(3) = sideTorque + bottomTorque;

            