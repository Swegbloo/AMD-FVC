function [difTrq] =  Fn_diffractionTorque(radius, depth, clearance, sigma)

% clc;
% clear all;
% 
% radius = 0.2;
% clearance = 0.95;
% depth = 1.0;
% m0a = 0.2;
% acg = 9.81;
% 
% sigma = (m0a/radius)*tanh(m0a/radius);



nEqs = 5;

eMatrix = zeros(nEqs, nEqs);
alpha = zeros(nEqs, 4);

gVector = zeros(nEqs, 1);
xVector = zeros(nEqs, 1);



mroots = dispersion_free_surface_vMikeM(sigma,nEqs,depth);
nRoots = size(mroots,2);
% disp(nRoots);
% return;
m0 = -1i*mroots(1);



alpha = zeros(nEqs,3);
alpha(:,3) = mroots(2:end).';


onesExps = zeros(nEqs,1);

lambda = zeros(nEqs,1);

jacobiSymbols = 2*ones(nEqs,1);
jacobiSymbols(1) = 1;

A0_star = zeros(nEqs,1);
vEigNs = zeros(nRoots, 1);
psiFuns = zeros(nEqs,2);
psiStarFuns = zeros(nEqs,1);

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
    onesExps(ik) = (-1)^(ik-1);
    lambda(ik) = (ik-1)*pi/clearance;
    
    mj = mroots(ik);
    
    rj = COMPRY(ik,radius,m0,alpha,1,0);
    rj_das = COMPRY(ik,radius,m0,alpha,2,0);
    if ik == 1
        rj_das = m0*rj_das;
        A0_star(ik) = -1.0*rj*1i*m0*dbesselj(0,m0*radius)*(1+0.5*sinh(2*m0)/m0)/(2*vEigNs(1)*rj_das*z0_das);
    else
        rj_das = mj*rj_das;
    A0_star(ik) = -1.0*rj*1i*m0*dbesselj(0,m0*radius)*(mj*cosh(m0)*sin(mj) + m0*sinh(m0)*cos(mj))/...
        (rj_das*z0_das*sqrt(vEigNs(1)*vEigNs(ik))*(m0^2 + mj^2));
    

    end
end



for ik = 1: nEqs
    if ik == 1
        psiFuns(ik,1) = 1.0;
        psiFuns(ik,2) = 0.;
    else
        psiFuns(ik,2) = lambda(ik)^2*dbesseli(0,lambda(ik)*radius);
        psiFuns(ik,1) = lambda(ik)*besseli(0,lambda(ik)*radius);

%         psiFuns(ik,1) = lambda(ik)*besseli(0,lambda(ik)*radius)* dbseelk(0,1.0e-06);
    end
end

% psiFuns(:,1) = 1;
% psiFuns(:,2) = 1;




% alpha_star = 2*pi*1i*sinh(mroots(1)*depth)*jacobiSymbols.*onesExps.*besselj(0,mroots(1)*radius)./(sinh(mroots(1))*...
%     (mroots(1)+lambda.^2));

alpha_star = zeros(nEqs,1);

for ik = 1:nEqs
    
    alpha_star(ik) = -1i*2*pi*besselj(0,m0*radius)*sinh(m0*clearance)*(-1)^(ik-1)/...
        (sinh(m0)*(m0^2 + lambda(ik)^2));
    
end


% eMatrix

for ik = 1:nEqs
    
    for jk = 1:nEqs
        
        sum = 0.0;
        for ij = 1:nEqs
            
            
            rj = COMPRY(ij,radius,m0,alpha,1,0);
            rj_das = COMPRY(ij,radius,m0,alpha,2,0);
  
            
            if ij == 1
                rj_das = m0*rj_das;
            else
                rj_das = mroots(ij)*rj_das;
            end
            clj = (clearance/2)*COMPCY(clearance,ij,jk-1,m0,N0,alpha);
            cnj = (clearance/2)*COMPCY(clearance,ij,ik-1,m0,N0,alpha);
            sum = sum + (jacobiSymbols(jk)/clearance)*(rj/rj_das)*(psiFuns(jk,2)/psiFuns(jk,1))*cnj*clj;
        end
        
%         sum
%         return
        eMatrix(ik,jk) = sum;
    end
end




% gVector

for ik = 1:nEqs
    
    sum = 0.0;
    
    for ij = 1:nEqs
        cnj = (clearance/2)*COMPCY(clearance,ij,ik-1,m0,N0,alpha);
        sum = sum + A0_star(ij)*cnj;
    end
    
    gVector(ik) = -2.0*pi*sum + alpha_star(ik);
    
end

xVector = gmres(eye(nEqs)-eMatrix, gVector, nEqs,1.0e-06, nEqs);
        



for ik = 1:nEqs
    
    if ik == 1
    psiStarFuns(ik) = 0.5*radius^2;
    else
        psiStarFuns(ik) = radius*besseli(1,lambda(ik)*radius)/(lambda(ik)*besseli(0,lambda(ik)*radius));
    end
    
end

sum = 0;
sum2 = 0;
 
for ik = 1:nEqs
    sum = sum + (-1)^(ik-1)*xVector(ik,1)*jacobiSymbols(ik)*psiStarFuns(ik);
    sum2 = sum2 + alpha(1,ik-1)*eps(ik-1)*(-1)^(ik)*phi_star(radius,ik-1,clearance); %add fun alpha, eps and consequences
end

difTrq(1) = (sigma*depth/(pi*clearance*radius^2))*sum;
difTrq(2) = -1i*v/(pi*a^3*d)*sum2;
% abs(difForce)
% m0*radius
return