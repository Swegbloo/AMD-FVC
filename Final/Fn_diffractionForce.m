function [difForces, difMoment] =  Fn_diffractionForce(radius, depth, clearance, sigma, formulation,...
    modes)

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

difForces = 0.;
difMoment = 0.;

nEqs = 5;
% m_order = 1;

mroots = dispersion_free_surface_vMikeM(sigma,nEqs,depth);
nRoots = size(mroots,2);
m0 = -1i*mroots(1);
alpha = zeros(nEqs, 4);

alpha(:,3) = mroots(2:end).';


vEigNs = zeros(nRoots, 1);
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
jacobiSymbols = 2*ones(nEqs,1);
jacobiSymbols(1) = 1;

if (strcmp(formulation,'Nokob'))

eMatrix = zeros(nEqs, nEqs);


gVector = zeros(nEqs, 1);







onesExps = zeros(nEqs,1);

lambda = zeros(nEqs,1);



A0_star = zeros(nEqs,1);

psiFuns = zeros(nEqs,2);
psiStarFuns = zeros(nEqs,1);



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
        
% eMatrix
% gVector
% xVector
% return


for ik = 1:nEqs
    
    if ik == 1
    psiStarFuns(ik) = 0.5*radius^2;
    else
        psiStarFuns(ik) = radius*besseli(1,lambda(ik)*radius)/(lambda(ik)*besseli(0,lambda(ik)*radius));
    end
    
end

sum = 0;
 
for ik = 1:nEqs
    sum = sum + (-1)^(ik-1)*xVector(ik,1)*jacobiSymbols(ik)*psiStarFuns(ik);
end

difForces(1) = (sigma*depth/(pi*clearance*radius^2))*sum;
% abs(difForce)
% m0*radius
elseif (strcmp(formulation,'Garrett'))
    
    Lmatrix = zeros(nEqs, nEqs);
    for ik = 1:nEqs
        id = ik-1;
        for jk = 1:nEqs
            
            if (jk == 1)
            Lmatrix(ik,jk) = (-1)^id*m0*clearance*sinh(m0*clearance)/(sqrt(vEigNs(1))*((m0*clearance)^2 ... 
                + (id*pi)^2));
            else
                Lmatrix(ik,jk) = (-1)^id*mroots(jk)*clearance*sin(mroots(jk)*clearance)/(sqrt(vEigNs(jk))*...
                    ((mroots(jk)*clearance)^2 - (id*pi)^2));
            end
        end
    end
    
    m_order = 0;
    Fn_Garrett_diffraction_solver;
   
    
    %% heave diffraction force calculation
    sum = 0.;
    fVector = Lmatrix*xVector;
    for ik = 2:nEqs
        id = ik -1 ;
        gmn = Fn_gFun_Garrett(ik,id,radius, clearance, m0, alpha(id,3),m_order,1);
        sum = sum + (-1)^id*fVector(ik,1)/(gmn*(id*pi*radius/clearance)^2);
    end
    
    difForces(1) = 2*sigma*(0.5*fVector(1,1)+2*sum);
    

    
    if length(modes) > 1
    m_order = 1;
    Fn_Garrett_diffraction_solver;
    fVector = Lmatrix*xVector;
    
    % surge diffraction force calculation
    
    sum = 0.;
    
    for ik = 2:nEqs
        
        id = ik - 1;
        
        sum = sum + alpha(id,4)^(-0.5)*xVector(ik,1)*(alpha(id,3)*radius)^(-1)*(sin(alpha(id,3)*depth)...
            -sin(alpha(id,3)*clearance));
        
    end 
    

    
    difForces(2) = -2.0*1i*m0*depth*tanh(m0*depth)*(vEigNs(1)^(-0.5)*xVector(1,1)*(m0*radius)^(-1)*...
        (sinh(m0*depth)-sinh(m0*clearance)) + sum);
    
    %% pitch moment calculation
    sum = 0.;
    
    for ik = 2:nEqs
        sum = sum + (1.0/sqrt(vEigNs(ik)))*xVector(ik,1)*(radius*mroots(ik))^(-2)...
            *(mroots(ik)*(depth-clearance)*sin(mroots(ik)*depth)+cos(mroots(ik)*depth)-cos(mroots(ik)*clearance));
    end
    
    sideTorque = -2*1i*m0*depth*tanh(m0*depth)*(vEigNs(1)^(-0.5)*xVector(1,1)*(m0*radius)^(-2)*...
        (m0*(depth-clearance)*sinh(m0*depth)-cosh(m0*depth)+cosh(m0*clearance))+sum);
   
    
    sum = 0.;
    
    for ik = 2:nEqs
        id = ik - 1;
        gmn = Fn_gFun_Garrett(ik,id,radius, clearance, m0, alpha(id,3),m_order,1);
        gmn_inv = 1.0/gmn;
        sum = sum + (-1)^id*(id*pi*radius/clearance)^(-2)*fVector(ik,1)*(gmn_inv-1);
    end
    
%     sum
    
    bottomTorque = -2*1i*m0*depth*tanh(m0*depth)*(0.25*fVector(1,1)+2*sum);
    
    difMoment = sideTorque + bottomTorque; %
    
    

%     xVector = gmres(eMatrix, bVector, nEqs,1.0e-06, nEqs);
 
    
    



%     tempSum = sum;


%     
%     tempSum = sum;
    end
end

return