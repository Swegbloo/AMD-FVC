
function [FZ_SING,EXTC,INTC] = cylinderSolverHeave_vYeung(a,h,mroots)     %(a,d,dpt,mroots,g)



% a=0.1;d=1.0;h=0.95;ka=1.1;
mop = 0;
ENDN = numel(mroots)-1;
INDN = 0:1:ENDN;

NT = size(INDN,2);


AMAT(1:NT,1:NT) = 0.0;

BVEC(1:NT,1) = 0.0;

SVEC(1:NT,1) = 0.0;

% ALPHAN(1:NT,1) = 0.0;



n_alp = ENDN;

alpha(1:n_alp,4)=0.0;


dn = h;

% vars = defGlobalVars;
% dpt = vars.waterDepth;

mroots = mroots.';

waveNumber = -1i*mroots(1,1);

alpha_0 = waveNumber;
arg = 2.0*alpha_0;
N0 = 0.5*(1.0 + sinh(arg)/arg);
% 
for ik=2:NT
    
	n = 0.51+(ik-2)*1.0;
    alpha(ik-1,1) = n*pi;
    alpha(ik-1,2) = (n+ 0.9)*pi;
    
    alpha(ik-1,3) = mroots(ik,1);
%     rm = real(mroots(ik,1));im=imag(mroots(ik,1));
%     alpha(ik-1,3) = rm - 1i*im;
    
	arg = 2.0*alpha(ik-1,3);
    alpha(ik-1,4) = 0.5*(1.0 + sin(arg)/arg);
    
end

% alpha = mroots(2:N,1);

% A_star(1:NT,1) = 0.0;
% 
% for in2=1:NT
%     
%     in = INDN(in2);
%     
% 	if(in==0)
%             
%             %%%% for heave %%%%
%             A_star = (a/(2.0*dn))*(1.0/alpha_0)*N0^-0.5;
%             A_star = -1.0*A_star*sinh(alpha_0*dn);
% 
%             %%%% end for heave %%%%
%             
%             %%%% for roll %%%%
% %             arg = alpha_0*dn;
% %             fac1 = 0.375*a^2 - 1.5*dn^2 + dn;
% %             fac2 = cosh(alpha_0) - sinh(arg)/arg - 2.0*cosh(arg);
% %             A_star(in2,1) = fac1*sinh(arg)/arg - fac2/alpha_0^2;
% %             A_star(in2,1) = A_star(in2,1)*N0^-0.5;
%             
%             %%%% end for roll %%%
%             
% %             mk = alpha_0;
% 	else
%             
%             mk = alpha(in,3);Nk = alpha(in,4);
%             
%             %%% for heave %%%
%             A_star = (a/(2.0*dn))*(1.0/alpha(ik-1,3))*alpha(ik-1,4)^-0.5;
%             A_star = -1.0*A_star*sin(alpha(ik-1,3)*dn);
%             %%% end for heave %%%
%             
%             %%% for roll %%%
%             
% 
% %             arg = mk*dn;
% %             fac1 = 0.375*a^2 - 1.5*dn^2 + dn;
% %             fac2 = 2.0*cos(arg) - cos(mk) - sin(arg)/arg;
% %             A_star(in2,1) = fac1*sin(arg)/arg - fac2/mk^2;
% %             A_star(in2,1) = A_star(in2,1)*Nk^-0.5;
%             
%             %%% end for roll %%%
%             
%             
% 
% 	end
% end


%%% begin of Yeung's formulation %%%

for in2 = 1:NT
    
    in = INDN(in2);
    
    %%% for heave %%%
    
    if (in==0)
        
        alpha_star = dn*(1.0/3.0 - 0.5*(a/dn)^2);

        % disp(alpha_star);
    else
        
%         lambda = in*pi/dn;
%         
%         arg = lambda*dn;
%         I1 = sin(arg)*(arg)^2 + 2.0*arg*cos(arg) - 2.0*sin(arg);
%         I1 = I1/lambda^3;
%         
%         I2 = 0.5*a^2/lambda;
%         I2 = I2*sin(arg);
%         
%         alpha_star = (I1 - I2)/dn^2;

        alpha_star = 2.0*dn/(in*pi)^2;
        alpha_star = alpha_star*(-1.0)^in;
        % disp(alpha_star);
    end
    
    %%% end for heave %%%%
    
    %%% for pitch %%%
    
%     if(in==0)
% %         mk = alpha_0;
% 
% %         alpha_star = -1.0*a*dn*(0.3333 - 0.25*(a/dn)^2);
% 
%           alpha_star = 0.0;
%     else
%         
% %         mk = alpha(ik-1,3);
% %         alpha_star = (2.0*a*dn*(-1.0)^(in+1))/(in*pi)^2;
% 
%           alpha_star = 0.0;
%         
%     end
    
    %%% end for pitch %%%
    
    sum = 0.0;
    
   for ik=1:n_alp
       
        CNK = COMPCY(dn,ik,in,alpha_0,N0,alpha);
        
        option = 1;
        RK = COMPRY(ik,a,alpha_0,alpha,option,mop);
        
        option = 2;
        RK_das = COMPRY(ik,a,alpha_0,alpha,option,mop);
        
        if(ik==1)

            mk = alpha_0;
            arg = alpha_0*dn;
%             fac1 = 0.375*a^2 - 1.5*dn^2 + dn;
%             fac2 = cosh(alpha_0) - sinh(arg)/arg - 2.0*cosh(arg);
%             A_star = fac1*sinh(arg)/arg - fac2/alpha_0^2;
%             A_star = A_star*N0^-0.5;

%             A_star = (a/(2.0*dn))*(1.0/alpha_0)*N0^-0.5;
%             A_star = -1.0*A_star*sinh(alpha_0*dn);

%             A_star = sinh(mk) - sinh(arg);
%             A_star = (N0^-0.5)*A_star/mk;

            A_star = (a/(2.0*dn))*(1.0/alpha_0)*N0^-0.5;
            A_star = -1.0*A_star*sinh(alpha_0*dn);
%             disp(A_star);

        else
            mk = alpha(ik-1,3);Nk = alpha(ik-1,4);
            

        
        
            arg = mk*dn;
%             fac1 = 0.375*a^2 - 1.5*dn^2 + dn;
%             fac2 = 2.0*cos(arg) - cos(mk) - sin(arg)/arg;
%             A_star = fac1*sin(arg)/arg - fac2/mk^2;
%             A_star = A_star*Nk^-0.5;

%             A_star = (a/(2.0*dn))*(1.0/alpha(ik-1,3))*alpha(ik-1,4)^-0.5;
%             A_star = -1.0*A_star*sin(alpha(ik-1,3)*dn);

%             A_star = sin(mk) - sin(arg);
%             A_star = (Nk^-0.5)*A_star/mk;

            A_star = (a/(2.0*dn))*(1.0/mk)*Nk^-0.5;
            A_star = -1.0*A_star*sin(mk*dn);
            % disp(A_star);
        
        end
        

        sum = sum + (CNK/mk)*(RK/RK_das)*A_star;
       
   end
        
%         in4 = INDN(in3);
        


        
	
        
        


    
    BVEC(in2,1) = sum - alpha_star; %
    
    SY = COMPSY(in,a,dn);

%     SY = compsy_roll(in,a,dn,alpha);
        
    SVEC(in2,1) = SY;
    
end  

% stop;

for ini2 = 1:NT
    
    ini = INDN(ini2);
    
    for inj2 = 1:NT
        
        
        inj = INDN(inj2);
                    
        sum = 0.0;
        
        for ik=1:n_alp
            
            if(ik==1)
                mk = alpha_0;
            else
                mk = alpha(ik-1,3);
            end
            

            CNK = COMPCY(dn,ik,ini,alpha_0,N0,alpha);
        

            CJK = COMPCY(dn,ik,inj,alpha_0,N0,alpha);
            
            
            option = 1;
            RK = COMPRY(ik,a,alpha_0,alpha,option,mop);
        
            option = 2;
            RK_das = COMPRY(ik,a,alpha_0,alpha,option,mop);
            
            sum = sum + (CNK*CJK/mk)*(RK/RK_das);
        
        end
        

        SY = SVEC(inj2,1);
        e = sum*SY;
        
        if (ini2==inj2)
            
            AMAT(ini2,inj2) = 1.0 - e;
            
        else
            
            AMAT(ini2,inj2) = -1.0*e;
            
        end
        
    end
    
end

% AMAT
% BVEC

%% end of Yeung's formulation %%%


%% begin of MH Nokob's formulation %%%
% epsn(1:NT,1) = 0.0;
% 
% A_k_star(1:NT,1) = 0.0;
% alpha_n_star(1:NT,1) = 0.0;
% 
% CNKM(1:NT,1:NT) = 0.0;
% RKD(1:NT,2) = 0.0;
% PSID(1:NT,2) = 0.0;
% 
% m0 = alpha_0;
% 
% for ii=1:NT
%     indi = INDN(ii);
%     if(indi == 0)
%         epsn(ii,1) = 1.0;
%         RKD(ii,1) = besselh(1,1,m0*a);
%         RKD(ii,2) = m0*dbesselh(1,1,m0*a);
%         
%         PSID(ii,1) = a;
%         
%         PSID(ii,2) = 1.0;
%     else
%         epsn(ii,1) = 2.0;
%         mk = alpha(indi,3);
%         RKD(ii,1) = besselk(1,mk*a);
%         RKD(ii,2) = mk*dbesselk(1,mk*a);
% 
%         PSID(ii,1) = besseli(1,indi*pi*a/dn);
%         PSID(ii,2) = indi*pi*dbesseli(1,indi*pi*a/dn)/dn;
% %         PSID(ii,2) = 0.5*pi*dbesseli(1,mk*a);
% 
%     end
%     
% end
% 
% for ii=1:NT
%     
%     indi = INDN(ii);
%     
%     if(indi==0)
%         z = -1/6.0 + 0.125*(a/dn)^2;
%         z = z*a*dn*epsn(ii,1);
%         alpha_n_star(ii,1) = z;
%         
%     else
%         
%         z = -1.0*epsn(ii,1)*a*dn*(-1.0)^indi;
%         z = z/(indi*pi)^2;
%         
%         alpha_n_star(ii,1) = z;
%         
%     end
%     
% end
% 
% for ii=1:NT
%     
%     indi = INDN(ii);
%     
%     if(indi==0)
%         
% %         m0 = alpha_0;
%         
% %         sum1 = 0.375*a^2 - 1.5*dn^2 - 1/m0^2;
% %         
% %         sum2 = cosh(m0) - 2.0*cosh(m0*dn) - m0*sinh(m0);
% % 
% % %         sum1 = 0.375*a^2 - 1.5*dn^2 + dn;
% % %         sum2 = cosh(m0) - sinh(m0*dn)/(m0*dn) - 2.0*cosh(m0*dn);
% %         
% %         arg = m0*dn;
% %         
% %         A_k_star(ii,1) = sum1*sinh(arg)/arg - sum2/m0^2;
% %         
% %         A_k_star(ii,1) = A_k_star(ii,1)*N0^-0.5;
% 
%         A_k_star(ii,1) = (16.0*m0*dn*cosh(m0*dn) + 8.0*m0^2*dn*sinh(m0) - 8.0*m0*dn*cosh(m0)...
%             -(12.0*m0^2*dn^2 - 3.0*m0^2*a^2 + 8.0)*sinh(m0*dn))/(8.0*m0^3*dn*sqrt(N0));
%         
%     else
%         
%         mk = alpha(indi,3);
%         Nk = alpha(indi,4);
% %         
% %         sum1 = 0.375*a^2 - 1.5*dn^2 + 1/mk^2;
% %         
% %         sum2 = 2.0*cos(mk*dn) - cos(mk) - mk*sin(mk);
% % 
% % %         sum1 = 0.375*a^2 - 1.5*dn^2 + dn;
% % %         sum2 = 2.0*cos(mk*dn) - cos(mk) - sin(mk*dn)/(mk*dn);
% %         
% %         arg = mk*dn;
% %         
% %         A_k_star(ii,1) = sum1*sin(arg)/arg - sum2/mk^2;
% %         
% %         A_k_star(ii,1) = A_k_star(ii,1)*Nk^-0.5;
% 
%         A_k_star(ii,1) = (-16.0*mk*dn*cos(mk*dn) + 8.0*mk^2*dn*sin(mk) + 8.0*mk*dn*cos(mk)...
%             -(12.0*mk^2*dn^2 - 3.0*mk^2*a^2 - 8.0)*sin(mk*dn))/(8.0*mk^3*dn*sqrt(Nk));
%         
%     end
%     
% end
% 
% for ii=1:NT
%     indi = INDN(ii);
%     
%     for jj=1:NT
%         
%         indj = INDN(jj);
%         
%         if(indj == 0)
% 
%             
%             option = 1;
%             
%             
%             Nn = N0;mn = m0;ik2=indi;
%             
%             CNKM(ii,jj) = compc_mhn(Nn,mn,ik2,dn,option);
%         else
%             mk = alpha(indj,3);
%             Nk = alpha(indj,4);
%             
%             option = 2;
%             Nn = Nk;mn = mk;ik2=indi;
%             CNKM(ii,jj) = compc_mhn(Nn,mn,ik2,dn,option);
%         end
%         
%     end
%     
% end
% 
% for ii=1:NT
%     
%     indi = INDN(ii);
%     
% %     if(ik==0)
% %     
% %             m0 = alpha_0;
% %             
% %             RK = besselh(1,1,m0*a);
% %             
% %             RK_das = m0*dbesselh(1,1,m0*a);
% %             
% %     else
% %             
% %             mk = alpha(ik,3);
% %             Nk = alpha(ik,4);
% %             
% %             RK = besselk(1,mk*a);
% %             RK_das = mk*dbesselk(1,mk*a);
% %             
% %     end
%     
%     sum = 0.0;
%     
%     for jj=1:NT
%         
% 
%         
%         sum = sum + A_k_star(jj,1)*RKD(jj,1)*CNKM(ii,jj)/RKD(jj,2);
%         
%     end
%     
%     BVEC(ii,1) = epsn(ii,1)*sum/dn - alpha_n_star(ii,1);
%     
% end
% % 
% for ii = 1:NT
%     
%     indi = INDN(ii);
%     
%     for jj = 1:NT
%         
%         
%         indj = INDN(jj);
%                     
%         sum = 0.0;
%       
%             
%         
%         for ik=1:NT
%             
% 
% 
%             
%             sum = sum + PSID(jj,2)*CNKM(ii,ik)*CNKM(jj,ik)*RKD(ik,1)/(PSID(jj,1)*RKD(ik,2));
%             
% 
%         
%         end
%         
% 
%         e = epsn(ii,1)*sum/dn;
%         
%         if (ii==jj)
%             
%             AMAT(ii,jj) = 1.0 - e;
%             
%         else
%             
%             AMAT(ii,jj) = -1.0*e;
%             
%         end
%         
%     end
%     
% end
%     
%     
% 
% 
% %% end of MH Nokob's formulation %%%
% 



% RESTART  = 20;
% TOL = 1.0e-06;
% MAXIT = 20;
ALPHAN= AMAT\BVEC; %gmres(AMAT,BVEC,RESTART,TOL,MAXIT)


% ALPHAN(1,1) = 11.83429 + 1i*2.61280;
% ALPHAN(2,1) = 0.22991 - 1i*0.02399;
% ALPHAN(3,1) = -0.07093 + 1i*0.00776;
% ALPHAN(4,1) = 0.03547 - 1i*0.00396;

% % FORCE EVALUATION %

%%% for heave %%%

sum = 0.0;

for in2=2:NT
    
    in = INDN(in2);
    
    lambda = in*pi/dn;
    add = ALPHAN(in2,1)*SVEC(in2,1)/(in*pi)^2;

    sum = sum + add*(-1)^in;    %
    
end

FZ = 0.25*(dn/a) - a/(16.0*dn) + (0.25*ALPHAN(1,1) + 2.0*(dn/a)*sum)/a;

FZ = 2.0*FZ;

% AM = real(FZ)/a
% DC = imag(FZ)/a

% AM(im,1) = real(FZ)/a;
% DC(im,1) = imag(FZ)/a;

FZ_SING = FZ;



%%% end for heave %%%

%%% for roll %%%

%% begin of Yeung's formulation %%%

% Ek(1:n_alp,1) = 0.0;
% 
% for ik=1:n_alp
%     
%     if(ik == 1)
%         
%         Ek(ik,1) = (dn-1)*sinh(alpha_0*dn)/alpha_0 + (cosh(alpha_0) - cosh(alpha_0*dn))/alpha_0^2;
%         Ek(ik,1) = Ek(ik,1)*N0^-0.5;
%         
%     else
%         mk = alpha(ik-1,3);nk = alpha(ik-1,4);
%         
%         Ek(ik,1) = (dn-1)*sin(mk*dn)/mk - (cos(mk) - cos(mk*dn))/mk^2;
%         Ek(ik,1) = Ek(ik,1)*nk^-0.5;
%     end
%     
% end
% 
% sum = 0.0;
% 
% for in2=2:NT
%     
%     in = INDN(in2);
%     
%     lambda = in*pi/dn;
%    
%     add = ALPHAN(in2,1)*(besseli(0,lambda*a)/(in*pi*besseli(1,lambda*a)) - ...
%         2.0*dn/(a*(in*pi)^2));
% 
%     sum = sum + add*(-1)^in;    %
%     
% end
% % 
% sum1 = 0.0;
% akrk_das(1:n_alp,1) = 0.0;
% 
% for ik=1:n_alp
%     if(ik==1)
%         mk = alpha_0;
%         
%         arg = alpha_0*dn;
%         fac1 = 0.375*a^2 - 1.5*dn^2 + dn;
%         fac2 = cosh(alpha_0) - sinh(arg)/arg - 2.0*cosh(arg);
%         A_star = fac1*sinh(arg)/arg - fac2/alpha_0^2;
%         A_star = A_star*N0^-0.5;
%     else
%         mk = alpha(ik-1,3);
%         
%         arg = mk*dn;
%         fac1 = 0.375*a^2 - 1.5*dn^2 + dn;
%         fac2 = 2.0*cos(arg) - cos(mk) - sin(arg)/arg;
%         A_star = fac1*sin(arg)/arg - fac2/mk^2;
%         A_star = A_star*Nk^-0.5;
%         
%     end
%     
% 	option = 1;
% 	RK = COMPRY(ik,a,alpha_0,alpha,option,mop);
%     
% 	option = 2;
% 	RK_das = COMPRY(ik,a,alpha_0,alpha,option,mop);
%     
%     sum3 =0.0;
% 	for in2=1:NT
%         
%         in = INDN(in2);
%         CNK = COMPCY(dn,ik,in,alpha_0,N0,alpha);
%         
%         sum3 = sum3 + SVEC(in2,1)*CNK*ALPHAN(in2,1);
%         
%     end
%         
% 
%         
% 
%     
%     Ak = (sum3 +  A_star)/(mk*RK_das);
%     
%     akrk_das(ik,1) = Ak*RK_das;
%     
%     sum1 = sum1 + Ak*RK*Ek(ik,1);
%     
% end
% 
% sum2 = 0.125*(a*dn - (a^3)/(6.0*dn) - ALPHAN(1,1));
% 
% FZ = sum1*(1/a)^2 ; %sum2 - (dn/a)*sum + 
% 
% % AM = real(FZ)/a
% % DC = imag(FZ)/a
% 
% AM(im,1) = real(FZ)/a;
% DC(im,1) = imag(FZ)/a;

%% end of Yeung's formulation %%%


%% begin of NH Nokob's formulation %%%

% psi_n_star(1:NT,1) = 0.0;
% Ak(1:NT,1) = 0.0;
% zk_star(1:NT,1) = 0.0;
% % tao(1:NT,1) = 0.0;
% 
% 
% for ik=1:NT
%     
%     indi = INDN(ik);
%     
%     if(indi==0)
%         
% %         m0 = alpha_0;
% %         sum1 = sinh(m0) - dn*sinh(m0*dn);
% %         sum2 = cosh(m0) - cosh(m0*dn);
% %         
% %         z = sum1/m0 - sum2/m0^2;
% 
%         z = (cosh(m0*dn) + m0*sinh(m0) - cosh(m0) - dn*m0*sinh(m0*dn))/(sqrt(N0)*m0^2);
% %         z = z/m0^2;
% 
% %         z = z*N0^-05;
%         zk_star(ik,1) = z;
%         
%         psi_n_star(ik,1) = 0.25*a^3;
% %         tao(ik,1) = (sinh(m0)-sinh(m0*dn))/(m0*sqrt(N0));
%     else
%         
%         mk = alpha(indi,3);
%         Nk = alpha(indi,4);
%         
% %         sum1 = sin(mk) - dn*sin(mk*dn);
% %         sum2 = cos(mk*dn) - cos(mk);
% %         
% %         z = sum1/mk - sum2/mk^2;
% 
%         z = (cos(mk) + mk*sin(mk) - dn*mk*sin(mk*dn) - cos(mk*dn))/(sqrt(Nk)*mk^2);
% %         z = z/mk^2;
% 
% 
% %         z = z*Nk^-0.5;
%         zk_star(ik,1) = z;
%         
%         lambda = indi*pi/dn;
%         
%         psi_n_star(ik,1) = a*(lambda*a*besseli(0,lambda*a) - 2.0...
%             *besseli(1,lambda*a))/(besseli(1,lambda*a)*lambda^2);
% 
% %           psi_n_star(ik,1) = a^2*besseli(2,lambda*a)/(lambda*besseli(1,lambda*a));
% 
% %         tao(ik,1) = (sin(mk)-sin(mk*dn))/(mk*sqrt(Nk));
%         
%     end
%     
% end
% 
% for ii=1:NT
%     
%     indi = INDN(ii);
%     
%     sum = 0.0;
%     
%     
%     for jj=1:NT
%         
% 
%         sum = sum + ALPHAN(jj,1)*PSID(jj,2)*CNKM(jj,ii)/PSID(jj,1);
% 
%         
%         
%     end
%     
%     Ak(ii,1) = sum + A_k_star(ii,1);
%     
% end
% 
% sum1 = 0.0;
% 
% for ii=1:NT
%     indi = INDN(ii);
%     
%     sum1 = sum1 + psi_n_star(ii,1)*ALPHAN(ii,1)*(-1.0)^indi;
%     
% end
% 
% sum2 = 0.0;
% 
% for ik = 1:NT
%     
% 
%     sum2 = sum2 + Ak(ik,1)*RKD(ik,1)*zk_star(ik,1)/RKD(ik,2);
% %     sum2 = sum2 + Ak(ik,1)*(RK/RK_das)*zk_star(ik,1);
% 
% %     sum2 = sum2 + Ak(ik,1)*RKD(ik,1)*tao(ik,1)/RKD(ik,2);
%     
% end
% 
% FZ = (6.0*a*dn^2 - a^3)/(48.0*dn) - sum1/a^3 - sum2/a^2;
% 
% % AM = real(FZ)
% % DC = imag(FZ)
%         
% % FZ = -1.0*sum2/a^2;
% 
% FZ_S(1,1) =  real(FZ);       
% FZ_S(1,2) = imag(FZ);        
     



%%% end of MH Nokob's formulation %%%

%%% end for roll %%%


%% begin for sway %%

% A_k(1:n_alp,1) = 0.0;
% A_k_star(1:n_alp,1) = 0.0;
% 
% for ik=1:n_alp
%     
%     if(ik==1)
%        A_k_star(ik,1) = sinh(alpha_0) - sinh(alpha_0*dn);
%        A_k_star(ik,1) = (N0^-0.5)*A_k_star(ik,1)/alpha_0;
%        
%     else
%        mk = alpha(ik-1,3);Nk = alpha(ik-1,4);
%        A_k_star(ik,1) = sin(mk) - sin(mk*dn);
%        A_k_star(ik,1) = (Nk^-0.5)*A_k_star(ik,1)/mk;
%        
%     end
%     
% end
% 
% for ik=1:n_alp
%     
%     sum = 0.0;
%     
%     if(ik==1)
%         mk = alpha_0;
%     else
%         mk = alpha(ik-1,3);
%     end
%     
%     for in2 = 1:NT
%     
%         in = INDN(in2);
%         
%         CNK = COMPCY(dn,ik,in,alpha_0,N0,alpha);
%         
%         SY = SVEC(in2,1);
%         
%         sum = sum + SY*CNK*ALPHAN(in2,1);
%         
%     end
%     
%     option = 2;
%     Rk_das  = COMPRY(ik,a,alpha_0,alpha,option,mop);
%     A_k(ik,1) = (sum + A_k_star(ik,1))/(mk*Rk_das);
%     
% end
% 
% sum = 0.0;
% 
% for ik=1:n_alp
%     
%     option=1;
%     Rk = COMPRY(ik,a,alpha_0,alpha,option,mop);
%     
%     sum = sum + A_k(ik,1)*A_k_star(ik,1)*Rk;
%     
% end
% 
% FZ = -1.0*sum/(a*(1.0 - dn));

%%% uncomment following for pitch %%%

% AM = real(FZ)
% DC = imag(FZ)                 


%% end for sway %%


        

% AM(im,1) = real(FZ);
% DC(im,1) = imag(FZ);

% end
% % 
% subplot(1,2,1)
% plot(INK,AM,'k.');
% 
% subplot(1,2,2)
% 
% plot(INK,DC,'k.');
% 
% [MM,II] = max(AM(:,1));

% z1 = AM;

% 
% % END OF FORCE EVALUATION %

Ak(1:NT,1) = 0.0;


for ik=1:NT
    
    if(ik==1)
        
        A_star = (a/(2.0*dn))*(1.0/alpha_0)*N0^-0.5;
        A_star = -1.0*A_star*sinh(alpha_0*dn);
        mk = alpha_0;
        
    else
        
        A_star = (a/(2.0*dn))*(1.0/alpha(ik-1,3))*alpha(ik-1,4)^-0.5;
        A_star = -1.0*A_star*sin(alpha(ik-1,3)*dn);
        mk = alpha(ik-1,3);
    end
    
    sum = 0.0;
    
    for in2=1:NT
        
        in = INDN(in2);
        CNK = COMPCY(dn,ik,in,alpha_0,N0,alpha);
        
        sum = sum + SVEC(in2,1)*CNK*ALPHAN(in2,1);
        
    end
    
	option = 2;
	RK_das = COMPRY(ik,a,alpha_0,alpha,option,mop);
    
    Ak(ik,1) = (sum + A_star)/(mk*RK_das);
    
% 	Ak(ik,1) = (sum + A_star)/mk;
% 
%     
    if(ik==1)
        fac = besselh(0,1,mk*a); %
    else
        fac = besselk(0,mk*a);
    end
    
    Ak(ik,1) = Ak(ik,1)*fac;

    
end



M(:,1) = real(Ak(:,1));
M(:,2) = imag(Ak(:,1));
M(:,3) = real(ALPHAN(:,1));
M(:,4) = imag(ALPHAN(:,1));


EXTC = M(:,1) + 1i*M(:,2);
INTC = M(:,3) + 1i*M(:,4);


% 
% alpha_w(1:NT,2) = 0.0;
% 
% for in = 1:NT
%     
%     if(in==1)
%         alpha_w(in,1) = alpha_0;
%         arg = 2.0*alpha_0*d;
%         N0 = 0.5*(1.0+sinh(arg))/arg;
%         alpha_w(in,2) = N0;
% %         alpha_w(in,2) = N0;
%     else
%         alpha_w(in,1) = alpha(in-1,3);
%         arg = 2.0*alpha(in-1,3)*d;
%         NN = 0.5*(1.0+sin(arg))/arg;
%         alpha_w(in,2) = NN;
% %         alpha_w(in,2) = alpha(in-1,4);
%     end
%     
% end
% 
% % print -dpdf sway_test_1.pdf
% % print('FillPageFigure','-dpdf','-fillpage')
% 
% 
% dlmwrite('H:\newraphs\coeffsy.txt',M,'delimiter','\t', 'precision', '%.6f');
% dlmwrite('C:\Users\sdechowdhury\Documents\MATLAB\HEAF.txt',FZ,'delimiter','\t', 'precision', '%.6f');
% dlmwrite('H:\newraphs\SVEC.txt',SVECIN,'delimiter','\t', 'precision', '%.6f');
% 
% dlmwrite('disp_roots_r.txt',alpha_w,'delimiter','\t', 'precision', '%.6f');

%%% write coefficients for external regions in pitch for multiple scattering %%%%%

% AKM(1:NT,1) = 0.0;
% 
% for ik=1:NT
%     indi = INDN(ik);
%     if(indi==0)
%         fac = besselh(1,1,m0*a)/dbesselh(1,1,m0*a);
%     else
%         mk = alpha(indi,3);
%         fac = besselk(1,mk*a)/dbesselk(1,mk*a);
%         
%     end
%     
%     AKM(ik,1) = Ak(ik,1)*fac;
%     
% end
% 
% M(:,1) = real(AKM(:,1));
% M(:,2) = imag(AKM(:,1));
% M(:,3) = real(ALPHAN(:,1));
% M(:,4) = imag(ALPHAN(:,1));
% 
% dlmwrite('H:\newraphs\coeffsy_pitch.txt',M,'delimiter','\t', 'precision', '%.6f');
% dlmwrite('H:\newraphs\PSINSTR.txt',psi_n_star,'delimiter','\t', 'precision', '%.6f');
% dlmwrite('H:\newraphs\ZKSTR.txt',zk_star,'delimiter','\t', 'precision', '%.6f');
% dlmwrite('H:\newraphs\RKDS.txt',RKD,'delimiter','\t', 'precision', '%.6f');
% 
% dlmwrite('H:\newraphs\FORCE_SINGLE.txt',FZ_S,'delimiter','\t', 'precision', '%.6f');
