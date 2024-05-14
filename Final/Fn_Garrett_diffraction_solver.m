    %% solve unknown Fourier coefficients for different m_orders %%

    eMatrix = zeros(nEqs, nEqs);
    
    for ik = 1:nEqs
        for jk = 1:nEqs
            
            sum = 0.;
            for ij = 1:nEqs
                id = ij-1;
                
                if (id == 0)
                    gmn_inv = 0.;
                else
                gmn = Fn_gFun_Garrett(ik,id,radius, clearance, m0, alpha(ik,3),m_order,1);
                gmn_inv = 1.0/gmn;
                end
                
                sum = sum + (clearance/depth)*jacobiSymbols(ij)*Lmatrix(ij,jk)*Lmatrix(ij,ik)*gmn_inv;
            end
            eMatrix(ik,jk) = sum;
            
            if (ik==jk)
                    
                gmn = Fn_gFun_Garrett(ik,id,radius, clearance, m0, alpha(ik,3),m_order,2);
  
                eMatrix(ik,jk) = eMatrix(ik,jk) + 1.0/gmn;

            end
            
        end
        
    end
    
    bVector = zeros(nEqs,1);
    bVector(1,1) = -2*1i/(pi*besselh(m_order, 1, m0*radius)*depth*z0_das);
    
    xVector = eMatrix\bVector;