function[z] = Fn_gFun_Garrett(row_id, id,radius, clearance, alpha_0, alp, m_order,option)

if (option == 1)
    arg = id*pi*radius/clearance;
    z = besseli(m_order,arg)/(arg*dbesseli(m_order,arg));
elseif (option == 2)
    if (row_id == 1)
        arg = alpha_0*radius;
        z = -besselh(m_order,arg)/(arg*dbesselh(m_order,1,arg));
    else
    arg = alp*radius;
    z = -besselk(m_order,arg)/(arg*dbesselk(m_order,arg));
    end
end

return
    