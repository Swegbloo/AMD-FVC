function[z] = COMPSY(in,a,dn)

%% for heave %%
a2 = 0.5*in*pi;

arg = in*pi*a/dn;
b2 = dbesseli(0,arg);
c2 = besseli(0,arg);

z = a2*b2/c2;

% if(isnan(z))
%     z = 0.0;
% end

%% end for heave %%



