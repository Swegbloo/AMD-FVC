clc;
clear all;
% addpath ../GrafInteraction_GitHubDesktop/Fns/;
addpath("C:\Users\Swegbloo\Documents\GitHub\AMD-FVC\Radiation\Sway");
g = 9.81;
% a = 0.2;
dpt = 1;
a=0.5;
d = 0.3333;
load("num_heave_am.mat");
load("num_heave_dp.mat");
load("num_heave_df.mat");
load("AM_Heave_d0.1a0.5.mat");
load("CW_greens.mat")
load("num_heave_am_2x.mat");
load("num_heave_df_2x.mat");
load("num_heave_dp_2x.mat");
load("num_heave_am_4x.mat");
load("num_heave_dp_4x.mat");
load("num_heave_df_4x.mat");
load("num_heave_df_4x_S2x.mat");
load("num_heave_df_4x_R2x.mat");
load("num_heave_am_8x.mat");
load("num_heave_dp_8x.mat");
load("num_heave_df_8x.mat");
load("num_heave_am_16x.mat");
load("num_heave_dp_16x.mat");
load("num_heave_df_16x.mat");
%load("num_surge_am.mat");
load("num_surge_am_1x.mat");
load("num_surge_am_2x.mat");
load("num_surge_am_4x.mat");
load("num_surge_am_8x.mat");
load("num_surge_dp.mat");
load("num_surge_dp_2x.mat");
load("num_surge_dp_4x.mat");
load("num_surge_dp_8x.mat");
load("num_surge_df.mat");
load("num_surge_df_2x.mat");
load("num_surge_df_4x.mat");
load("num_surge_df_8x.mat");
load("num_pitch_am.mat");
load("num_pitch_am_2x.mat");
load("num_pitch_am_4x.mat");
load("num_pitch_am_8x.mat");
load("num_pitch_dp.mat");
load("num_pitch_dp_2x.mat");
load("num_pitch_dp_4x.mat");
load("num_pitch_dp_8x.mat");
load("num_pitch_df.mat");
load("num_pitch_df_2x.mat");
load("num_pitch_df_4x.mat");
load("num_pitch_df_8x.mat");


AM = num_heave_am(:,2);
num_freqS = num_heave_am(:,1);
DP = num_heave_dp(:,2);
DF = num_heave_df(:,2);
HAMYeung = Data001;
AM2 = num_heave_am_2x(:,2);
DF2 = num_heave_df_2x(:,2);
DP2 = num_heave_dp_2x(:,2);
AM4 = num_heave_am_4x(:,2);
DF4 = num_heave_df_4x(:,2);
DP4 = num_heave_dp_4x(:,2);
DF4S2 = num_heave_df_4x_S2x(:,2);
DF4R2 = num_heave_df_4x_R2x(:,2);
AM8 = num_heave_am_8x(:,2);
DP8 = num_heave_dp_8x(:,2);
DF8 = num_heave_df_8x(:,2);
AM16 = num_heave_am_16x(:,2);
DP16 = num_heave_dp_16x(:,2);
DF16 = num_heave_df_16x(:,2);
%AM_S = num_surge_am(:,2);
AM_S = num_surge_am_1x(:,2);
AM_S_2 = num_surge_am_2x(:,2);
AM_S_4 = num_surge_am_4x(:,2);
AM_S_8 = num_surge_am_8x(:,2);
DP_S = num_surge_dp(:,2);
DP_S_2 = num_surge_dp_2x(:,2);
DP_S_4 = num_surge_dp_4x(:,2);
DP_S_8 = num_surge_dp_8x(:,2);
DF_S = num_surge_df(:,2);
DF_S_2 = num_surge_df_2x(:,2);
DF_S_4 = num_surge_df_4x(:,2);
DF_S_8 = num_surge_df_8x(:,2);
AM_P = num_pitch_am(:,2);
AM_P_2 = num_pitch_am_2x(:,2);
AM_P_4 = num_pitch_am_4x(:,2);
AM_P_8 = num_pitch_am_8x(:,2);
DP_P = num_pitch_dp(:,2);
DP_P_2 = num_pitch_dp_2x(:,2);
DP_P_4 = num_pitch_dp_4x(:,2);
DP_P_8 = num_pitch_dp_8x(:,2);
DF_P = num_pitch_df(:,2);
DF_P_2 = num_pitch_df_2x(:,2);
DF_P_4 = num_pitch_df_4x(:,2);
DF_P_8 = num_pitch_df_8x(:,2);

nEqs = 20;
formulation = 'Garrett';


% define
for i = 1:50
    DP(i) = DP(i)/num_heave_dp(i,3);
    DP2(i) = DP2(i)/num_heave_dp(i,3);
    DP4(i) = DP4(i)/num_heave_dp(i,3);
    DP8(i) = DP8(i)/num_heave_dp(i,3);
    DP16(i) = DP16(i)/num_heave_dp(i,3);

    DP_S(i) = DP_S(i)/num_heave_dp(i,3);
    DP_S_2(i) = DP_S_2(i)/num_heave_dp(i,3);
    DP_S_4(i) = DP_S_4(i)/num_heave_dp(i,3);
    DP_S_8(i) = DP_S_8(i)/num_heave_dp(i,3);

    DP_P(i) = DP_P(i)/num_heave_dp(i,3);
    DP_P_2(i) = DP_P_2(i)/num_heave_dp(i,3);
    DP_P_4(i) = DP_P_4(i)/num_heave_dp(i,3);
    DP_P_8(i) = DP_P_8(i)/num_heave_dp(i,3);
end


freqScale = 0.4:0.01:4;


nFreqs = size(freqScale,2);

% first column for heave, second for pitch.
enrg_h = zeros(1,nFreqs);
enrg_d = zeros(1,nFreqs);
enrg_h_num = zeros(1,50);
enrg_d_num = zeros(1,50);
l_num_heave = zeros(1,50);
enrg_h_num_2x = zeros(1,50);
enrg_d_num_2x = zeros(1,50);
l_num_heave_2x = zeros(1,50);
enrg_h_num_4x = zeros(1,50);
enrg_h_num_8x = zeros(1,50);
enrg_d_num_8x = zeros(1,50);
l_num_heave_8x = zeros(1,50);
enrg_d_num_4x = zeros(1,50);
l_num_heave_4x = zeros(1,50);
freqinv = zeros(1,nFreqs);
l_a_heave = zeros(1,nFreqs);
p = 1025;
A = 1;
B = pi*a^2*dpt*(1-d)*p*g;



aMass = zeros(nFreqs,3);   % 1st col: heave, 2nd col: surge, 3rd col: pitch 
damping = zeros(nFreqs,3); % same as damping




difForces = zeros(nFreqs,4); 
diffMoment = zeros(nFreqs,3);

% tempSums = zeros(nFreqs,1);


[aMass(:,2),damping(:,2)] = sway(a,d,nEqs,freqScale);


for ik = 1:nFreqs

    freq = freqScale(ik)/a;
sigma = (freq)*tanh(freq*dpt);
lambda(ik,1) = 2*pi*dpt/freq;

% pitchDiffMoment = Fn_diffractionTorque(nEqs, a,d,sigma);
% return
% diffMoments(ik,:) = pitchDiffMoment(:).';
% return

[difForce, difMoment] = Fn_diffractionForce(a, dpt, d, sigma,formulation,{'heave',...
    'surge', 'pitch'});
% % % return
% % 
% % heave diffraction force components
% 
difForces(ik,1) = real(difForce(1));
difForces(ik,2) = imag(difForce(1));
difForces(ik,3) = abs(difForce(1));
% 
% % surge diffraction force components
difForces(ik,4) = real(difForce(2));
difForces(ik,5) = imag(difForce(2));
difForces(ik,6) = abs(difForce(2));
% 
% % pitch moment components
difForces(ik,7) = real(difMoment);
difForces(ik,8) = imag(difMoment);
difForces(ik,9) = abs(difMoment);



% tempSums(ik,1) = tempSum;

 mroots = dispersion_free_surface_vMikeM(sigma,nEqs,dpt);
% return

 [fz_cmplx_iso,dvec,cvec] = cylinderSolverHeave_vYeung(a,d,mroots);

% mroots(1) = -1i*mroots(1);
disp(a);
disp(d);
disp(nEqs);
disp(mroots);


pitchMoment = solcyl_pitch(a,dpt,d,mroots);


 %fz_cmplx_iso = fz_cmplx_iso *a^3/(a^2*d);




% omega = sqrt(sigma*g);


aMass(ik,1) = real(fz_cmplx_iso);
damping(ik,1) = imag(fz_cmplx_iso);



aMass(ik,3) = real(pitchMoment);
damping(ik,3) = imag(pitchMoment);

%% do added mass and damping for surge here and store as aMass(ik,2) and damping(ik,2) %% 
% disp('Type 1 for Heave ONLY: ');
% disp('Type 2 for Surge ONLY: ');
% disp('Type 3 for Pitch ONLY: ');
% disp('Type 4 for Heave, Surge & Pitch UNCOUPLED: ');
% "\n";
% ip = input();
fprintf('%s\n',['computed for frequency ', num2str(freqScale(ik)),'..']);
% if ip == 1
    
    enrg_h(1,ik) = energy_harnessed(damping(ik,1)*1025*pi*(a*dpt)^3*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,3)*B*A,0,0,0,0); % correction needed with amplitude dimensionalization
    %damping(ik,1)*1025*pi*(a*dpt)^3*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,3)*B*A,
    enrg_d(1,ik) = energy_delivered(dpt,freq,A);
    l_a_heave(1,ik) = enrg_h(1,ik)/enrg_d(1,ik); %capture width analytical
    freqinv(ik) = pi/freqScale(ik);
% elseif ip == 2
%     enrg_h(1,ik) = energy_harnessed(0,0,0,0,damping(ik,2)*1025*pi*(a*dpt)^2*(1-d)*dpt*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,6)*B*A); % correction needed with amplitude dimensionalization
%     %damping(ik,1)*1025*pi*(a*dpt)^3*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,3)*B*A,
%     enrg_d(1,ik) = energy_delivered(dpt,freq,A);
%     l(1,ik) = enrg_h(1,ik)/enrg_d(1,ik); %capture width
%     freqinv(ik) = pi/freqScale(ik);
% elseif ip == 3
%     enrg_h(1,ik) = energy_harnessed(0,0,0,0,damping(ik,2)*1025*pi*(a*dpt)^2*(1-d)*dpt*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,6)*B*A); % correction needed with amplitude dimensionalization
%     %damping(ik,1)*1025*pi*(a*dpt)^3*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,3)*B*A,
%     enrg_d(1,ik) = energy_delivered(dpt,freq,A);
%     l(1,ik) = enrg_h(1,ik)/enrg_d(1,ik); %capture width
%     freqinv(ik) = pi/freqScale(ik);
% elseif ip == 4
%     enrg_h(1,ik) = energy_harnessed(0,0,0,0,damping(ik,2)*1025*pi*(a*dpt)^2*(1-d)*dpt*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,6)*B*A); % correction needed with amplitude dimensionalization
%     %damping(ik,1)*1025*pi*(a*dpt)^3*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,3)*B*A,
%     enrg_d(1,ik) = energy_delivered(dpt,freq,A);
%     l(1,ik) = enrg_h(1,ik)/enrg_d(1,ik); %capture width
%     freqinv(ik) = pi/freqScale(ik);
% end
end
for ik = 1:50
    enrg_h_num(1,ik) = energy_harnessed(DP(ik),DF(ik),0,0,0,0);
    enrg_d_num(1,ik) = energy_delivered(dpt,num_freqS(ik)/a,A);
    l_num_heave(1,ik) = enrg_h_num(1,ik)/enrg_d_num(1,ik); %capture width numerical

    enrg_h_num_2x(1,ik) = energy_harnessed(DP2(ik),DF2(ik),0,0,0,0);
    enrg_d_num_2x(1,ik) = energy_delivered(dpt,num_freqS(ik)/a,A);
    l_num_heave_2x(1,ik) = enrg_h_num_2x(1,ik)/enrg_d_num_2x(1,ik); %capture width numerical 2x

    enrg_h_num_4x(1,ik) = energy_harnessed(DP4(ik),DF4(ik),0,0,0,0);
    enrg_d_num_4x(1,ik) = energy_delivered(dpt,num_freqS(ik)/a,A);
    l_num_heave_4x(1,ik) = enrg_h_num_4x(1,ik)/enrg_d_num_4x(1,ik); %capture width numerical 4x

    enrg_h_num_8x(1,ik) = energy_harnessed(DP8(ik),DF8(ik),0,0,0,0);
    enrg_d_num_8x(1,ik) = energy_delivered(dpt,num_freqS(ik)/a,A);
    l_num_heave_8x(1,ik) = enrg_h_num_8x(1,ik)/enrg_d_num_8x(1,ik); %capture width numerical 8x
end

% return

% plot(freqScale, abs(diffMoments(:,3)),'.');
% return

% outFileName = 'test_solutions/heaveCoefficients.txt';
% dlmwrite(outFileName,[freqScale' difForces(:,3) aMass(:,1)/a damping(:,1)/a],'delimiter','\t','precision',4)

figure(1)
hold on;
% plot(freqScale, aMass(:,1)/a,'x')
% plot(num_freqS,AM/(1025*pi*a^4),'.');
% plot(HAMYeung(:,1),HAMYeung(:,2));
% plot(num_freqS,AM2/(1025*pi*a^4),'+');
% plot(num_freqS,AM4/(1025*pi*a^4),'-');
% plot(num_freqS,AM8/(1025*pi*a^4),'X');
% plot(num_freqS,AM16/(1025*pi*a^4),'-');
% title('Heave Added Mass (a = 0.5 & d = 0.1)');
% ylabel('µ_3_3/a');
% xlabel('m_0.a');
% legend('Yeung 1981','NEMOH 231 Panels','Yeung Screenshot','NEMOH 483 Panels','NEMOH 987 Panels','NEMOH 1995 Panels','NEMOH 4032 Panels');

% plot(freqScale, aMass(:,2),'x');
% plot(num_freqS,AM_S/(1025*pi*a^2*(1-d)),'.');
% plot(num_freqS,AM_S_2/(1025*pi*a^2*(1-d)),'x');
% plot(num_freqS,AM_S_4/(1025*pi*a^2*(1-d)),'+');
% plot(num_freqS,AM_S_8/(1025*pi*a^2*(1-d)),'-');
% title('Surge Added Mass (a = 0.5 & d = 0.1)');
% legend('Yeung 1981','NEMOH 231 Panels','NEMOH 483 Panels','NEMOH 987 Panels','NEMOH 1995 Panels');

plot(freqScale, aMass(:,3)/a,'x');
plot(num_freqS,AM_P/(1025*pi*a^4*dpt^2),'.');
plot(num_freqS,AM_P_2/(1025*pi*a^4*dpt^2),'x');
plot(num_freqS,AM_P_4/(1025*pi*a^4*dpt^2),'+');
plot(num_freqS,AM_P_8/(1025*pi*a^4*dpt^2),'-');
title('Pitch Added Mass (a = 0.5 & d = 0.1)');
legend('Yeung 1981','NEMOH 231 Panels','NEMOH 483 Panels','NEMOH 987 Panels','NEMOH 1995 Panels');


hold off;

% refData = csvread('test_solutions\Garrett_torques_dbya_0.75_hbya_0.25.csv');
% 
% plot(freqScale, difForces(:,9))
% hold on;
% plot(refData(:,1), refData(:,2),'x');   


% outFileName = 'test_solutions/heaveDiffraction_Garrett.txt';
% dlmwrite(outFileName,[freqScale' difForces],'delimiter','\t','precision',4)

figure(2)
hold on;
% plot(freqScale, damping(:,1)/a,'x');
% plot(num_freqS,DP/(1025*pi*a^4),'.');
% plot(num_freqS,DP2/(1025*pi*a^4),'+');
% plot(num_freqS,DP4/(1025*pi*a^4),'-');
% plot(num_freqS,DP8/(1025*pi*a^4),'X');
% plot(num_freqS,DP16/(1025*pi*a^4),'-');
% title('Heave Damping (a = 0.5 & d = 0.1)');
% ylabel('λ_3_3/a');
% xlabel('m_0.a');
% legend('Yeung 1981','NEMOH 231 Panels','NEMOH 483 Panels','NEMOH 987 Panels','NEMOH 1995 Panels','NEMOH 4032 Panels');

% plot(freqScale, damping(:,2),'x');
% plot(num_freqS,DP_S/(1025*pi*a^2*(1-d)),'.');
% plot(num_freqS,DP_S_2/(1025*pi*a^2*(1-d)),'+');
% plot(num_freqS,DP_S_4/(1025*pi*a^2*(1-d)),'-');
% plot(num_freqS,DP_S_8/(1025*pi*a^2*(1-d)),'X');
% plot(num_freqS,DP16/(1025*pi*a^4),'-');
% title('Surge Damping (a = 0.5 & d = 0.1)');
% ylabel('λ_3_3/a');
% xlabel('m_0.a');
% legend('Yeung 1981','NEMOH 231 Panels','NEMOH 483 Panels','NEMOH 987 Panels','NEMOH 1995 Panels');

plot(freqScale, damping(:,3)/a,'x');
plot(num_freqS,DP_P/(1025*pi*a^4*dpt^2),'.');
plot(num_freqS,DP_P_2/(1025*pi*a^4*dpt^2),'+');
plot(num_freqS,DP_P_4/(1025*pi*a^4*dpt^2),'-');
plot(num_freqS,DP_P_8/(1025*pi*a^4*dpt^2),'X');
% plot(num_freqS,DP16/(1025*pi*a^4),'-');
title('Pitch Damping (a = 0.5 & d = 0.1)');
ylabel('λ_3_3/a');
xlabel('m_0.a');
legend('Yeung 1981','NEMOH 231 Panels','NEMOH 483 Panels','NEMOH 987 Panels','NEMOH 1995 Panels');

hold off;
% figure(2)
% plot(freqScale, difForces(:,2))


figure(3)
hold on;
% plot(freqScale, difForces(:,3),'-');
% plot(num_freqS,DF/(pi*a^2*dpt*(1-d)*1025*9.81)),'.';
% plot(num_freqS,DF2/(pi*a^2*dpt*(1-d)*1025*9.81),'+');
% plot(num_freqS,DF4/(pi*a^2*dpt*(1-d)*1025*9.81),'.');
% plot(num_freqS,DF4S2/(pi*a^2*dpt*(1-d)*1025*9.81),'x');
% plot(num_freqS,DF4R2/(pi*a^2*dpt*(1-d)*1025*9.81),'+');
% plot(num_freqS,DF8/(pi*a^2*dpt*(1-d)*1025*9.81),'-');
% plot(num_freqS,DF16,'-');

% title('Heave Diffraction Force (a = 0.5 & d = 0.1)');
% ylabel('F_3_3/a');
% xlabel('m_0.a');
% legend('Yeung 1981','NEMOH 987','NEMOH 987 + 2xS Panels','NEMOH 987 + 2xR Panels','NEMOH 4023 Panels');

% plot(freqScale, difForces(:,6),'-');
% plot(num_freqS,DF_S/(pi*a^2*dpt*(1-d)*1025*9.81)),'.';
% plot(num_freqS,DF_S_2/(pi*a^2*dpt*(1-d)*1025*9.81),'+');
% plot(num_freqS,DF_S_4/(pi*a^2*dpt*(1-d)*1025*9.81),'.');
% plot(num_freqS,DF_S_8/(pi*a^2*dpt*(1-d)*1025*9.81),'.');
%plot(num_freqS,DF4S2/(pi*a^2*dpt*(1-d)*1025*9.81),'x');
%plot(num_freqS,DF4R2/(pi*a^2*dpt*(1-d)*1025*9.81),'+');
%plot(num_freqS,DF8/(pi*a^2*dpt*(1-d)*1025*9.81),'-');
%plot(num_freqS,DF16,'-');

% title('Surge Diffraction Force (a = 0.5 & d = 0.1)');
% ylabel('F_3_3/a');
% xlabel('m_0.a');
% legend('Yeung 1981','NEMOH 231 Panels','NEMOH 483 Panels','NEMOH 987 Panels','NEMOH 1995 Panels');

% plot(freqScale, difForces(:,6),'-');
% plot(num_freqS,DF_S/(pi*a^2*dpt*(1-d)*1025*9.81)),'.';
% plot(num_freqS,DF_S_2/(pi*a^2*dpt*(1-d)*1025*9.81),'+');
% plot(num_freqS,DF_S_4/(pi*a^2*dpt*(1-d)*1025*9.81),'.');
% plot(num_freqS,DF_S_8/(pi*a^2*dpt*(1-d)*1025*9.81),'.');
%plot(num_freqS,DF4S2/(pi*a^2*dpt*(1-d)*1025*9.81),'x');
%plot(num_freqS,DF4R2/(pi*a^2*dpt*(1-d)*1025*9.81),'+');
%plot(num_freqS,DF16,'-');

% title('Surge Diffraction Force (a = 0.5 & d = 0.1)');
% ylabel('F_3_3/a');
% xlabel('m_0.a');
% legend('Yeung 1981','NEMOH 231 Panels','NEMOH 483 Panels','NEMOH 987 Panels','NEMOH 1995 Panels');

plot(freqScale, difForces(:,9),'-');
plot(num_freqS,DF_P/(pi*a^3*dpt*(1-d)*1025*9.81)),'.';
plot(num_freqS,DF_P_2/(pi*a^3*dpt*(1-d)*1025*9.81),'+');
plot(num_freqS,DF_P_4/(pi*a^3*dpt*(1-d)*1025*9.81),'.');
plot(num_freqS,DF_P_8/(pi*a^3*dpt*(1-d)*1025*9.81),'.');
%plot(num_freqS,DF4S2/(pi*a^2*dpt*(1-d)*1025*9.81),'x');
%plot(num_freqS,DF4R2/(pi*a^2*dpt*(1-d)*1025*9.81),'+');
%plot(num_freqS,DF8/(pi*a^2*dpt*(1-d)*1025*9.81),'-');
%plot(num_freqS,DF16,'-');

title('Pitch Diffraction Force (a = 0.5 & d = 0.1)');
ylabel('F_3_3/a');
xlabel('m_0.a');
legend('Yeung 1981','NEMOH 231 Panels','NEMOH 483 Panels','NEMOH 987 Panels','NEMOH 1995 Panels');
hold off;

figure(4)
hold on;
plot(freqScale,l_a_heave/(2*a),'.');
plot(CW_greens(:,1),CW_greens(:,2),'+');
plot(num_freqS,l_num_heave/(2*a),'-');
plot(num_freqS,l_num_heave_2x/(2*a),'-');
plot(num_freqS,l_num_heave_4x/(2*a),'x');
plot(num_freqS,l_num_heave_8x/(2*a),'-');

title('Heave Capture Width (a = 0.5 & d = 0.1)');
ylabel('CW/(2.a)');
xlabel('m_0.a');
legend('Yeung 1981','Greens Theorem','NEMOH 231 Panels','NEMOH 483 Panels','NEMOH 987 Panels','NEMOH 1995 Panels');
hold off;

