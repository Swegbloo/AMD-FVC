clc;
clear all;
% addpath ../GrafInteraction_GitHubDesktop/Fns/;

g = 9.81;
% a = 0.2;
dpt = 1.0;
a=0.5;



d = 0.9;


nEqs = 5;
formulation = 'Garrett';


% define



freqScale = 0.01:0.01:8;



nFreqs = size(freqScale,2);

% first column for heave, second for pitch.
enrg_h = zeros(1,nFreqs);
enrg_d = zeros(1,nFreqs);
freqinv = zeros(1,nFreqs);
l = zeros(1,nFreqs);
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


 fz_cmplx_iso = fz_cmplx_iso *a^3/(a^2*d);




% omega = sqrt(sigma*g);


aMass(ik,1) = real(fz_cmplx_iso);
damping(ik,1) = imag(fz_cmplx_iso);



aMass(ik,3) = real(pitchMoment);
damping(ik,3) = imag(pitchMoment);

%% do added mass and damping for surge here and store as aMass(ik,2) and damping(ik,2) %% 
disp('Type 1 for Heave ONLY: ');
disp('Type 2 for Surge ONLY: ');
disp('Type 3 for Pitch ONLY: ');
disp('Type 4 for Heave, Surge & Pitch UNCOUPLED: ');
"\n";
ip = input();
fprintf('%s\n',['computed for frequency ', num2str(freqScale(ik)),'..']);
if ip == 1
    enrg_h(1,ik) = energy_harnessed(0,0,0,0,damping(ik,2)*1025*pi*(a*dpt)^2*(1-d)*dpt*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,6)*B*A); % correction needed with amplitude dimensionalization
    %damping(ik,1)*1025*pi*(a*dpt)^3*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,3)*B*A,
    enrg_d(1,ik) = energy_delivered(dpt,freq,A);
    l(1,ik) = enrg_h(1,ik)/enrg_d(1,ik); %capture width
    freqinv(ik) = pi/freqScale(ik);
elseif ip == 2
    enrg_h(1,ik) = energy_harnessed(0,0,0,0,damping(ik,2)*1025*pi*(a*dpt)^2*(1-d)*dpt*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,6)*B*A); % correction needed with amplitude dimensionalization
    %damping(ik,1)*1025*pi*(a*dpt)^3*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,3)*B*A,
    enrg_d(1,ik) = energy_delivered(dpt,freq,A);
    l(1,ik) = enrg_h(1,ik)/enrg_d(1,ik); %capture width
    freqinv(ik) = pi/freqScale(ik);
elseif ip == 3
    enrg_h(1,ik) = energy_harnessed(0,0,0,0,damping(ik,2)*1025*pi*(a*dpt)^2*(1-d)*dpt*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,6)*B*A); % correction needed with amplitude dimensionalization
    %damping(ik,1)*1025*pi*(a*dpt)^3*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,3)*B*A,
    enrg_d(1,ik) = energy_delivered(dpt,freq,A);
    l(1,ik) = enrg_h(1,ik)/enrg_d(1,ik); %capture width
    freqinv(ik) = pi/freqScale(ik);
elseif ip == 4
    enrg_h(1,ik) = energy_harnessed(0,0,0,0,damping(ik,2)*1025*pi*(a*dpt)^2*(1-d)*dpt*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,6)*B*A); % correction needed with amplitude dimensionalization
    %damping(ik,1)*1025*pi*(a*dpt)^3*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,3)*B*A,
    enrg_d(1,ik) = energy_delivered(dpt,freq,A);
    l(1,ik) = enrg_h(1,ik)/enrg_d(1,ik); %capture width
    freqinv(ik) = pi/freqScale(ik);
end
end

% return

% plot(freqScale, abs(diffMoments(:,3)),'.');
% return

% outFileName = 'test_solutions/heaveCoefficients.txt';
% dlmwrite(outFileName,[freqScale' difForces(:,3) aMass(:,1)/a damping(:,1)/a],'delimiter','\t','precision',4)


figure(1)
plot(freqScale, aMass(:,3)/a,'x')


% refData = csvread('test_solutions\Garrett_torques_dbya_0.75_hbya_0.25.csv');
% 
% plot(freqScale, difForces(:,9))
% hold on;
% plot(refData(:,1), refData(:,2),'x');   


% outFileName = 'test_solutions/heaveDiffraction_Garrett.txt';
% dlmwrite(outFileName,[freqScale' difForces],'delimiter','\t','precision',4)

figure(2)
plot(freqScale, damping(:,3)/a,'x')
% figure(2)
% plot(freqScale, difForces(:,2))


figure(3)
plot(freqScale, difForces(:,3),'.');

