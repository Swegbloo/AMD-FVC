clc;
%clear all;
addpath ../GrafInteraction_GitHubDesktop/Fns/;

g = 9.81;
a = 0.5;
dpt = 1.0;
d = 0.1;

nEqns = 2;


% define

freqScale = 0.4:0.01:4;
% freqScale = 0.1;


nFreqs = size(freqScale,2);

% first column for heave, second for pitch.


aMass = zeros(nFreqs,2);    
damping = zeros(nFreqs,2);
difForces = zeros(nFreqs,3); 
diffMoments = zeros(nFreqs,2);
enrg_h = zeros(1,nFreqs);
enrg_d = zeros(1,nFreqs);
freqinv = zeros(1,nFreqs);
l = zeros(1,nFreqs);
p = 1025;
g = 9.81;
A = 1;
B = pi*a^2*dpt*(1-d)*p*g;
for ik = 1:nFreqs

freq = freqScale(ik)/a;
sigma = (freq)*tanh(freq*dpt);

% difForce = Fn_diffractionForce(a, dpt, d, sigma);

%pitchDiffMoment = Fn_diffractionTorque(a,d,sigma);
% return
%diffMoments(ik,1) = abs(pitchDiffMoment);
% return
 difForce = Fn_diffractionForce(a, dpt, d, sigma);


 difForces(ik,1) = real(difForce);
 difForces(ik,2) = imag(difForce);
 difForces(ik,3) = sqrt(difForces(ik,1)^2+difForces(ik,2)^2);
 mroots = dispersion_free_surface_vMikeM(sigma,nEqns,dpt);
% return

[fz_cmplx_iso,dvec,cvec] = cylinderSolverHeave_vYeung(a,d,mroots);

% pitchMoment = solcyl_pitch(a,dpt,d,mroots);

% fz_cmplx_iso = fz_cmplx_iso *a^3/(a^2*d);




% omega = sqrt(sigma*g);

aMass(ik,1) = real(fz_cmplx_iso);
damping(ik,1) = imag(fz_cmplx_iso);



% aMass(ik,2) = real(pitchMoment);
% damping(ik,2) = imag(pitchMoment);


fprintf('%s\n',['computed for frequency ', num2str(freqScale(ik)),'..']);
if ik == 1
    m(1) = -1i*mroots(1);
end
enrg_h(1,ik) = energy_harnessed(damping(ik,1)*1025*pi*(a*dpt)^3*sqrt(g/dpt*freq*tanh(freq)),difForces(ik,3)*B*A,0,0); % correction needed with amplitude dimensionalization
enrg_d(1,ik) = energy_delivered(dpt,freq,A);
l(1,ik) = enrg_h(1,ik)/enrg_d(1,ik); %capture width
freqinv(ik) = pi/freqScale(ik);
end
%disp(diffMoments);
%plot(freqScale, real(difForce(:,1)),'.');
%figure(1)
hold on;
plot(freqScale*dpt, l/(2*a*dpt))
%plot(datapoints2_paper(:,1),datapoints2_paper(:,2));
plot(capture_width_pap(:,1),capture_width_pap(:,2));
% figure(1)
% plot(freqScale, damping(:,1)/a)
%plot(freqScale, abs(difForces(:,3)));

% figure(2)
% plot(freqScale, aMass(:,1)/a)
% figure(2)
% plot(freqScale, difForces(:,2))