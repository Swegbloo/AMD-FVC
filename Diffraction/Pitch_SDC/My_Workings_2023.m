clc;
clear all;

g = 9.81;
a = 0.2;
dpt = 1.0;
d = 0.5;

nEqns = 4;


% define

freqScale = 0.01:0.01:3;
%freqScale = 0.01;
nFreqs = size(freqScale,2);

% first column for heave, second for pitch.

aMass = zeros(nFreqs,2);    
damping = zeros(nFreqs,2);
difTrq = zeros(nFreqs,4); 

for ik = 1:nFreqs

    freq = freqScale(ik)/a;
sigma = (freq)*tanh(freq*dpt);

% difForce = Fn_diffractionForce(a, dpt, d, sigma);
% difForces(ik,1) = real(difForce);   %heave
% difForces(ik,2) = imag(difForce);
% difForces(ik,3) = sqrt(difForces(ik,1)^2+difForces(ik,2)^2);

difTrq = Fn_diffractionForce(a, dpt, d, sigma);
difTrq(ik,1) = real(difTrq);   %pitch
difTrq(ik,2) = imag(difTrq);
difTrq(ik,3) = difTrq(ik,1)+difTrq(ik,2);

mroots = dispersion_free_surface_vMikeM(sigma,nEqns,dpt);
% [fz_cmplx_iso,dvec,cvec] = cylinderSolverHeave_vYeung(a,d,mroots);

%pitchMoment = solcyl_pitch(a,dpt,d,mroots);

% fz_cmplx_iso = fz_cmplx_iso *a^3/(a^2*d);

omega = sqrt(sigma*g);

% aMass(ik,1) = real(fz_cmplx_iso);
% damping(ik,1) = imag(fz_cmplx_iso);
%return;

% aMass(ik,2) = real(pitchMoment);
% damping(ik,2) = imag(pitchMoment);

fprintf('%s\n',['computed for frequency ', num2str(freqScale(ik)),'..']);



end

figure(1)
plot(freqScale, difTrq(:,3))

% figure(2)
% plot(freqScale, difForces(:,2))
