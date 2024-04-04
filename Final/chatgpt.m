clc;
clear all;
addpath ../GrafInteraction_GitHubDesktop/Fns/;

g = 9.81;
a = 0.2;
dpt = 1.0;
d = 0.95;

nEqns = 2;

% Define frequency scale
freqScale = 0.01:0.1:4;
nFreqs = size(freqScale, 2);

% Arrays to store results
pitchMoments = zeros(nFreqs, 1);

for ik = 1:nFreqs
    freq = freqScale(ik) / a;
    sigma = (freq) * tanh(freq * dpt);

    pitchMoment = Fn_pitchingMoment(a, dpt, d, sigma);

    pitchMoments(ik) = pitchMoment;

    fprintf('%s\n', ['computed for frequency ', num2str(freqScale(ik)), '..']);
end

figure;
plot(freqScale, abs(pitchMoments));
xlabel('Frequency (Hz)');
ylabel('Pitching Moment');
title('Pitching Moments vs Frequency');

function [pitchMoment] = Fn_pitchingMoment(radius, depth, clearance, sigma)
    nEqs = 5;

    mroots = dispersion_free_surface_vMikeM(sigma, nEqs, depth);
    m0 = -1i * mroots(1);

    alpha = zeros(nEqs, 4);
    alpha(:, 3) = mroots(2:end).';

    A0_star = zeros(nEqs, 1);
    vEigNs = zeros(nEqs, 1);

    for ik = 1:nEqs
        if ik == 1
            vEigNs(ik) = 0.5 + sinh(2 * m0) / (4 * m0);
        else
            vEigNs(ik) = 0.5 + sin(2 * mroots(ik)) / (4 * mroots(ik));
        end
    end

    z0_das = m0 * sinh(m0) / sqrt(vEigNs(1));
    N0 = vEigNs(1);
    alpha(:, 4) = vEigNs(2:end).';

    for ik = 1:nEqs
        lambda = (ik - 1) * pi / clearance;
        mj = mroots(ik);

        rj = COMPRY(ik, radius, m0, alpha, 1, 0);
        rj_das = COMPRY(ik, radius, m0, alpha, 2, 0);

        if ik == 1
            rj_das = m0 * rj_das;
            A0_star(ik) = -1.0 * rj * 1i * m0 * dbesselj(0, m0 * radius) * (1 + 0.5 * sinh(2 * m0) / m0) / ...
                (2 * vEigNs(1) * rj_das * z0_das);
        else
            rj_das = mj * rj_das;
            A0_star(ik) = -1.0 * rj * 1i * m0 * dbesselj(0, m0 * radius) * (mj * cosh(m0) * sin(mj) + m0 * sinh(m0) * cos(mj)) / ...
                (rj_das * z0_das * sqrt(vEigNs(1) * vEigNs(ik)) * (m0^2 + mj^2));
        end
    end

    pitchMoment = -2.0 * pi * A0_star(1);
end

function [z] = COMPRY(dn, ik, in, alpha_0, N0, alpha)
    if (ik == 1)
        arg = alpha_0 * dn;
        z = 2.0 * sinh(arg) / arg;
        z = z * (-1)^in;

        z2 = 1.0 + (in * pi / arg)^2;
        z2 = z2 * N0^0.5;
        z = z / z2;
    else
        mk = alpha(ik - 1, 3);
        nk = alpha(ik - 1, 4);

        arg1 = mk * dn - in * pi;
        arg2 = mk * dn + in * pi;

        z = 2.0 * sin(arg1) / arg1;
        z = z * mk * dn / arg2;
        z = z * nk^-0.5;
    end
end
