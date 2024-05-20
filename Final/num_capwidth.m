a=0.5;
clc;
B3 = zeros(50);
DF3 = zeros(50);
B5 = zeros(50);
DF5 = zeros(50);
B1 = zeros(50);
DF1 = zeros(50);

load("num_heave_dp.mat");
B3 = num_heave_dp;
load("num_heave_df.mat");
DF3 = num_heave_df;
% load("num_surge_dp.mat");
% B1 = num_surge_dp;
% load("num_surge_df.mat");
% DF1 = num_surge_df;
% load("num_pitch_dp.mat");
% B5 = num_pitch_dp;
% load("num_pitch_df.mat");
% DF5 = num_pitch_df;
eh = zeros(1,50);
ed = zeros(1,50);
cw = zeros(1,50);
for ik = 1:50
    eh(1,ik) = energy_harnessed(B3(ik,2),DF3(ik,2),B5(ik,2),DF5(ik,2),B1(ik,2),DF1(ik,2));
    ed(1,ik) = energy_delivered(1,B3(ik,1)/0.5,1);
    cw(1,ik) = eh(ik)/ed(ik);
end

plot(B3(:,1),cw/(2*a));