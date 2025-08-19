clear variables
close all
clc

load('data_r1.mat'); load('data_r2.mat');
load('data_r3.mat'); load('data_r4.mat');
load('data_r5.mat');

figure()
loglog(DOF_L_r1,err_M_r1,DOF_L_r2,err_M_r2,DOF_L_r3,err_M_r3,DOF_L_r4,err_M_r4,DOF_L_r5,err_M_r5)
legend("r1","r2","r3","r4","r5"),xlabel("Total DOF")
xlim([10 1000])

figure()
loglog(N0_L_r1,err_M_r1,N0_L_r2,err_M_r2,N0_L_r3,err_M_r3,N0_L_r4,err_M_r4,N0_L_r5,err_M_r5)
legend("r1","r2","r3","r4","r5"),xlabel("DOF (no interpoation)")
xlim([10 200])
