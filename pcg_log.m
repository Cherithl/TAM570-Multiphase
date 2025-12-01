clc;close all
pcg_itp_12000 = load("pcg_iter_12000.mat").itp_vec;
pcg_itp_120000 = load("pcg_iter_120000.mat").itp_vec;
pcg_itp_1200000 = load("pcg_iter_1200000.mat").itp_vec;

hdr;
dt=0.5e-03;
Tfinal = 60;
nsteps = ceil(Tfinal/dt);

time_vec = dt*(1:nsteps);

figure
scatter(time_vec,pcg_itp_12000,'o','filled')
hold on
scatter(time_vec,pcg_itp_120000,'o','filled')
scatter(time_vec,pcg_itp_1200000,'o','filled')
legend('Gr = 12000','Gr = 120000','Gr = 1200000',intp,ltx,fs,16)
xlabel('time',fs,16)
ylabel('iterations',fs,16)
title('\textbf{PCG (Pressure solve)}',intp,ltx,fs,20)
box on
