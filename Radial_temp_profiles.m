clc; close all
hdr;

rad_vec = (Y_glob(:,1) - 1)/2;

figure
plot(rad_vec, T_glob(:,1),'-',lw,2,'Color','k')
hold on
plot(rad_vec, T_glob(:,9),'--',lw,2,'Color','k')
plot(rad_vec, T_glob(:,17),'-.',lw,2,'Color','k')
plot(rad_vec, T_glob(:,25),':',lw,2,'Color','k')
legend('$\theta = 0^\circ$','$\theta = 60^\circ$','$\theta = 120^\circ$','$\theta = 180^\circ$',intp,ltx,fs,13)
xlabel('$(R-R_i)/(R_0-R_i)$',intp,ltx,fs,16)
ylabel('$(T-T_i)/(T_0-T_i)$',intp,ltx,fs,16)
box on
% plot(T_glob(:,33),'o')
% plot(T_glob(:,41),'o')
