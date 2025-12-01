figure(100); clf; hold on;

colors = lines(size(KE_N_array,1));            % color array
legend_entries = cell(size(KE_N_array,1),1);   % store legend strings


%% Convective time

k_steady     = KE_N_array(end,end); %% Kinetic Energy of the whole domain
U_mag_steady = sqrt(2*k_steady);
D            = 6;

time_vec = dt*(1:nsteps);


%% PLOTTING
for j = 1:size(KE_N_array,1)
    
    % Extract the KE curve for this j
    KEj  = KE_N_array(j, :);

    % Plot with unique color
    plot(time_vec, KEj, '-' ,'LineWidth', 3, 'Color', colors(j,:));

    % Store legend text
    % legend_entries{j} = sprintf('$N_{p} = %d$', N_vec(j));
    legend_entries{j} = sprintf('$%d \\times %d$', Nelx, Ny_vec(j));
end


xlabel('t','Interpreter','latex','FontSize',16)
ylabel('$\int_{\Omega}(u^2 + v^2)\, dV$','Interpreter','latex','FontSize',16)

s = sprintf('Gr = %.2f, Pr = %.2f, N = %d', Gr, Pr, N);
title({'\textbf{h-refinement}', s}, 'Interpreter', 'latex', 'FontSize', 20);

% Legend
legend(legend_entries, 'Interpreter','latex', 'Location','best', 'FontSize',14, 'Box','on');

box on
hold off
