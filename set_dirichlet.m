function [U,V,T] = set_dirichlet(Uo,Vo,To,Mu,Mv,Mt,X,Y);

eps = 1e-14;
R1 = 1;
R2 = 3;
R = sqrt(X.^2 + Y.^2);
Min  = ((R < R1+eps) & (R > R1-eps));
Mout = ((R < R2+eps) & (R > R2-eps));
Mt = 1 - (Min + Mout);

T1 = 1;
T2 = 0;

U = 0 + 0*X;
V = 0 + 0*X;  %% Desired field at inflow
T = 0 + 0*X;  %% Desired field at inflow

T(Min) = T1;
T(Mout) = T2;

U = Mu.*Uo + (1-Mu).*U;   %% Old value in interior, new value on inflow
V = Mv.*Vo + (1-Mv).*V;
T = Mt.*To + (1-Mt).*T;

% figure(5);
% for e=1:32
%     mesh(reshape(X(:,e,:), 9, 9), reshape(Y(:,e,:), 9, 9), reshape(Mt(:,e,:), 9, 9)); hold on
% end

% se_mesh(X,Y,U,'U bdry'); pause(1); pause



























