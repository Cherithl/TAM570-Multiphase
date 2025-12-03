%% 1D CH
close all; clear; clc;

hdr;

N = 100;

L = 2*pi;
Cn = 0.18;

tfinal = 9000;
% nsteps = 2000000;
nsteps = 1e7;
dt = tfinal / nsteps;

[Ah,Bh,Ch,Dh,z,w] = semhat(N);

nh = N+1; Ih = speye(nh);
Rx = Ih;%(:, 2:end); Rx(1, end) = 1; Rx = Rx';

x = L/2*(z+1); Lx = max(x) - min(x); Lx2=Lx/2; Lxi=2/Lx;

Dhx = Lxi * Dh;

Ax = Lxi*Rx*Ah*Rx'; Bx = Lx2*Rx*Bh*Rx'; Bxi = diag(diag(Bx).^(-1));

phi0 = Rx * (cos(2*x) + 1/100*exp(cos(x + 1/10)));

figure
plot(x, phi0)

psi0 = phi0.^3 - phi0;

mu0 = psi0 + (Bx \ (Ax * phi0));

phi = phi0; phi1 = 0*phi; phi2 = phi1; phi3 = phi2;
psi = psi0; psi1 = 0*psi; psi2 = psi1; psi3 = psi2;


LHS = [11/6*Bx, dt*Ax; -Ax*Cn^2, Bx ];

S = (6*dt*Cn^2)/11*Ax*Bxi*Ax + Bx;

figure
for istep = 1:nsteps; time = istep*dt;
    if istep==1; a1=1; a2=0; a3=0; b0=1; b1=1; b2=0; b3=0; end
    if istep==2; a1=2; a2=-1; a3=0; b0=1.5; b1=2; b2=-.5; b3=0; end
    if istep==3; a1=3; a2=-3; a3=1; b0=11/6; b1=3; b2=-1.5; b3=2/6; end
    
    psi3 = psi2; psi2 = psi1; psi1 = phi.^3-phi;
    phi3 = phi2; phi2 = phi1; phi1 = phi;
    
    % rhs1 = Bx * (b1*phi1 + b2*phi2 + b3*phi3);
    % rhs2 = Bx * (a1*psi1 + a2*psi2 + a3*psi3);
    % RHS = [rhs1; rhs2];
    % soln = LHS \ RHS;
    % phi = soln(1:nh);
    % mu = soln(nh+1:end);
    

    rhs_mu  = Bx * (a1*psi1 + a2*psi2 + a3*psi3) + (6*(Cn^2)/11)*Ax*(b1*phi1 + b2*phi2 + b3*phi3);
    mu      = S \ rhs_mu;
    phi     = (6/11)*(b1*phi1 + b2*phi2 + b3*phi3) - ((6/11)*Bxi)*( dt * (Ax* mu));

    if mod(istep, 5000) == 0
        plot(x, phi); xlim([0, 2*pi]); ylim([-1, 1]); title(num2str(time)); drawnow;
    end

end

