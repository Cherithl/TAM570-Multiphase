close all; clear; clc;

hdr;    % 2-D SEM multi-element

N=8;
Nelx = 6;  Nely = 2;

Gr = 1200000;
Pr = 0.8;
alpha = 1/(sqrt(Gr)*Pr);
Re = sqrt(Gr);
nu = 1 / Re;

% inner and outer radii
R1 = 1;
R2 = 3;

T1 = 1;
T2 = 0;

%% Set ICs, problem parameters as function of N
[U,V,T,z,w,Dh,X,Y,Grr,Grs,Gss,Bl,Xr,Rx,Jac,Q,glo_num,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,BC_all,dA]=set_sem_all(N, Nelx, Nely);


[X_glob] = stitch(X, Nelx,Nely, N);
[Y_glob] = stitch(Y, Nelx,Nely, N);

R_vec = Y_glob(:,1);
dR_in = R_vec(2)-R_vec(1);
dR_out = R_vec(end)-R_vec(end-1);
d = R_vec(end) - R_vec(1);

theta_vec = atan2(Y_glob(1,:),X_glob(1,:))*180/pi;
theta_vec = mod(theta_vec, 360);
[theta_sorted, idx] = sort(theta_vec);

%% Set dealiasing operators, JM,DM,BMh
Nd = floor(1.5*N);
[zd,wd]=zwgl(Nd); Bd=diag(wd);
JM=interp_mat(zd,z);
DM=deriv_mat (zd);
BMh=tensor3(Bd*JM,1,Bd*JM,1+0*Jac); %% effectively, Bd*JM*Jac*JM'*Bd'

%% Set plotting interpolation operator
Nf = floor(1.2*N); [zf,wf]=zwuni(Nf); Jf=interp_mat(zf,z);

%% Set timestep size
[lam_max_est]=est_lam_cfl(U,V,Rx); ldt_max = 0.5;
lam_max_est=max(1e-5,lam_max_est);
dt=ldt_max/lam_max_est;

dt=0.5e-03;
Tfinal = 60;
nsteps = ceil(Tfinal/dt);
%nsteps = 1;

%% Initialize BDFk/EXTk arrays
O=0*X; F=O; G=O; H=O;
U1=O;U2=O;U3=O; V1=O;V2=O;V3=O;
F1=O;F2=O;F3=O; G1=O;G2=O;G3=O;
f1=O;f2=O;f3=O; g1=O;g2=O;g3=O;
T1=T;T2=T;T3=T; H1=O;H2=O;H3=O;

P=O;  %% Initialize pressure to zero

%% System-solve parameters
ifnull=0; tol=1.e-6; max_iter=140;

%%%%% TIME STEPPING LOOP %%%%%%

ke_vec   = zeros(nsteps,1);
time_vec = zeros(nsteps,1);

itp_vec = [];
Nu_in_integral_arr = [];
Nu_out_integral_arr = [];


kk=0; k=0; time=0;
for istep =1:nsteps; k=k+1;

    %%   Set dt and time
    %    [lam_max_est]=est_lam_cfl(U,V,Rx); ldt_max = 0.5; % Variable dt?
    %    dt=ldt_max/lam_max_est;
    time = time+dt;

    %%   Set updated BDFk/EXTk coefficients
    ndt = nu*dt; adt = alpha*dt;
    if k==1; a1=1; a2=0; a3=0; b0=1; b1=1; b2=0; b3=0; end;
    if k==2; a1=1.5; a2=-.5; a3=0; b0=1.5; b1=2; b2=-.5; b3=0; end;
    if k==3; a1=3; a2=-3; a3=1; b0=11/6; b1=3; b2=-1.5; b3=2/6; end;
    d1=dt*a1; d2=dt*a2; d3=dt*a3;

    %%   Set dealiased advecting field
    [Cr,Cs]=set_advect_c(U,V,JM,BMh,Jac,Rx);

    %%   Set body force, volumetric heating
    if k==1; QT =  0*Y; end;  %% No forcing, no heating
    if k==1; FX =  0*Y; end;
    if k==1; FY =  0*Y; end;
    FY = T;


    %%   Evaluate curl-curl term (to be extrapolated)
    [curlcurlX,curlcurlY,Omega]=curlcurl(U,V,Bl,Rx,Dh);
    % curlcurlX = 0*X;
    % curlcurlY = 0*Y;
    % Omega = Lxi*(Dhx*V) - Lyi*(U*Dhy');
    % curlcurlX =  Bl.*(Lyi*(Omega*Dhy'));
    % curlcurlY = -Bl.*(Lxi*(Dhx*Omega));

    %
    %    Set Dirichlet conditions onto old fields
    %
    [Ub,Vb,Tb]=set_dirichlet(U,V,T,Mu,Mv,Mt,X,Y);


    %%   Compute u-hat and u-tilde
    U3=U2;U2=U1;U1=U;
    F3=F2;F2=F1;F1=-advectl(U,Cr,Cs,JM,DM)+Bl.*FX;
    f3=f2;f2=f1;f1=-nu*curlcurlX;
    Uh=Bl.*(b1*U1+b2*U2+b3*U3)+(d1*F1+d2*F2+d3*F3);
    Ut=Uh+(d1*f1+d2*f2+d3*f3);
    Uh=Uh-axl(Ub,b0,ndt,Bl,Grr,Grs,Gss,Dh);

    %%   Compute v-hat and v-tilde
    V3=V2;V2=V1;V1=V;
    G3=G2;G2=G1;G1=-advectl(V,Cr,Cs,JM,DM)+Bl.*FY;
    g3=g2;g2=g1;g1=-nu*curlcurlY;
    Vh=Bl.*(b1*V1+b2*V2+b3*V3)+(d1*G1+d2*G2+d3*G3);
    Vt=Vh+(d1*g1+d2*g2+d3*g3);
    Vh=Vh-axl(Vb,b0,ndt,Bl,Grr,Grs,Gss,Dh);


    %%   Compute t-hat
    T3=T2;T2=T1;T1=T;
    H3=H2;H2=H1;H1=-advectl(T,Cr,Cs,JM,DM)+Bl.*QT;
    Th=Bl.*(b1*T1+b2*T2+b3*T3)+(d1*H1+d2*H2+d3*H3);
    Th=Th-axl(Tb,b0,adt,Bl,Grr,Grs,Gss,Dh);

    %    Pressure correction
    divUt = weak_div(Ut,Vt,1.,Rx,Dh)/dt;
    
    %%   Add inhomogeneous Neumann data to divUT, if any. (Eq.(15) in split_slides.pdf)
    b0dt = b0/dt;
    divUt = divUt - b0dt*( (1-Mu).*unxa_v.*Ub - (1-Mv).*unya_v.*Vb );

    %%   Pressure-Poisson solve
    h1=1; h0=0;
    divUt = divUt -axl(P,h0,h1,Bl,Grr,Grs,Gss,Dh);
    [dP,itp,res,lamda_h]=...
        pcg_lambda(divUt,tol,max_iter,h0,h1,Mp,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull);
    s=['Pressure. Step/Iter: = ' int2str([istep itp])];
    %    hold off; se_mesh  (X,Y,dP,s);  drawnow;
    P = P+dP;
    itp_vec = [itp_vec, itp];

    [dPdx,dPdy]=grad(P,Rx,Dh);
    Uh = Uh - dt*Bl.*dPdx;
    Vh = Vh - dt*Bl.*dPdy;

    %    Viscous/diffusive solves (diagonally-preconditioned CG):

    %%   Implicit solve - diagonal preconditioner
    dAT=1./dA; dAT=1./(b0*qqt(Q,Bl)+adt*dAT);
    dAU=1./dA; dAU=1./(b0*qqt(Q,Bl)+ndt*dAU);
    [U,itu,res,lamda_h]=...
        pcg_lambda(Uh,tol,max_iter,b0,ndt,Mu,Q,Bl,Grr,Grs,Gss,Dh,dAU,ifnull);
    [V,itv,res,lamda_h]=...
        pcg_lambda(Vh,tol,max_iter,b0,ndt,Mv,Q,Bl,Grr,Grs,Gss,Dh,dAU,ifnull);

    [T,itt,res,lamda_h]=...
        pcg_lambda(Th,tol,max_iter,b0,adt,Mt,Q,Bl,Grr,Grs,Gss,Dh,dAT,ifnull);

    U=U+Ub;  %% Add back any prescribed Dirichlet conditions
    V=V+Vb;
    T=T+Tb;

    % time_vec(istep) = time;
    % ke_vec(istep)   = abs(sum(sum(sum( 0.5 * ( (U.*U).*Bl  + (V.*V).*Bl )))));

    %    Diagonostics
    if mod(istep,1000)==0 || istep==1;  kk=kk+1;
        
time
        [U_glob] = stitch(U, Nelx,Nely, N);
        [V_glob] = stitch(V, Nelx,Nely, N);
        [T_glob] = stitch(T, Nelx,Nely, N);

        %      disp([itp itu itv itt])

        hold off;

        % umax = max(max(max(abs(U))));
        % vmax = max(max(max(abs(V))));
        % tmax = max(max(max(abs(T))));
        % um(kk) = umax; vm(kk) = vmax;
        % tm(kk) = tmax; ti(kk) = time;
        % 
        % Uf = tensor3(Jf,1,Jf,U);  Uf=U;
        % Vf = tensor3(Jf,1,Jf,V);  Vf=V;
        % Xf = tensor3(Jf,1,Jf,X);  Xf=X;
        % Yf = tensor3(Jf,1,Jf,Y);  Yf=Y;
        % Tf = tensor3(Jf,1,Jf,T);  Tf=P;

        
        % se_mesh  (Xf,Yf,Tf,s);  hold on;
        % hold off; 
        % se_quiver(Xf,Yf,Uf,Vf,s,T); axis([-R2, R2, -R2, R2]); hold on; %saveas(gcf, sprintf('T_contour_%d.png', istep));
        % drawnow
        % disp([umax vmax tmax])

    figure(2);
        s=['time=' num2str(time) ,', pcg-iter= ' num2str(itp)'.'];
        contourf(X_glob, Y_glob, T_glob, 80, "Linecolor","None"); 
        c=colorbar; c.Label.String = 'Temperature'; c.Label.FontSize=13;   colormap("hot")
        hold on
        clim([0,1])
        % quiver(X_glob,Y_glob,U_glob,V_glob,1,'LineWidth',1.5,'Color', 'k'); hold on;
        axis equal
        title(s,fs,16);
        xlabel('',fs,16);
        ylabel('',fs,16);
    % 
    % saveas(gcf, sprintf('T_12000_contour_%d.png', istep));
    % 
        hold off
        drawnow
    % 
    % 
    % figure(3);
    %     s=['time=' num2str(time) ,', pcg-iter= ' num2str(itp)'.'];
    %     contourf(X_glob, Y_glob, T_glob, 80, "Linecolor","None"); 
    %     c=colorbar; c.Label.String = 'Temperature'; c.Label.FontSize=13;   colormap("hot")
    %     hold on
    %     clim([0,1])
    %     quiver(X_glob,Y_glob,U_glob,V_glob,"off",'LineWidth',1,'Color', 'k'); hold on;
    %     axis equal
    %     title(s,fs,16);
    %     xlabel('',fs,16);
    %     ylabel('',fs,16);
    % 
    % saveas(gcf, sprintf('T_12000_quiver_%d.png', istep));
    % 
    %     hold off
    %     drawnow
        

    % figure(10)
    % hold on
    % plot(time_vec(istep),ke_vec(istep),'or','LineWidth',2)

    end

[T_glob] = stitch(T, Nelx,Nely, N);   
Nu_in_vec = d*abs( ( T_glob(2,:) - T_glob(1,:) )/dR_in );
Nu_out_vec = d*abs( ( T_glob(end,:) - T_glob(end-1,:) )/dR_in );

Nu_in_sorted  = Nu_in_vec(idx);
Nu_out_sorted = Nu_out_vec(idx);

theta_rad = deg2rad(theta_sorted);
Nu_in_integral_arr  = [Nu_in_integral_arr, trapz(theta_rad, Nu_in_sorted)/(2*pi)];
Nu_out_integral_arr = [Nu_out_integral_arr, trapz(theta_rad, Nu_out_sorted)/(2*pi)];




end


% save('Nu_12000.mat','Nu_in_integral_arr','Nu_out_integral_arr');

return
%% plotting
figure(7);
E=size(X,2);
N1=size(X,1);

% hold off;
for e=1:E
    x=X(:,e,:); x=reshape(x,N1,N1);
    y=Y(:,e,:); y=reshape(y,N1,N1);
    u=U(:,e,:); u=reshape(u,N1,N1);
    v=V(:,e,:); v=reshape(v,N1,N1);
    %quiver(x,y,u,v,'off','Color', 'k', 'Linewidth', 2.5); hold on;
    contourf(x, y, reshape(T(:, e, :), N1, N1), 40, Linecolor='None'); axis equal; 
    h = colorbar; cm = colormap('gray'); %colormap(flipud(cm));
    %mesh(x, y, reshape(T(:, e, :), N1, N1)); hold on;
    quiver(x,y,u,v,'off','Color', 'g');
end
xlabel('X',fs,16);
ylabel('Y',fs,16);
% view(55,33);

