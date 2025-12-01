function [U,V,T,z,w,Dh,X,Y,Grr,Grs,Gss,Bl,Xr,Rx,Jac,Q,glo_num,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,BC_all,dA]...
             = set_sem_all(N, Nelx, Nely, annulus_flag);

hdr;    % 2-D SEM multi-element

%% Base Geometry
x0 =  0;  x1=3;  Lx=x1-x0;   % Domain coordinates and size
y0 =  0;  y1=1;  Ly=y1-y0;

zc = zwuni(Nelx); xc = x0 + Lx*(zc+1)/2;
zc = zwuni(Nely); yc = y0 + Ly*(zc+1)/2;

%% Problem parameters, as function of N

N1=N+1;
E=Nelx*Nely;

[z,w]=zwgll(N);                % Set basic operators
Dh=deriv_mat(z);

[R,S]=ndgrid(z,z);             % Build SEM mesh
X=zeros(N1,E,N1); Y=X;
e=0; 
for ey=1:Nely
    for ex=1:Nelx; e=e+1;
        xe0=xc(ex); xe1=xc(ex+1); xeL=xe1-xe0;
        ye0=yc(ey); ye1=yc(ey+1); yeL=ye1-ye0;
        X(:,e,:) = xe0 + xeL*(R+1)/2;
        Y(:,e,:) = ye0 + yeL*(S+1)/2;
    end
end

if(annulus_flag)
    [X,Y] = morph_circle(X, Y, Nelx, Nely, N);         % Morph mesh
end

%% Geometric factors, Jacobians and Mass-matrices
[Grr,Grs,Gss,Bl,Xr,Rx,Jac]=geom_elem(X,Y,Dh,w); % Terms for "A"
vol = sum(sum(sum(Bl)));
[Q,glo_num]=set_tp_semq(Nelx,Nely,N);

%% Mask-matrices for Boundary Opeartors
BC_all = [ 'P' 'P' 'D' 'D' ;     %% U
           'P' 'P' 'D' 'D' ;     %% V
           'P' 'P' 'N' 'N' ;     %% P
           'P' 'P' 'D' 'D' ];    %% T

[Mu,Q,glo_num]=set_mask(BC_all(1,:),Nelx,Nely,Q,glo_num); 
[Mv,Q,glo_num]=set_mask(BC_all(2,:),Nelx,Nely,Q,glo_num); 
[Mp,Q,glo_num]=set_mask(BC_all(3,:),Nelx,Nely,Q,glo_num); ifnull=1;
[Mt,Q,glo_num]=set_mask(BC_all(4,:),Nelx,Nely,Q,glo_num); 

[unxa_v,unya_v] = set_unxy(Mu,X,Y,Xr);

dA=diag_sem(Grr,Grs,Gss,Dh); dA=qqt(Q,dA); dA=1./dA;

U = 0 + 0*X;   %% Initial conditions
V = 0 + 0*X;

phi      = 0*X;
mu_field = 0*X;

if(annulus_flag)
    R1 = 1;
    R2 = 3;
    
    R = sqrt(X.^2 + Y.^2);
    T = (R2 - R) ./ (R2 - R1);
else
    T = (y1-Y)./(y1-y0);
end



