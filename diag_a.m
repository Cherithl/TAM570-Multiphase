function da = diag_a(Grr,Grs,Gss,Dr,Ds);

nr=size(Dr,1);
ns=size(Ds,1);

GRR = zeros(size(Grr,1));
GSS = zeros(size(Gss,1));

GRR(:,:,1) = Grr;
GSS(:,:,1) = Gss;

da = (Dr.*Dr)'*GRR + GSS*(Ds.*Ds);

da(  1,  1) = da(  1,  1)+Grs(  1,  1)*Dr(  1,  1)*Ds(  1,  1)*2;
da(end,  1) = da(end,  1)+Grs(end,  1)*Dr(end,end)*Ds(  1,  1)*2;
da(  1,end) = da(  1,end)+Grs(  1,end)*Dr(  1,  1)*Ds(end,end)*2;
da(end,end) = da(end,end)+Grs(end,end)*Dr(end,end)*Ds(end,end)*2;


