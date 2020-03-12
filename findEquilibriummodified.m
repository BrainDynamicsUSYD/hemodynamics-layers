function [consts,vecs,mats] = findEquilibrium(consts,vecs,mats)

Zxt = squeeze(mats.Xi_0(1,1,:))';
beta = squeeze(consts.beta(1,1,:));

mats.c_p =(consts.F_0*consts.flatC_p+(1-consts.flatC_p)*mats.Inflow.*consts.F_0)/consts.P_0average;
consts.c_pAverage = mean(mean(mean(mats.c_p)));
mats.P_0 = mats.Inflow.*consts.F_0./mats.c_p;

a = -consts.c1/consts.deltaz^2 + vecs.spacez*0;
b = consts.D*squeeze(mats.c_p(1,1,:))' + 2*consts.c1/consts.deltaz^2 + vecs.spacez*0;
c = a;
d = squeeze(consts.F_0*consts.D*mats.Inflow(1,1,:))';

b(1) = 1;
c(1) = -1;        
d(1) = 0;
a(consts.Nz) = 1;
b(consts.Nz) = -1;
d(consts.Nz) = 0;
for k = 2:consts.Nz
    c(k) = c(k)/(b(k)-a(k)*c(k-1));
    d(k) = (d(k)-a(k)*d(k-1))/(b(k)-a(k)*c(k-1));
end
P_0(consts.Nz) = d(consts.Nz);
for k = consts.Nz-1:-1:1
    P_0(k) = d(k) - c(k)*P_0(k+1);
end

mats.P_0 = squeeze(repmat(shiftdim(P_0,-2),consts.Nx,consts.Ny));
mats.c2 = (mats.P_0./(mats.Xi_0.^(mats.beta)));

mats.tau_0 = mats.Xi_0/consts.rho_f./mats.c_p./mats.P_0;

consts.tau_0average = mean(mean(mean(mats.tau_0))); % Hemodynamic transit time, units of s        

vp = -consts.c1/consts.D*gradientv2(P_0,consts.deltaz);
vecs.vp = vp;


kxM = mats.kxM;
kyM = mats.kyM;
kzM = mats.kzM;
Nx = consts.Nx;
Ny = consts.Ny;
Nz = consts.Nz;
deltax = consts.deltax;
deltay = consts.deltay;
deltaz = consts.deltaz;
deltat = consts.deltat;
spacex = vecs.spacex;
spacey = vecs.spacey;
spacez = vecs.spacez;

L = mats.L;

Q1 = consts.eta + 1./squeeze(mats.tau_0(1,1,:))' + consts.rho_f./squeeze(mats.Xi_0(1,1,:))'.*(squeeze(mats.Inflow(1,1,:)*consts.F_0)'-squeeze(mats.c_p(1,1,:))'.*P_0);
Q2 = -consts.rho_f*vp.*gradientv2(squeeze(mats.Xi_0(1,1,:))',consts.deltaz)./squeeze(mats.Xi_0(1,1,:))'.^2;
QInter = Q1 + Q2;
QGen = consts.psi*Zxt*consts.eta;
QGrad = consts.rho_f*vp./squeeze(mats.Xi_0(1,1,:))'/consts.deltaz/2;

a = QGrad;
b = QInter;
c = -QGrad;
d = QGen;

b(1) = 1;
c(1) = -1;        
d(1) = 0;
a(consts.Nz) = 1;
b(consts.Nz) = -1;
d(consts.Nz) = 0;
% c(1) = c(1)./b(1);
% d(1) = d(1)./b(1);

for k = 2:consts.Nz
    c(k) = c(k)/(b(k)-a(k)*c(k-1));
    d(k) = (d(k)-a(k)*d(k-1))/(b(k)-a(k)*c(k-1));
end
Q_0(consts.Nz) = d(consts.Nz);
for k = consts.Nz-1:-1:1
    Q_0(k) = d(k) - c(k)*Q_0(k+1);
end
Q_0 = squeeze(repmat(shiftdim(Q_0,-2),consts.Nx,consts.Ny));

[Xiy,Xix,Xiz] = gradientv2(mats.Xi_0,spacey,spacex,spacez);

for k = 1:consts.Nt*2
    [Qyderiv,Qxderiv,Qzderiv] = gradientv2(Q_0,spacey,spacex,spacez);
    S = mats.Inflow*consts.F_0 - mats.c_p.*mats.P_0;
    
%      b = -reshape(S(2:(Nx-1),2:(Ny-1),2:(Nz-1)),(Nx-2)*(Ny-2)*(Nz-2),1);
%         
%         Y = X\b;
%         
%         vfp = reshape(Y,Nx-2,Ny-2,Nz-2);
% 
%         vfp = padarray(vfp,[1 1 1]);
%         
%         [vyfp,vxfp,vzfp] = gradientv2(vfp,deltay,deltax,deltaz);
    
    Skxkyz = kxM.*(fft(kxM.*kyM.*fft(kyM.*S,[],2),[],1));
                v = 1/deltaz^2+zeros(Nx,Ny,Nz);
                b = -2*v-L;
                c = v;
                d = Skxkyz;
                  for r = 1:Nz
                    if r == 1
                        b(:,:,r) = 1;
                        c(:,:,r) = -1;
                        d(:,:,r) = 0;
                    elseif r == Nz
                        v(:,:,r) = 1;
                        b(:,:,r) = -1;
                        d(:,:,r) = 0;
                    else
                        c(:,:,r) = c(:,:,r)./(b(:,:,r)-v(:,:,r).*c(:,:,r-1));
                        d(:,:,r) = (d(:,:,r)-v(:,:,r).*d(:,:,r-1))./(b(:,:,r)-v(:,:,r).*c(:,:,r-1));
                    end
                end
                Tkxkyz(:,:,Nz) = d(:,:,Nz);
                for r = Nz-1:-1:1
                   Tkxkyz(:,:,r) = d(:,:,r)-c(:,:,r).*Tkxkyz(:,:,r+1);        
                end
        Y = real(kxM.*ifft(kxM.*kyM.*ifft(kyM.*Tkxkyz,[],2),[],1));
        
        [vyfp,vxfp,vzfp] = gradientv2(Y,spacey,spacex,spacez);                             
   
   Qxt1 = consts.eta+mats.Inflow*consts.F_0*consts.rho_f./mats.Xi_0;
   Qxt2 = -(vxfp.*Xix+vyfp.*Xiy+vzfp.*Xiz).*consts.rho_f./mats.Xi_0.^2.;
   QxtIntermediate = Qxt1 + Qxt2;
   QxtGeneration = consts.eta*consts.psi*mats.Xi_0;
   Qxtgradientv2 = consts.rho_f*(vxfp.*Qxderiv+vyfp.*Qyderiv+vzfp.*Qzderiv./mats.Xi_0);
   
   Q_0 = (Q_0.*(1-deltat*QxtIntermediate/2)+deltat*(QxtGeneration-Qxtgradientv2))./(1+deltat*QxtIntermediate/2);  
   Q_0(:,:,1) = Q_0(:,:,2);
   Q_0(:,:,Nz) = Q_0(:,:,Nz-1);
end

mats.Q_0 =Q_0;
end

