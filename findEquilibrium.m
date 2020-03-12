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

Q_0 = consts.psi*mats.Xi_0./(1+consts.rho_f.*mats.c_p.*mats.P_0./consts.eta./mats.Xi_0);

mats.Q_0 =Q_0;
end

