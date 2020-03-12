function [consts,vecs,mats] = mainHemodynamics(consts,vecs,mats)

Fxt = mats.Fxt;
Zxt = cat(4,mats.Xi_0,mats.Xi_0);
Zxtswap1 = Zxt(:,:,:,1);
Zxtswap2 = Zxt(:,:,:,2);
Zuse2 = Zxt(:,:,:,1);
Pnew = mats.P_0;
P_0 = Pnew;
Qxt = mats.Q_0;
tau = mats.tau_0;
Yxtout = 0*Qxt;
Pderivnew = 6*del2(Pnew,vecs.spacey,vecs.spacex,vecs.spacez);

Nt = consts.Nt;
Nx = consts.Nx;
Ny = consts.Ny;
Nz = consts.Nz;
N = consts.N;
F_0 = consts.F_0;
deltat = consts.deltat;
deltax = consts.deltax;
deltay = consts.deltay;
deltaz = consts.deltaz;
D = consts.D;
rho_f = consts.rho_f;
c1 = consts.c1;
eta = consts.eta;
psi = consts.psi;

gamma_lower = consts.gamma_lower;
kappa = consts.kappa;

spacex = vecs.spacex;
spacey = vecs.spacey;
spacez = vecs.spacez;
time = vecs.time;

    kxM = mats.kxM;
    kyM = mats.kyM;
    kzM = mats.kzM;

    L = mats.L;

c2 = mats.c2;
c_p = mats.c_p;
beta = mats.beta;

central_time = 5;
sigma_x = 0.001;
sigma_y = 0.001;
sigma_t = 0.5;

[x y z t] = ndgrid(vecs.spacex,vecs.spacey,vecs.spacez,vecs.time);

for a = 3:Nt;
   Zuse1 = Zuse2;
   Zuse2 = 0*Zuse2;
   Zxt(:,:,:,1) = Zxtswap1;
   Zxt(:,:,:,2) = Zxtswap2;
   for b = 1:N;         
        Fxt = cat(4, Fxt(:,:,:,2), F_0 + (mats.zxt(:,:,:,a)*deltat^2/N^2+(Fxt(:,:,:,2)-F_0)*(2-gamma_lower*deltat^2/N^2)+(Fxt(:,:,:,1)-F_0)*(kappa*deltat/N/2-1))/(1+kappa*deltat/N/2));    
        Pold = Pnew;
        Zxt(:,:,1,2) = Zxt(:,:,2,2);
        Zxt(:,:,Nz,2) = Zxt(:,:,Nz-1,2);
        Pnew = c2.*Zxt(:,:,:,2).^beta;
        Pderivold = Pderivnew;
        Pderivnew = 6*del2(Pnew,spacey,spacex,spacez);    
        Zxtintermediate = mats.Inflow.*(Fxt(:,:,:,2)*(deltat^2/N^2*D*3/2+rho_f*deltat/N)+Fxt(:,:,:,1)*(-D*deltat^2/N^2/2-rho_f*deltat/N))-(c_p.*(Pnew.*(D*deltat^2/N^2*3/2+rho_f*deltat/N)+Pold.*(-D*deltat^2/N^2/2-rho_f*deltat/N)));
        Zxt = cat(4,Zxt(:,:,:,2),(Zxtintermediate+c1*(3/2*Pderivnew-1/2*Pderivold)*deltat^2/N^2+2*Zxt(:,:,:,2)+Zxt(:,:,:,1)*(D*deltat/N/2/rho_f-1))./(1+D*deltat/N/2/rho_f));
        
        Zxt(:,:,1,2) = Zxt(:,:,2,2);
        Zxt(:,:,Nz,2) = Zxt(:,:,Nz-1,2);
        Zuse2 = Zuse2 + Zxt(:,:,:,2);
   end
   Zuse2 = Zuse2/N;
   Zxtswap1 = Zxt(:,:,:,1);
   Zxtswap2 = Zxt(:,:,:,2);
   Zxt(:,:,:,1) = Zuse1;
   Zxt(:,:,:,2) = Zuse2;      
        
   tau(:,:,:,a) = (Zxt(:,:,:,2)+Zxt(:,:,:,1))./Pnew/rho_f./c_p/2;
   
   QxtIntermediate = eta + mats.c_p.*Pnew.*rho_f./(Zxt(:,:,:,2)+Zxt(:,:,:,1))*2;
   QxtGeneration = eta.*psi.*(Zxt(:,:,:,2)+Zxt(:,:,:,1))/2;
   Qxt = (Qxt.*(1-deltat*QxtIntermediate/2)+deltat*QxtGeneration)./(1+deltat*QxtIntermediate/2); 
   %Qxt(:,:,1) = Qxt(:,:,2);
   %Qxt(:,:,Nz) = Qxt(:,:,Nz-1);
   Yxt = consts.linearFactor*mats.Xi_0/rho_f.*(consts.k1*(1-Qxt./mats.Q_0)+consts.k2*(1-Qxt./mats.Q_0./Zxt(:,:,:,2).*mats.Xi_0)+consts.k3*(1-Zxt(:,:,:,2)./mats.Xi_0));     

   if mod(a,consts.N2) == 0;
    Fxtout(:,:,:,a/consts.N2) = (Fxt(:,:,:,2)-F_0)*consts.linearFactor+F_0;
    Zxtout(:,:,:,a/consts.N2) = ((Zxt(:,:,:,2)-mats.Xi_0)*consts.linearFactor+mats.Xi_0);
    Qxtout(:,:,:,a/consts.N2) = ((Qxt-mats.Q_0)*(consts.linearFactor)+mats.Q_0)./consts.hemodensity;
    Yxtout(:,:,:,a/consts.N2) = Yxt; %(Pnew-c2.*mats.Xi_0.^(consts.beta-1).*(mats.Xi_0+consts.beta*(Zxt(:,:,:,2)-mats.Xi_0)))./mats.P_0;%consts.linearFactor*(Zxtintermediate+c1*(3/2*Pderivnew-1/2*Pderivold));
    %Yxt;
    Qxtterm1out(:,:,:,a/consts.N2) = Yxt;%Qxt.*(Zxt(:,:,:,2)-Zxt(:,:,:,1))/deltat*2./(Zxt(:,:,:,2)+Zxt(:,:,:,1));
    Qxtterm2out(:,:,:,a/consts.N2) = Yxt;%-rho_f*Qxt.*mats.Inflow.*Fxt(:,:,:,2)*2./(Zxt(:,:,:,2)+Zxt(:,:,:,1)) + (psi*(Zxt(:,:,:,2)+Zxt(:,:,:,1))/2-Qxt)*eta;
    Qxtterm3out(:,:,:,a/consts.N2) = Yxt;%Qxt.*(-Qxt2)-Qxtgradient;
    Qxtterm4out(:,:,:,a/consts.N2) = Yxt;%Qxt.*(Zxt(:,:,:,2)-Zxt(:,:,:,1))/deltat*2./(Zxt(:,:,:,2)+Zxt(:,:,:,1)) -rho_f*Qxt.*mats.Inflow.*Fxt(:,:,:,2)*2./(Zxt(:,:,:,2)+Zxt(:,:,:,1)) + (psi*(Zxt(:,:,:,2)+Zxt(:,:,:,1))/2-Qxt)*eta + Qxt.*(-Qxt2)-Qxtgradient;
   end
end


mats.Fxt = Fxtout;
mats.Zxt = Zxtout;
mats.Qxt = Qxtout;
mats.Yxt = Yxtout;
mats.Qxtterm1 = Qxtterm1out;
mats.Qxtterm2 = Qxtterm2out;
mats.Qxtterm3 = Qxtterm3out;
mats.Qxtterm4 = Qxtterm4out;

end