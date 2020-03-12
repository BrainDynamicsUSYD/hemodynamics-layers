function [consts,vecs,mats] = FNormalization(consts,vecs,mats)
Fxt = mats.Fxt;
Fxtmax = 0;
for a = 3:consts.Nt;
   Fxt = cat(4, Fxt(:,:,:,2), consts.F_0 + (mats.zxt(:,:,:,a)*consts.deltat^2+(Fxt(:,:,:,2)-consts.F_0)*(2-consts.gamma_lower*consts.deltat^2)+(Fxt(:,:,:,1)-consts.F_0)*(consts.kappa*consts.deltat/2-1))/(1+consts.kappa*consts.deltat/2));
   Fxtmax = max(Fxtmax,max(max(max(max(Fxt)))));
end
Fnormalization = Fxtmax/consts.F_0;
mats.zxt = (consts.FWantedMagnitude-1)/consts.linearFactor/Fnormalization*mats.zxt;
end

