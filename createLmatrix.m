function [consts,mats] = createLmatrix(consts)
deltakx = 1/(consts.xmax-consts.xmin);
deltaky = 1/(consts.ymax-consts.ymin);
kx = linspace(-(consts.Nx-1)/2*deltakx,(consts.Nx-1)/2*deltakx,consts.Nx);
ky = linspace(-(consts.Ny-1)/2*deltaky,(consts.Ny-1)/2*deltaky,consts.Ny);
[kxM, kyM, kzM] = ndgrid((-1).^[1:consts.Nx],(-1).^[1:consts.Ny],(-1).^[1:consts.Nz]);
mats.kxM = kxM;
mats.kyM = kyM;
mats.kzM = kzM;
L = zeros(consts.Nx,consts.Ny,consts.Nz);
for f = 1:consts.Nx
    for g = 1:consts.Ny
        L(f,g,:) = 4*pi^2*(kx(f)^2+ky(g)^2);
    end
end
mats.L = L;
end

