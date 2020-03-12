function [Inflow] = inflowTerm(consts,mats)
        Inflow = squeeze(exp(-(mats.zz(1,1,:,1)-mats.zz(floor(consts.Nx/2),ceil(consts.Ny/2),consts.inflowCenter,1)).^2/2/consts.inflowSpread^2));
        
        Inflow = (Inflow./sum(Inflow)*consts.Nz)';
        Inflow = squeeze(repmat(shiftdim(Inflow,-2),consts.Nx,consts.Ny));
end

