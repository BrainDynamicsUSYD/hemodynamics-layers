function [consts,vecs,mats] = mainMatrix(consts,mats)

vecs.time = linspace(0,consts.timefinal,consts.Nt)-2*consts.deltat;
vecs.spacex = linspace(consts.xmin,consts.xmax,consts.Nx);
vecs.spacey = linspace(consts.ymin,consts.ymax,consts.Ny);
vecs.spacez = linspace(consts.zmin,consts.zmax,consts.Nz);
[mats.xx, mats.yy, mats.zz, mats.tt] = ndgrid(vecs.spacex, vecs.spacey, vecs.spacez, vecs.time); 

end

