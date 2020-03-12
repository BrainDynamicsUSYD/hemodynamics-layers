function generateFigures
% Messing Around
Fnum = 1
D = 200;
Yxtminmin = 0;
Yxtmaxmax = 0;
v_wave = 0;   
consts = hemodynamicConstants(v_wave,D); 
consts.Nt = consts.timepoints/consts.N;
consts.deltat = consts.timefinal/(consts.Nt-1);
consts.deltax = (consts.xmax-consts.xmin)/(consts.Nx-1);
consts.deltay = (consts.ymax-consts.ymin)/(consts.Ny-1);
consts.deltaz = (consts.zmax-consts.zmin)/(consts.Nz-1);
consts.linearFactor = 1;
[consts,mats] = createLmatrix(consts);
[consts,vecs,mats] = mainMatrix(consts,mats);
mats.Inflow = inflowShape(consts,mats);
[consts,mats] = defineXi0(consts,mats);
[consts,vecs,mats] = findEquilibrium(consts,vecs,mats);
mats.v_b = sqrt(consts.c1.*mats.c2.*mats.beta.*mats.Xi_0.^(mats.beta-1));
consts.v_baverage = mean(mean(mean(mats.v_b)));
consts.C = max(max(max(consts.deltat/consts.N*mats.v_b*(1/consts.deltax+1/consts.deltay+1/consts.deltaz))));
[consts,vecs,mats] = appliedDrive(consts,vecs,mats);
[consts,vecs,mats] = FNormalization(consts,vecs,mats);
[consts,vecs,mats] = mainHemodynamics(consts,vecs,mats);
[consts,vecs,mats,Yxtminmin,Yxtmaxmax] = plottingResults(consts,vecs,mats,Fnum,Yxtminmin,Yxtmaxmax);

 end