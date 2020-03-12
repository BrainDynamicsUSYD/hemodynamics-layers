function Mainfunction
fig6 = figure;
set(fig6,'units','normalized','outerposition',[0 0 1 1]);
hold on
Fnum = 9;
% mags = [];
% maxes = [];
%for v_wave = [0]
D = 200;
Yxtminmin = 0;
Yxtmaxmax = 0;
for v_wave = [0.5:0.1:2];   
%hold on
consts = hemodynamicConstants(v_wave,D); 
for linearFactor = [1]  
consts.linearFactor = linearFactor;
consts.divFactor = 0;
consts.dimRed = 0;
[consts,mats] = createLmatrix(consts);

% Matrix definition

[consts,vecs,mats] = mainMatrix(consts,mats);

mats.Inflow = inflowShape(consts,mats);

[consts,mats] = defineXi0(consts,mats);


if consts.divFactor == 0
    [consts,vecs,mats] = findEquilibrium(consts,vecs,mats);
else
    [consts,vecs,mats] = findEquilibriummodified(consts,vecs,mats);
end

mats.v_b = sqrt(consts.c1.*mats.c2.*mats.beta.*mats.Xi_0.^(mats.beta-1));
consts.v_baverage = mean(mean(mean(mats.v_b)));
consts.C = max(max(max(consts.deltat/consts.N*mats.v_b*(1/consts.deltax+1/consts.deltay+1/consts.deltaz))));

if consts.dimRed == 0

[consts,vecs,mats] = appliedDrive(consts,vecs,mats);

[consts,vecs,mats] = FNormalization(consts,vecs,mats);

if consts.divFactor == 0
    [consts,vecs,mats] = mainHemodynamics(consts,vecs,mats);
else
    [consts,vecs,mats] = mainHemodynamicsmodified(consts,vecs,mats);
end

[consts,vecs,mats,Yxtminmin,Yxtmaxmax] = plottingResults(consts,vecs,mats,Fnum,Yxtminmin,Yxtmaxmax);

else
    
[consts,vecs,mats] = dimReduction(consts,vecs,mats);    
    
end

end
end
end