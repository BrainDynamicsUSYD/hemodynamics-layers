function [consts,mats] = defineXi0(consts,mats)

Nx = consts.Nx;
Ny = consts.Ny;
Nz = consts.Nz;

Xi_lower = consts.Xi_lower;
Xi_upper = consts.Xi_upper;
Xi_middle = consts.Xi_middle;
Xi_peak = consts.Xi_peak;

Xi_0(1) = Xi_lower;
Xi_0(Nz) = Xi_upper;
Xi_0(2:Xi_peak) = Xi_lower + (0:Xi_peak-2)/(Xi_peak-2)*(Xi_middle-Xi_lower);
Xi_0(Xi_peak:Nz-1) = Xi_middle + (0:Nz-1-Xi_peak)/(Nz-1-Xi_peak)*(Xi_upper-Xi_middle);

mats.Xi_0 = squeeze(repmat(shiftdim(Xi_0,-2),Nx,Ny));
mats.v_0 = mats.Xi_0/consts.rho_f;

mats.beta = consts.beta - 0*mats.Xi_0;

consts.Xi_0average = mean(mean(mean(mats.Xi_0))); % Resting mass density contributed by blood, units of kg m^(-3)
consts.v_0average = consts.Xi_0average/consts.rho_f; % Resting blood volume fraction
        
consts.hemodensity = 0.056;

consts.psi = consts.hemodensity/consts.Xi_0average; % Ratio of hemoglobin concentration to blood density, units of mol kg^(-1)
  

consts.Q_0 = consts.E_0*consts.hemodensity;  

consts.eta = consts.rho_f*consts.F_0/consts.Xi_0average/(consts.hemodensity/consts.Q_0-1); % Oxygen consumption rate, units of s^(-1), theoretically rho_hb/tau



end

