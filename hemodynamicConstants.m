function consts = hemodynamicConstants(v_wave,D)

        % Values mainly from Aquino et al 2013, some updated
        
        consts.c1 = 30e-10; % Pressure coupling constant, changed in code to changed v_b range is 1-300e-11
        
        consts.alpha = 0.31; % Grubb's exponent
        
        consts.beta = 1/consts.alpha; % Mean elasticity exponent of cortical vessels
        
        consts.kappa = 0.65; % Blood flow signal decay rate, units of s^(-1)
        
        consts.gamma_lower = 0.4192; % Flow-dependent elimination constant, units of s^(-2)
        
        consts.rho_f = 1062; % Blood mass density, units of kg m^(-3)
        
        consts.D = D; % Effective blood viscosity, range is 106-850, units of kg m^(-3) s^(-1)
        
        consts.rho_hb = 0.4; % Resting blood oxygen extraction fraction
        
        consts.k1 = 4.2; % Magnetic field parameters in the signal equation at 3T, TE = 30ms
        
        consts.k2 = 1.7; % Magnetic field parameters in the signal equation at 3T, TE = 30ms
        
        consts.k3 = 0.41; % Magnetic field parameters in the signal equation at 3T, TE = 30ms
        
        consts.F_0 = 0.01; % Baseline cerebral blood flow, units of s^(-1)
        
        consts.P_0average = 10000;
        
        consts.E_0 = 0.4; % Resting oxygen consumption rate   
        
        consts.v_wave = v_wave;
        
        consts.signal_time = 4;
        consts.Nx = 2^6;
        consts.Ny = 2^1;
        consts.Nz = 2^6;
        consts.timepoints = 2^11;
        consts.N = 2^0;
        consts.N2 = 2^0;
        consts.Nt = consts.timepoints/consts.N;
        consts.timefinal = 20;
        consts.xmax = 0.015;
        consts.xmin = -0.015;
        consts.ymax = 0.005;
        consts.ymin = -0.005;
        consts.zmin = -0.0000;
        consts.zmax = 0.0032;
        consts.deltat = consts.timefinal/(consts.Nt-1);
        consts.deltax = (consts.xmax-consts.xmin)/(consts.Nx-1);
        consts.deltay = (consts.ymax-consts.ymin)/(consts.Ny-1);
        consts.deltaz = (consts.zmax-consts.zmin)/(consts.Nz-1);
        
        consts.l = consts.zmax-consts.zmin; % Cortical thickness, units of m
        
        consts.spread = 0.001; % FWHM of Gaussian is 2.355 times this
        consts.FWantedMagnitude = 2;
        consts.flatC_p = 0;
        
        consts.YDirPresent = 0;
        
        consts.inflowCenter = ceil(0.25*consts.Nz);
        consts.inflowSpread = 0.001;
               
        consts.Xi_peak = ceil(0.433*consts.Nz);
        consts.Xi_lower = 34;
        consts.Xi_middle = 38;
        consts.Xi_upper = 28; % Check this!!! 34 38 24-28
        
end

