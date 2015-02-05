function STATES = InitCond()
    % Below the initial conditions of the differential equation are given.
    % They are chosen, such that the system is in steady state at t=0
    global vivi
    all_indices()
    
    STATES(ind.R_k)     = 0.061e-6;     %'wi in component wi (metre)'
    STATES(ind.N_Na_k)  = 0.99796e-3;   %'N_Nai in component N_Nai (micromolar_metre)'
    STATES(ind.N_K_k)   = 5.52782e-3;   %'N_Ki in component N_Ki (micromolar_metre)'
    STATES(ind.N_HCO3_k)= 0.58804e-3;   %'N_HCO3i in component N_HCO3i (micromolar_metre)'
    STATES(ind.N_Cl_k)  = 0.32879e-3;   %'N_Cli in component N_Cli (micromolar_metre)'
    STATES(ind.N_Na_s)  = 4.301041e-3;  %'N_Nao in component N_Nao (micromolar_metre)'
    STATES(ind.N_K_s)   = 0.0807e-3;    %'N_Ko in component N_Ko (micromolar_metre)'
    STATES(ind.N_HCO3_s)= 0.432552e-3;  %'N_HCO3o in component N_HCO3o (micromolar_metre)'
    STATES(ind.K_p)     = 3e3;         % uM,  [K+] in de perivascular space
    STATES(ind.w_k)     = 0.1815e-3;    % [-]  BK-Channel open probability
   
   
    
    STATES(ind.Ca_i)    = 0.1649;   %*0.1;            % calcium concentration in cytosol
    STATES(ind.s_i)     =  1.361;  %* 0.1          % calcium concentration in sacroplasmatic reticulum
    STATES(ind.v_i)     = -50.3; %*-60;    vivi; %        % mV celmembrane of SMC
    STATES(ind.w_i)     = 0.234; %*0.1;            % open state probability of calcium-activated K channels
    STATES(ind.I_i)     = 0.45; %*0.1;            % IP3 concentration
    
    STATES(ind.K_i)     = 100e3;            %uM [K+] in SMC
    
    STATES(ind.Ca_j)    = 0.5628; %*0.1;            % calcium concentration in EC cytosol
    STATES(ind.s_j)     = 0.8392; %*0.1;            % calcium concentration in endoplasmatic reticulum
    STATES(ind.v_j)     = -65.24; %*-75;            % mV celmembrane of EC
    STATES(ind.I_j)     = 1.35; %*0.1;            % IP3 concentration in EC
    
    STATES(ind.Mp)      = 0.25;  %Mp + M + AMp + AM = 1 !
    STATES(ind.AMp)     = 0.25;
    STATES(ind.AM)      = 0.25;
    
    STATES(ind.R)       = 24.8e-6; % *20e-6; % 15e-6;

          %Hannah:
    STATES(ind.Ca_k)      =0.05e-3;       % uM Bennet 2008
    STATES(ind.s_k)      =0.1e-3;         % uM
    STATES(ind.h_k)      =0.1e-3;         % [-]
    STATES(ind.I_k)      =0.01e-3;        % uM Bennet 2008 
    STATES(ind.EET_k)    =0.1e-3;         % uM
%     
        
    %% NO pathway 
    global NO_switch
    STATES(ind.NOi)     = 0.05;
    STATES(ind.NOj)     = 0.05;
    STATES(ind.NOn)     = 0.1;
    STATES(ind.eNOS_act)= NO_switch*0.7;%3;%!!!!!!!!!!!!!!!!
    STATES(ind.nNOS_act)= NO_switch*0.3;%0.3;%!!!!!!!!!!!!!!!!
    STATES(ind.Ca_n)    = 0.0001;
    STATES(ind.E_b)     = 1/3; % E_b + E_6c + E_5c = 1 !
    STATES(ind.E_6c)     = 1/3;
    STATES(ind.E_5c)    = 1/3;
    STATES(ind.cGMP)    = 8;
    STATES(ind.M_Y)     = 0.5; %M_Y + Mp_Y = 1 !
    STATES(ind.Mp_Y)    = 0.5; 
    STATES(ind.NOa)     = 0.1;
end