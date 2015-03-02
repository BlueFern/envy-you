
%% ODE indices

ind.R_k     = 1;  % m - astrocyte volume-area ratio (AVAR)
ind.N_Na_k  = 2;  % uMm - Astrocyte sodium concentration (times AVAR)
ind.N_K_k   = 3;  % uMm - Astrocyte potassium concentration (times AVAR)
ind.N_HCO3_k= 4;  % uMm - Astrocyte HCO3 concentration (times AVAR)
ind.N_Cl_k  = 5;  % uMm - Astrocyte chloride concentration (times AVAR)
ind.N_Na_s  = 6;  % uMm - synaptic cleft sodium concentration (times AVAR)
ind.N_K_s   = 7;  % uMm - synaptic cleft potassium concentration (times AVAR)
ind.N_HCO3_s= 8;  % uMm - synaptic cleft HCO3 concentration (times AVAR)

ind.K_p     = 9;  % uM -  perivascular space potassium concentration
ind.w_k     = 10; % [-] -  BK-Channel open probability

ind.Ca_i    = 11; % uM - calcium concentration in SMC cytosol
ind.s_i     = 12; % uM - calcium concentration in SMC sacroplasmatic reticulum
ind.v_i     = 13; % mV - celmembrane potential of SMC
ind.w_i     = 14; % [-] - open state probability of calcium-activated K channels
ind.I_i     = 15; % uM - IP3 concentration ins SMC cytosol

ind.K_i     = 16; % uM - SMC potassium concentration

ind.Ca_j    = 17; % uM - calcium concentration in EC cytosol
ind.s_j     = 18; % uM - calcium concentration in EC endoplasmatic reticulum
ind.v_j     = 19; % mV - celmembrane potential of EC
ind.I_j     = 20; % uM - IP3 concentration in EC

ind.Mp      = 21; % [-] - fraction oif free phosphorylated bridges
ind.AMp     = 22; % [-] - fraction of attached phosphorylated bridges
ind.AM      = 23; % [-] - fraction of attached dephosphorylated bridges

ind.R       = 24; % uM - Vessel radius
%indices for astr. Ca2+
ind.Ca_k    = 25; % uM - calcium concentration in the Astrocyte cytosol
ind.s_k     = 26; % uM - calcium concentration in the Astrocyte endoplasmic reticulum 
ind.h_k     = 27; % [-] - probability of Calcium occupying the inhibitory binding site ? :-) 
ind.I_k     = 28; % uM - IP3 concentration in the Astrocyte cytosol
ind.EET_k   = 29; % uM - EET concentration in the Astrocyte cytosol 
%% NO pathway - ODE
i = ind.EET_k ;
ind.Ca_n            = i+1;        % the calcium concentration in the neuron in n
ind.nNOS_act        = i+2;        % activation and deactivation of nNOS
ind.NOn             = i+3;        % NO concentration in the neuron
ind.nNOS_act_max    = i+4;        % maximum activation and deactivation of nNOS
ind.NOn_max         = i+5;        % maximum NO concentration in the neuron
ind.NOk             = i+6;        % NO concentration in the AC
ind.NOk_max         = i+7;        % maximum NO concentration in the AC
ind.NOi             = i+8;        % NO concentration in the SMC
ind.NOi_max         = i+9;        % maximum NO concentration in the SMC
ind.E_b             = i+10;  
ind.E_6c            = i+11;  
ind.E_5c            = i+12; 
ind.cGMP            = i+13;  
ind.eNOS_act        = i+14;        % activation and deactivation of eNOS
ind.NOj             = i+15;        % NO concentration in the EC
ind.eNOS_act_max    = i+16;        % maximum activation and deactivation of eNOS
ind.NOj_max         = i+17;        % maximum NO concentration in the EC
ind.M_Y             = i+18; 
ind.Mp_Y            = i+19; 

%% NO pathway - fluxes NE
j = 1;
flu.tau_w           = j;
flu.P_NR2AO         = j+1 ;
flu.P_NR2BO         = j+2 ; 
flu.openProbTerm    = j+3 ;  
flu.I_Ca            = j+4 ;  
flu.phi_N           = j+5 ;
flu.dphi_N          = j+6 ; 
flu.N               = j+7 ;  
flu.CaM             = j+8 ; 
%flu.NOni              = j+9 ; % NO diffusion flux from the neuron to the SMC

%% Astrocyte indices
flu.R_s     = 1; % m -
flu.N_Cl_s  = 2; % uMm -
flu.Na_k    = 3; % uM - 
flu.K_k     = 4; % uM - 
flu.HCO3_k  = 5; % uM - 
flu.Cl_k    = 6; % uM - 
flu.Na_s    = 7; % uM - 
flu.K_s     = 8; % uM - 
flu.HCO3_s  = 9; % uM - 
flu.Cl_s    = 10; % uM - 

flu.E_Na_k  = 11; % V - 
flu.E_K_k   = 12; % V - 
flu.E_Cl_k  = 13; % V - 
flu.E_NBC_k = 14; % V - 

flu.v_k = 15; % V -      

flu.J_KCC1_k =16; %uMm s-1 - 
flu.J_NBC_k  =17; %uMm s-1 - 
flu.J_NKCC1_k  =18; %uMm s-1 - 
flu.J_NaK_k  =19;  %uMm s-1 - 
flu.J_Na_k   =20; %uMm s-1 - 
flu.J_K_k    =21; %uMm s-1 - 

flu.J_BK_k  =22; %uMm s-1 -  
flu.E_BK_k  =23; % V - 
flu.w_inf   =24; % [-] -
flu.phi_w   =25; % s-1 -

%astrocyte fluxes by Hannah
flu.J_ip3    =26; % uM/s -  calcium flux from ER to cytosolic by IP3 receptors!
flu.J_ERleak =27; % uM/s - calcium leak flux from ER to the cytosol
flu.J_pump   =28; % uM/s - ATP dependent calcium flux from cytoplasm to ER

flu.vh_3    = 29; % mV - calcium dependent open state voltage shift
flu.G_pr    = 30; % [-] - glutamate induced ratio of active to total G-protein receptors
flu.B_cyt   = 31; % [-] - calcium uffer parameter

%% SMC-pointers

flu.v_coup_i        = 1; 
flu.Ca_coup_i       = 2;
flu.IP3_coup_i      = 3;
flu.rho_i           = 4;
flu.J_IP3_i         = 5;
flu.J_SRuptake_i    = 6;
flu.J_CICR_i        = 7;
flu.J_extrusion_i   = 8;
flu.J_leak_i        = 9;
flu.J_VOCC_i        = 10;
flu.J_NaCa_i        = 11;
flu.J_NaK_i         = 12;
flu.J_Cl_i          = 13;
flu.J_K_i           = 14;
flu.Kactivation_i   = 15;
flu.J_degrad_i      = 16;
flu.J_stretch_i     = 17;

flu.v_KIR_i         = 18;
flu.G_KIR_i         = 19;
flu.J_KIR_i         = 20;

flu.M               = 21;
flu.h_r             = 22;
flu.E_K_i           = 23;

flu.K1_c            = 24;
flu.K6_c            = 25;
%% NO pathway - fluxes SMC
k = flu.K6_c;
flu.k4              = k+1;   
flu.R_cGMP1         = k+2;   
flu.R_NO            = k+3;   
flu.v_Ca3           = k+4;  
flu.P_O             = k+5;  
flu.R_cGMP2         = k+6;  
flu.K2_c            = k+7;  
flu.K5_c            = k+8; 
flu.kmlcp           = k+9; 
flu.Act_eNOS_Ca     = k+10;
flu.Act_eNOS_wss    = k+11;
flu.P_NOj_eNOS      = k+12;

%% EC-pointers

flu.v_coup_j         = 1;
flu.Ca_coup_j        = 2;
flu.IP3_coup_j       = 3;
flu.rho_j           = 4;
flu.J_0_j           = 5;
flu.J_IP3_j         = 6;
flu.J_ERuptake_j    = 7;
flu.J_CICR_j        = 8;
flu.J_extrusion_j   = 9;
flu.J_leak_j        = 10;
flu.J_cation_j      = 11;
flu.J_BKCa_j        = 12;
flu.J_SKCa_j        = 13;
flu.J_K_j           = 14;
flu.J_R_j           = 15;
flu.J_degrad_j      = 16;
flu.J_stretch_j     = 17;

%% NO pathway - fluxes EC
l = flu.J_stretch_j;
flu.W_tau_w         = l+1;  
flu.F_tau_w         = l+2;  
% flu.NOji            = l+3; % NO diffusion flux from the EC into the SMC
% flu.NOjl            = l+4; % NO diffusion flux from the EC into the lumen


