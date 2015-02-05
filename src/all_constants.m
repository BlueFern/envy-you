global NVU
%% General constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Farad       = 96500         ;% [C mol-1]      Faradays constant
R_gas       = 8.315         ;% [J mol-1K-1]
Temp        = 300           ;% [K]

unitcon    = 10^(3)         ;% [-]            Factor to convert equations to another unit

%% Constants of Ostby, the values come from CELL ML
L_p         = 2.1e-9        ;% [m uM-1s-1]
R_tot       = 8.79e-8       ;% [m]
X_k         = 12.41e-3      ;% [uMm]
z_Na        = 1             ;% [-]
z_K         = 1             ;% [-]
z_Cl        = -1            ;% [-]
z_NBC       = -1            ;% [-]
g_K_k       = 40            ;% [ohm-1m-2]
g_KCC1_k    = 1e-2          ;% [ohm-1m-2]
g_NBC_k     = 7.57e-1       ;% [ohm-1m-2]
g_Cl_k      = 8.797e-1      ;% [ohm-1m-2]
g_NKCC1_k   = 5.54e-2       ;% [ohm-1m-2]
g_Na_k      = 1.314         ;% [ohm-1m-2]
J_NaK_max   = 1.42e-3       ;% [uMm s-1]
K_Na_k      = 10e3          ;% [uM]
K_K_s       = 1.5e3         ;% [uM]
k_C         = 7.35e-5       ;% [muM s-1]



%% Constants of the BK-channel

A_ef_k      = 3.7e-9        ;% m2       Area of an endfoot of an astrocyte, equal to Area astrocyte at synaptic cleft
v_6         = 22e-3         ;% V        Constant
v_4         = 14.5e-3       ;% V        A measure of the spread of the distrubution
psi_w       = 2.664         ;% s-1      A characteristic time
G_BK_k      = 4.3e3         ; % pS      Constant estimation based on Ermentrout
g_BK_k      = G_BK_k*10^(-12)/A_ef_k ;% ohm-1m-2  Specific capacitance of the BK-Channel in units of Ostby
VR_pa       = 0.001       ;% [-]       The estimated volume ratio of perivascular space to astrocyte: Model estimation
VR_ps       = 0.001       ;% [-]       The estimated volume ratio of perivascular space to SMC: Model Estimation	




%% SMC constants


F_il = 7.5e2            ;%[-] scalingsfactor to fit the experimental data of Filosa
z_1 =4.5                ;%[-] parameter fitted on experimental data of Filosa
z_2 =-1.12e2            ;%[-] parameter fitted on experimental data of Filosa
z_3 =4.2e-1             ;%[-] parameter fitted on experimental data of Filosa
z_4 =-1.26e1            ;%[-] parameter fitted on experimental data of Filosa
z_5 =-7.4e-2            ;%[-] parameter fitted on experimental data of Filosa

% mVmicroM-1 The change in membrane potential by a scaling factor

% Koeningsberger et al.

Fmax_i		= 0.23;		% [microM/s]
Kr_i 		= 1; 		% [microM] ; Half saturation constant for agonist-dependent Ca$^{2+}$ entry
G_Ca		= 0.00129;	% [microM/mV/s]
v_Ca1		= 100;		% [mV]
v_Ca2		= -24;	    % [mV]  why did we change it to -35 temporarely???
R_Ca		= 8.5;		% [mV]
G_NaCa		= 0.00316;	% microM/mV/s
c_NaCa		= 0.5;		% microM
v_NaCa		= -30;
B_i			= 2.025;
cb_i		= 1;
C_i			= 55;
sc_i		= 2;
cc_i		= 0.9;
D_i			= 0.24;
vd_i		= -100;
Rd_i		= 250;
L_i			= 0.025;
gam			= 1970; % mVmicroM-1 The change in membrane potential by a scaling factor
F_NaK		= 0.0432;
G_Cl		= 0.00134;
v_Cl		= -25;
G_K			= 0.00446;
vK_i		= -94; 
lab 		= 45;
%c_w			= 0;
bet			= 0.13;

R_K			= 12;
k_i			= 0.1;
% K_d         = 1;            % = 1e3 nM Gonzalez
% B_T         = 100;          % = 1e5 nM Gonzalez
% Kinf_i      = 1e5;          % 100 mM K+ concentration in SMC

G_stretch   = 0.0061;       % uM mV-1 s-1
P_str       = 30;
Esac        = -18;          % mV
alpha1      = 0.0074;
sig0        = 500;


%% EC constants
% EC Koeningsberger et al.

Fmax_j		= 0.23;		% [microM/s]
Kr_j		= 1;
B_j 		= 0.5;
cb_j		= 1;
C_j			= 5;
sc_j		= 2;
cc_j		= 0.9;
D_j			= 0.24;
L_j			= 0.025;
G_cat 		= 0.66e-3;
E_Ca		= 50;
m3cat		= -0.18; %-6.18; %changed value!!! 
m4cat 		= 0.37;
J0_j 		= 0.029; %constant Ca influx (EC)
C_m 		= 25.8;
G_tot		= 6927; %e-6; %changed
vK_j 		= -80;
a1			= 53.3;
a2			= 53.3;
b			= -80.8;
c 			= -0.4; %-6.4; %changed value!!! 
m3b			= 1.32e-3;
m4b			= 0.3;
m3s			= -0.28;
m4s			= 0.389;
G_R			= 955;
v_rest		= -31.1;
k_j			= 0.1;

%Constants by Hannah's model

VR_ERcyt = 0.185;   %DL        ratio of ER volume to cytosolic volume
J_max    = 2880;    %uMs-1      max IP3 production rate
K_I      = 0.03;    %uM         IP3 dissociation constant
K_act   = 0.17;     %uM         dissociation constant for Ca2+
k_on     = 2;       %uMs-1         Rate of Ca2+ binding to IP3 receptors
k_inh    = 0.1;     %uM         dissociation constant for Ca2+ binding
V_max    = 20;      %uMs-1        max pumping constant for ATP dependent calcium pump
k_pump   = 0.24;    %uM         dissociation constant for the pump
r_h      = 4.8;     %uMs-1       max rate of IP3 production 
k_deg    = 1.25;    % s-1       rate constant for IP3 degrdation
K_G      = 8.82;    %uM         G-protein dissociation constant
sig      = 1.235*10e-3; %uM  ratio of activities of bound and unbound receptors
ck_min   = 0.1;     %uM          minimum Ca2+ concentration for EET production
V_eet    = 72;      %uM           rate constant for EET production
k_eet    = 7.2;     %uM          rate constant for EET degrdation
eet_shift =  2e-3;    %VuM-1
vh_4     = 14.5e-3;    %V
vh_5     = 8e-3;       %V       
vh_6     = -15e-3;     %V
Ca_3     = 0.4;     %uM
Ca_4     = 0.15;    %uM
psi_h    = 2.664;   %s-1
BK_end   = 40;      %DL Fink    L&E
K_ex     = 0.26;    %uM Fink    L&E
B_ex     = 11.35;   %uM Fink    L&E

if NVU == 1
    P_L      = 5.2; %  uM; (Farr)    concentration factor based on steady state calcium balance
else
    P_L      = 0.0804;  %uMs-1L&E      5.2uM (Farr)    concentration factor based on steady state calcium balance
end



global CASE J_PLC g_hat p_hat p_hatIP3

if CASE==0
g_hat 		= 0;
p_hat 		= 0;
p_hatIP3 	= 0;
elseif CASE==1
g_hat 		= 0.5;
p_hat 		= 0;
p_hatIP3 	= 0.05;
elseif CASE==2
g_hat 		= 0.5;
p_hat 		= 0.05;
p_hatIP3 	= 0.05;
elseif CASE==3
g_hat 		= 0;
p_hat 		= 0;
p_hatIP3 	= 0.05;
elseif CASE==4 % was CASE 5 before!
g_hat 		= 0.5;
p_hat 		= 0.05;
p_hatIP3 	= 0;
elseif CASE==5 % was CASE 6 before!
g_hat 		= 0.5;
p_hat 		= 0;
p_hatIP3 	= 0;
elseif CASE==6 % was CASE 7 before!
g_hat 		= 0;
p_hat 		= 0.05;
p_hatIP3 	= 0;
elseif CASE==7 % was CASE 8 before!
g_hat 		= 0;
p_hat 		= 0.05;
p_hatIP3 	= 0.05;
elseif CASE==8     % not really part of the set
g_hat 		= 50;
p_hat 		= 0.05;
p_hatIP3 	= 0.05;
end


%% Myosin crossbridge model
global C_Hillmann
% K2_c        = 0.5 * C_Hillmann;  % NO pathway: cGMP-dependend!
K3_c        = 0.4 * C_Hillmann;
K4_c        = 0.1 * C_Hillmann;
% K5_c        = 0.5 * C_Hillmann;  % NO pathway: cGMP-dependend!
K7_c        = 0.1 * C_Hillmann;
gam_cross   = 17 * C_Hillmann;

%% Koningsberger

% nu_r        = 100*inv(0.0075);
% sigp0_r     = 0.0191*inv(0.0075);
% kp_r        = 0.15;
% r0_r        = 20;
% siga0_r     = 1.8e5;
% ka_r        = 0.0006;
% ra_r        = 12;
% rb_r        = 15;
% hb_r        = 3;
% P_r         = 30*inv(0.0075);


%% Kelvin Voigt

P_r         = 4000;  % Pa
rb_r        = 20e-6; % m
h0_r        = 3e-6;  % m 
R0pas_r     = 20e-6;
R0act_r     = 12e-6;
Eact_r      = 233e3;
Epas_r      = 66e3;
nu_r        = 1e4;
%% NO pathway

% NE********************
LArg		= 100;			% Nochmal genauer ausrechnen! --> Ping paper! - wird spaeter variabel!
%F           = 96500;        % [-] ; Faraday's constant (see above)
v_spine     = 8e-8;         % [fL] ; the volume of the neuronal dendritic spine
k_ex        = 1600;         % [s^{-1}] ; the decay rate constant of internal calcium concentration
Ca_rest     = 0.1; 			% [\muM] ; the resting calcium concentration (in Comerford+David2008: 2.830 mM; in Santucci2008P: 0.1 \muM)
lambda      = 20;           % [-] ; the buffer capacity

V_nNOS      = 1.435;        % [-] ; NO production - nNOS concentration ratio (Chen+Popel2007)

K_actNOS    = 9.27e-2;      % [microM]

%tau_ni      = 0.01;       	% [s]; time for NO to diffuse the distance between the neuron and SMC ; estimation
%tau_ji      = 0.01;        % [s]; time for NO to diffuse the distance between the EC and SMC ; estimation
D_NO 		= 3300;			% [\muM^2 s^{-1}] ; NO diffusion coefficient 
% dist_ni		= 50;			% [\mu m] ; estimation

dist_na     = 25;
dist_ai     = 25;
dist_ij		= 3.75;			% [\mu m] ; Kavdia 2002

% tau_ni      = dist_ni^2/(2*D_NO);  % Einstein-Smoluchowski equation, Lancaster1997
tau_na      = dist_na^2/(2*D_NO);  % Einstein-Smoluchowski equation, Lancaster1997
tau_ai      = dist_ai^2/(2*D_NO);  % Einstein-Smoluchowski equation, Lancaster1997
tau_ij      = dist_ij^2/(2*D_NO);  % Einstein-Smoluchowski equation, Lancaster1997

k_O2        = 9.6e-6;       % [microM^{-2} s^{-1}] 
On          = 200;         	% [microM] ; the tissue O2 concentration in the neuron
v_n         = -40;          % [mV] ; the neuronal membrane potential , assumed to be approx constant in this model
G_M         = 46;        	% [pS] ; the conductance of the NMDA channel to Ca2+ compaired  
P_Ca_P_M    = 3.6;         	% [-] ; the relative conductance of the NMDA channel to Ca2+ compared to monovalent ions
Ca_ex       = 2e3;          % [microM] ; the external calcium concentration (in Comerford+David2008: 1.5 mM!)
M           = 1.3e5;        % [microM] ; the concentration of monovalent ions in the neuron
R           = 8.314;      	% [Jmol^{-1}K^{-1}]
T           = 310.65;     	% [K] ; temperature
betA        = 650 ;       	% [uM] ; changed it!
betB        = 2800 ;        % [uM] ; changed it!
Q1          = 1.9e5;      	% [-]
Q2          = 2.1e5;      	% [-]
Q3          = 0.4e5;      	% [-]
Q4          = 0.26e5;      	% [-]
CaM_thresh = Ca_rest/( Ca_rest*(Q1 + 2*Q1*Q2*Ca_rest + 3*Q1*Q2*Q3*Ca_rest^2 + 4*Q1*Q2*Q3*Q4*Ca_rest^3)/ (1 + Q1*Ca_rest + Q1*Q2*Ca_rest^2 + Q1*Q2*Q3*Ca_rest^3 + Q1*Q2*Q3*Q4*Ca_rest^4) );
% CaM_thresh  = 2.7764e-5;       %  [mM] 

% EC**********************
Oj          = 200;         	% [\muM]; the O2 concentration in the EC
 K_dis       = 9e-2;    		% [\muM s^{-1}]          % = 0.09 [\mu M s^{-1}]
K_eNOS      = 4.5e-1;    	% [\muM]            % = 0.45 [\mu M] ; Michaelis constant for dx(eNOS_act)
 mu2         = 0.0167;       % [s^{-1}] ; the rate constant at which the eNOS is deactivated 
g_max       = 0.3;    		% [microM s^{-1}], maximal wss activation - fitted to Kavdia2003; in Comerford2008: 0.06, in Hannahs thesis 17.6, because she mixed up qmax and gmax in Comerford2008 
alp         = 2;            % [-] (in Wiesner1997: 3)
W_0         = 1.4;        	% [Pa^{-1}]
delt_wss    = 2.86 ;        % [Pa] ; the membrane shear modulus


% SMC*********************
k_dno       = 0.01;         % [s^{-1}]  (Hannahs code)
k1          = 2e3 ;    		% [\muM^{-1}s^{-1}] == 2000 muM^{-1}s^{-1}; a rate constant
k2          = 0.1;          % [s^{-1}]; a rate constant
k3          = 3;        	% [\muM^{-1}s^{-1}] == 3 muM^{-1}s^{-1}; a rate constant
k_1         = 100;          % [s^{-1}]; a rate constant


% m =2; % cGMP influence (0 - lowest influence)
global m

if m==0
    V_max_sGC = 1.26;  % \muM s^{-1}; the maximum cGMP production rate
    k_pde = 0.0695; % s^{-1}
	C_4 = 0.4; % [s^{-1}] (note: the changing units are correct!)
elseif m==1
    V_max_sGC = 1.09;  % \muM s^{-1};
	k_pde = 0.032;% s^{-1}
	C_4 = 0.098; % [s^{-1} microM^{-1}] (note: the changing units are correct!)
elseif m==2
    V_max_sGC = 0.8520;  % \muM s^{-1};
	k_pde = 0.0195;% s^{-1}
	C_4 = 0.011; % [s^{-1} microM^{-2}] (note: the changing units are correct!)
end
%%

K_m_pde = 2;           		% [microM]
R_Kfit = 30.8; % 55;                  % [mV]
V_cGMP = 66.9; % 68 ;               % [mV]
V_NO = 100;                 % [mV]
V_b = 283.7; % 215; (Hannah)                 % [mV]
K_m_cGMP = 0.55;       		% [microM]
K_m_NO = 0.2;     		    % [microM]
k_mlcp_b = 0.0086;          % [s^{-1}]
k_mlcp_c = 0.0327;          % [s^{-1}]
K_m_mlcp = 5.5;        		% [microM]
k_mlck = 1180;

v_Ca3		= -27; % cGMP and NO dependent

c_wi = 0; % translation factor for Ca dependence of KCa channel activation sigmoidal [microM] 
bet_i= 0.13; % translation factor for membrane potential dependence of KCa channel activation sigmoidal [microM^2] 

K_mArg = 1.5; % in EC and neuron, maybe we want to make it cell-specific like O2! --> Chen2006 & 2007
K_mO2_j = 7.7;
K_mO2_n = 243; % 140-145 uM Chen2007

    V_maxNOS    = 25e-3;   		% [\muM] 1.683;%0.025; %2.5e-8;%2.5e-8; %1.683e-4; %1.324; %2.5e-5; % 0.025;% 2.5e-8; %2.4925e-8;  %0.025; %2.0265e-9; % 1.054e-7;  % 1.683;[mM s^{-1}] (Hayashi1999)
    %V_eNOS     = 0.24;         % [-] ; NO production - eNOS concentration ratio (Chen+Popel2006)


V_NOj_max = 1.22;% 0.24; % 1.22; % s^-1, maximum catalytic rate of NO production (Chen2006) - obtained from fig 6 & equ 17 & 18 - we might want to change this value! - Liu: 0.585 micromol / min / (mg of protein)
V_NOn_max = 4.22; % s^-1, maximum catalytic rate of NO production (Chen2006) - obtained from fig 6 & equ 17 & 18 - we might want to change this value! - Liu: 0.585 micromol / min / (mg of protein)
