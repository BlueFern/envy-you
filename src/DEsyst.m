function [dy] = DEsyst (time,state)
% Below all conservation equation are calculated. 
% Note that "getRef" is the function that contains the input stimulus of 
% the model. 
global J_PLC Ca_switch

dy = zeros(size(state));
all_constants(); % All constants used in this model
all_indices() ; % All indices used in this model

% All additional equation are calculated using the state variables. They are stored in three matrices; AC, SMC and EC
% Note that EC is empty at the moment.

[NE,AC,SMC,EC] = all_fluxes(time, state); 

% Astrocyte
dy(ind.R_k     ) = L_p * (AC(flu.Na_k) + AC(flu.K_k) + AC(flu.Cl_k) + AC(flu.HCO3_k)...
             - AC(flu.Na_s) - AC(flu.K_s) - AC(flu.Cl_s) - AC(flu.HCO3_s) + X_k / state(ind.R_k));  % m s-1
dy(ind.N_Na_k  ) = -AC(flu.J_Na_k) - 3 * AC(flu.J_NaK_k) + AC(flu.J_NKCC1_k) + AC(flu.J_NBC_k );    % uMm s-1
dy(ind.N_K_k   ) = -AC(flu.J_K_k ) + 2 * AC(flu.J_NaK_k) + AC(flu.J_NKCC1_k) + AC(flu.J_KCC1_k)...
                -AC(flu.J_BK_k);                                                    % uMm s-1
dy(ind.N_HCO3_k) = 2 * AC(flu.J_NBC_k);                                                 % uMm s-1
dy(ind.N_Cl_k  ) = dy(ind.N_Na_k) + dy(ind.N_K_k) - dy(ind.N_HCO3_k);                           % uMm s-1, modified equation compared to the one of Ostby
dy(ind.N_Na_s  ) = - k_C * getRef(time,'ft') - dy(ind.N_Na_k);                         % uMm s-1
dy(ind.N_K_s   ) = k_C * getRef(time,'ft') - dy(ind.N_K_k) +AC(flu.J_BK_k) ;                          % uMm s-1
dy(ind.N_HCO3_s) = - dy(ind.N_HCO3_k);                                                  % uMm s-1
dy(ind.K_p     ) = AC(flu.J_BK_k) / (VR_pa*state(ind.R_k)) + (SMC(flu.J_KIR_i))/(VR_ps);     % uM s-1
dy(ind.w_k     ) = AC(flu.phi_w) * (AC(flu.w_inf) - state(ind.w_k));                            % s-1


%astr. DE for Ca2+ by Farr
dy(ind.Ca_k)=AC(flu.B_cyt)*(AC(flu.J_ip3)-AC(flu.J_pump)+AC(flu.J_ERleak)); %L&E astrocytic calcium concentration
dy(ind.s_k)=-1/VR_ERcyt*dy(ind.Ca_k); % ER calcium concentration
dy(ind.h_k)=k_on*(k_inh-(state(ind.Ca_k)+k_inh)*state(ind.h_k)); 
dy(ind.I_k)=r_h*AC(flu.G_pr)-k_deg*state(ind.I_k);  
if state(ind.Ca_k)> ck_min
    dy(ind.EET_k)=V_eet*(state(ind.Ca_k)-ck_min)-k_eet*state(ind.EET_k);  
else 
    dy(ind.EET_k)= -k_eet*state(ind.EET_k);              
end

% Smooth muscle cell
dy(ind.Ca_i)    = Ca_switch * (SMC(flu.Ca_coup_i) + SMC(flu.rho_i) * (SMC(flu.J_CICR_i) + SMC(flu.J_IP3_i) + SMC(flu.J_leak_i) - SMC(flu.J_SRuptake_i) - SMC(flu.J_extrusion_i)- SMC(flu.J_VOCC_i) + SMC(flu.J_NaCa_i) + 0.1*SMC(flu.J_stretch_i)));
dy(ind.s_i) 	= - SMC(flu.J_CICR_i) - SMC(flu.J_leak_i) + SMC(flu.J_SRuptake_i) ;
dy(ind.v_i)     = SMC(flu.v_coup_i) + gam * (-SMC(flu.J_NaK_i) - SMC(flu.J_Cl_i) - 2*SMC(flu.J_VOCC_i) - SMC(flu.J_NaCa_i) - SMC(flu.J_K_i) - SMC(flu.J_stretch_i) - SMC(flu.J_KIR_i));
dy(ind.w_i) 	= lab * (SMC(flu.Kactivation_i) - state(ind.w_i));
dy(ind.I_i) 	= SMC(flu.IP3_coup_i) - SMC(flu.J_degrad_i)  ; 

dy(ind.K_i)     = - SMC(flu.J_KIR_i) - SMC(flu.J_K_i) + SMC(flu.J_NaK_i);                                            % uM s-1

% Endothelium cell
dy(ind.Ca_j)	= EC(flu.Ca_coup_j) + EC(flu.rho_j) * (EC(flu.J_IP3_j) - EC(flu.J_ERuptake_j) + EC(flu.J_CICR_j) - EC(flu.J_extrusion_j) + EC(flu.J_leak_j)...
    + EC(flu.J_cation_j) + EC(flu.J_0_j) + EC(flu.J_stretch_j));
dy(ind.s_j) 	= EC(flu.J_ERuptake_j) - EC(flu.J_CICR_j) - EC(flu.J_leak_j) ;
dy(ind.v_j) 	=  - 1/C_m * ( EC(flu.J_K_j)	+ EC(flu.J_R_j) ) + EC(flu.v_coup_j);	
dy(ind.I_j) 	= EC(flu.IP3_coup_j) + J_PLC - EC(flu.J_degrad_j)  ;

% Myosin crossbridge model


%dy(ind.Mp)      = K4_c*state(ind.AMp) + SMC(flu.K1_c) *SMC(flu.M) - (K2_c + K3_c)*state(ind.Mp);
%dy(ind.AMp)     = K3_c*state(ind.Mp) + SMC(flu.K6_c) *state(ind.AM) - (K4_c + K5_c)*state(ind.AMp);  % K7_c was corrected to K4_c
%dy(ind.AM)      = K5_c*state(ind.AMp) - (K7_c + SMC(flu.K6_c) )*state(ind.AM);


% NO-Pathway:
dy(ind.Mp)      = K4_c*state(ind.AMp) + SMC(flu.K1_c)*SMC(flu.M) - (SMC(flu.K2_c) + K3_c)*state(ind.Mp);
dy(ind.AMp)     = K3_c*state(ind.Mp) + SMC(flu.K6_c)*state(ind.AM) - (K4_c + SMC(flu.K5_c))*state(ind.AMp); % K4 here was K7!!!!!!!!!!!!!!!!!!!
dy(ind.AM)      = SMC(flu.K5_c)*state(ind.AMp) - (K7_c + SMC(flu.K6_c))*state(ind.AM);


% Radius change

%Kevin Voigt

F_r=state(ind.AMp) + state(ind.AM);

E_r = Epas_r + F_r*(Eact_r -Epas_r);
R0_r= R0pas_r + F_r*R0pas_r*(0.6 - 1);


dy(ind.R)= R0pas_r/nu_r *(state(ind.R)*P_r/SMC(flu.h_r) - E_r * ((state(ind.R) - R0_r)/R0_r));

% if F_r1 <= 0.4
%     F_r = 0.4;
% else
%     F_r = F_r1;
% end

%dy(ind.R)       = 1/nu_r *( R0pas_r*state(ind.R)*P_r /SMC(flu.h_r)  - Epas_r * ( state(ind.R) - R0pas_r)...
%    - F_r/0.8 * ((Etot_r - Epas_r) * state(ind.R) + Epas_r*R0pas_r - Etot_r*R0act_r));


%Koningsberger
% siga_r = siga0_r *((state(AMp) + state(AM))/0.8) * exp(-ka_r *(state(R)+h_r-ra_r)^2);
% 
% if state(R) >= r0_r
%     sigp_r = sigp0_r * (exp(kp_r * (state(R)-r0_r)) -1);
% else
%     sigp_r = sigp0_r * kp_r * (1 - (state(R)^2/r0_r^2)^(-3/2));
% end
% 
% 
% dy(R)       = 1/nu_r * ( P_r*state(R)/h_r - sigp_r - siga_r);

%% NO pathway 
global nNOS_switch eNOS_switch
% NE           
       dy(ind.Ca_n)       = (((NE(flu.I_Ca))/(2*Farad*v_spine)-(k_ex*(state(ind.Ca_n)-Ca_rest)))/(1+lambda));                                                   % cytosolic [Ca2+] in the NE in \muM
       dy(ind.nNOS_act)   = nNOS_switch * (V_maxNOS*(NE(flu.CaM))/(K_actNOS+NE(flu.CaM))-mu2*state(ind.nNOS_act));                                    % activated nNOS in \muM
       dy(ind.NO_n)        = ((state(ind.nNOS_act) * V_NOn_max * (On/(K_mO2_n+On)) * (LArg_n/(K_mArg+LArg_n)) )) + ((state(ind.NO_k)-state(ind.NO_n))/tau_nk) - (k_O2*state(ind.NO_n)^2*On);   % NO concentration in the neuron ; (95)
       %maximal production - FIG 6:  
       dy(ind.nNOS_act_max)   = nNOS_switch * (V_maxNOS - mu2 * state(ind.nNOS_act_max));                                    % activated nNOS in \muM
       % absolutes Limit des Models: (unrealistisch)
       %        dy(ind.NOn_max)        = V_NOn_max*state(ind.nNOS_act_max) + ((state(ind.NOk_max)-state(ind.NOn_max))/tau_nk) - (k_O2*state(ind.NOn_max)^2*On);   % NO concentration in the neuron ; (95)
       dy(ind.NOn_max)        = (state(ind.nNOS_act) * V_NOn_max ) + ((state(ind.NOk_max)-state(ind.NOn_max))/tau_nk) - (k_O2*state(ind.NOn_max)^2*On);   % NO concentration in the neuron ; (95)

% AC
       dy(ind.NO_k)        = ((state(ind.NO_n)-state(ind.NO_k))/tau_nk) + ((state(ind.NO_i)-state(ind.NO_k))/tau_ki) - (k_O2*state(ind.NO_k)^2*Ok);
       % maximal production - FIG 6: 
       dy(ind.NOk_max)    = ((state(ind.NOn_max)-state(ind.NOk_max))/tau_nk) + ((state(ind.NOi_max)-state(ind.NOk_max))/tau_ki) - (k_O2*state(ind.NOk_max)^2*Ok);

% SMC              
       dy(ind.NO_i)        = ((state(ind.NO_k)-state(ind.NO_i))/tau_ki) + ((state(ind.NO_j)-state(ind.NO_i))/tau_ij) - (k_dno * state(ind.NO_i));                             % NO concentration in the SMC ; 
       % maximal production - FIG 6: 
       dy(ind.NOi_max)    = ((state(ind.NOk_max)-state(ind.NOi_max))/tau_ki) + ((state(ind.NOj_max)-state(ind.NOi_max))/tau_ij) - (k_dno * state(ind.NOi_max));                             % NO concentration in the SMC ; 
       
       dy(ind.E_b)        = -k1*state(ind.E_b)*state(ind.NO_i)+k_1*state(ind.E_6c) + SMC(flu.k4)*state(ind.E_5c);     
       dy(ind.E_6c)       = k1*state(ind.E_b)*state(ind.NO_i)-k_1*state(ind.E_6c)-k2*state(ind.E_6c)- k3*state(ind.E_6c)*state(ind.NO_i);
       dy(ind.E_5c)       = k3*state(ind.E_6c)*state(ind.NO_i)+k2*state(ind.E_6c) - SMC(flu.k4)*state(ind.E_5c);
       dy(ind.cGMP)       = V_max_sGC*state(ind.E_5c)-(k_pde*state(ind.cGMP)^2)/(K_m_pde+state(ind.cGMP)); % ^ 2 is correct! (Yang2005)
       
% EC
       gam_eNOS = 0.1;
       dy(ind.eNOS_act)   = eNOS_switch * (gam_eNOS * SMC(flu.Act_eNOS_Ca)  + (1-gam_eNOS) * SMC(flu.Act_eNOS_wss)- mu2*state(ind.eNOS_act)) ;          % (104)

       dy(ind.NO_j)        = (V_NOj_max * (state(ind.eNOS_act)) * (Oj/(K_mO2_j+Oj)) * (LArg_j/(K_mArg+LArg_j)) ) + ((state(ind.NO_i)-state(ind.NO_j))/tau_ij)   - k_O2*(state(ind.NO_j))^2*Oj - state(ind.NO_j)*4*3300/(25^2);
       % maximal production - FIG 6: 
       Act_eNOS_Ca_max = 0.09; % at infinitively high Ca_j concentration
       Act_eNOS_wss_max = 0.04; % at infinitively high wss concentration
       dy(ind.eNOS_act_max)   = eNOS_switch * (gam_eNOS * Act_eNOS_Ca_max + (1-gam_eNOS) * Act_eNOS_wss_max - mu2*state(ind.eNOS_act_max)) ;          % (104)
       % absolutes Limit des Models: (unrealistisch)
%        dy(ind.NOj_max)        = V_NOj_max*state(ind.eNOS_act_max) + ((state(ind.NOi_max)-state(ind.NOj_max))/tau_ij)   - k_O2*(state(ind.NOj_max))^2*Oj - state(ind.NOj_max)*4*3300/(25^2);
       dy(ind.NOj_max)        = V_NOj_max*state(ind.eNOS_act) + ((state(ind.NOi_max)-state(ind.NOj_max))/tau_ij)   - k_O2*(state(ind.NOj_max))^2*Oj - state(ind.NOj_max)*4*3300/(25^2);

% Yang2005 - simplified Hai & Murphy model
       dy(ind.M_Y)     = - k_mlck * state(ind.M_Y) + SMC(flu.kmlcp) * state(ind.Mp_Y);
       dy(ind.Mp_Y)    = + k_mlck * state(ind.M_Y) - SMC(flu.kmlcp) * state(ind.Mp_Y);
end