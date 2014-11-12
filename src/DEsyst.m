function [dy] = DEsyst (time,state)
% Below all conservation equation are calculated. 
% Note that "getRef" is the function that contains the input stimulus of 
% the model. 

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
dy(ind.N_Na_s  ) = -  k_C *0* getRef(time,'ft') - dy(ind.N_Na_k);                         % uMm s-1
dy(ind.N_K_s   ) = k_C *0*getRef(time,'ft') - dy(ind.N_K_k) + AC(flu.J_BK_k); % uMm s-1 ;     

dy(ind.N_HCO3_s) = - dy(ind.N_HCO3_k);                                                  % uMm s-1
dy(ind.K_p     ) = AC(flu.J_BK_k) / (VR_pa*state(ind.R_k)) + (SMC(flu.J_KIR_i))/(VR_ps);     % uM s-1
dy(ind.Ca_p     ) =  -AC(flu.J_TRPV_k)/(VR_pa)-Ca_dec*(AC(flu.Ca_p)-Ca_p0)  ;
%dy(ind.Ca_p     ) =  -AC(flu.J_TRPV_k)- SMC(flu.J_Ca)  - Ca_dec*(AC(flu.Ca_p)-Ca_p0) ;
%;
%
%
dy(ind.w_k     ) = AC(flu.phi_w) * (AC(flu.w_inf) - state(ind.w_k)); 
%dy(ind.z_k     ) = 1/(AC(flu.t_Ca)*AC(flu.Ca_p))*(AC(flu.z_inf)-state(ind.z_k)); % s-1
dy(ind.z_k     ) = 1/(2000*AC(flu.t_Ca))*(AC(flu.z_inf)-state(ind.z_k));


% %astr. DE for Ca2+ (Hannah)
%dy(ind.ck)=B_cyt*(AC(flu.J_ip3)-AC(flu.J_pump)+AC(flu.J_ERleak)); % FARR astrocytic calcium concentration
dy(ind.ck)=AC(flu.B_cyt)*(AC(flu.J_ip3)-AC(flu.J_pump)+AC(flu.J_ERleak)+1000*AC(flu.J_TRPV_k)); %L&E astrocytic calcium concentration
%
%dy(ind.sk)=-1/VR_ERcyt*dy(ind.ck); % ER calcium concentration
dy(ind.sk)=-1/VR_ERcyt*(dy(ind.ck)-AC(flu.J_TRPV_k)); % ER calcium concentration
%
dy(ind.hk)=k_on*(k_inh-(AC(flu.ck)+k_inh)*AC(flu.hk));  %the action of  the IP3 receptors that have not been inactivated by Calcium
dy(ind.ik)=r_h*AC(flu.G_pr)-k_deg*AC(flu.ik);  %IP3 concentration 
%dy(ind.eetk)=V_eet*(AC(flu.ck)-ck_min)-k_eet*AC(flu.eetk);  % FARR EET concentration
if AC(flu.ck)> ck_min
    dy(ind.eetk)=V_eet*(AC(flu.ck)-ck_min)-k_eet*AC(flu.eetk);  % L&E EET concentration
else 
    dy(ind.eetk)= -k_eet*AC(flu.eetk);              %L&E
end

% Smooth muscle cell
dy(ind.Ca_i)    = (SMC(flu.Ca_coup_i) + SMC(flu.rho_i) * (SMC(flu.J_CICR_i) + SMC(flu.J_IP3_i) + SMC(flu.J_leak_i) - SMC(flu.J_SRuptake_i) - SMC(flu.J_extrusion_i)...
    - SMC(flu.J_VOCC_i) + SMC(flu.J_NaCa_i) + 0.1*SMC(flu.J_stretch_i)));
% - SMC(flu.J_Ca)
dy(ind.s_i) 	= - SMC(flu.J_CICR_i) - SMC(flu.J_leak_i) + SMC(flu.J_SRuptake_i) ;
dy(ind.v_i)     = SMC(flu.v_coup_i) + gam * (-SMC(flu.J_NaK_i) - SMC(flu.J_Cl_i) - 2*SMC(flu.J_VOCC_i) - SMC(flu.J_NaCa_i) - SMC(flu.J_K_i) - SMC(flu.J_stretch_i) - SMC(flu.J_KIR_i));
%-SMC(flu.J_Ca)
dy(ind.w_i) 	= lab * (SMC(flu.Kactivation_i) - state(ind.w_i));
dy(ind.I_i) 	= SMC(flu.IP3_coup_i) - SMC(flu.J_degrad_i)  ; 

dy(ind.K_i)     = - SMC(flu.J_KIR_i) - SMC(flu.J_K_i) + SMC(flu.J_NaK_i);                                            % uM s-1

% Endothelium cell
dy(ind.Ca_j)	= EC(flu.Ca_coup_j) + EC(flu.rho_j) * (EC(flu.J_IP3_j) - EC(flu.J_ERuptake_j) + EC(flu.J_CICR_j) - EC(flu.J_extrusion_j) + EC(flu.J_leak_j)...
    + EC(flu.J_cation_j) + EC(flu.J_0_j) + EC(flu.J_stretch_j));
dy(ind.s_j) 	= EC(flu.J_ERuptake_j) - EC(flu.J_CICR_j) - EC(flu.J_leak_j) ;
dy(ind.v_j) 	=  - 1/C_m * ( EC(flu.J_K_j)	+ EC(flu.J_R_j) ) + EC(flu.v_coup_j);	
%dy(ind.I_j) 	= EC(flu.IP3_coup_j) + J_PLC - EC(flu.J_degrad_j)  ;
dy(ind.I_j) 	= 0;
% Myosin crossbridge model


dy(ind.Mp)      = K4_c*state(ind.AMp) + SMC(flu.K1_c) *SMC(flu.M) - (K2_c + K3_c)*state(ind.Mp);
dy(ind.AMp)     = K3_c*state(ind.Mp) + SMC(flu.K6_c) *state(ind.AM) - (K4_c + K5_c)*state(ind.AMp);  % K7_c was corrected to K4_c
dy(ind.AM)      = K5_c*state(ind.AMp) - (K7_c + SMC(flu.K6_c) )*state(ind.AM);

% Radius change

%Kevin Voigt

F_r=state(ind.AMp) + state(ind.AM);

E_r = Epas_r + F_r*(Eact_r -Epas_r);
R0_r= R0pas_r + F_r*R0pas_r*(0.6 - 1);


%dy(ind.R)= R0pas_r/nu_r *(state(ind.R)*P_r/SMC(flu.h_r) - E_r * ((state(ind.R) - R0_r)/R0_r));
dy(ind.R)= R0pas_r/nu_r *(AC(rad.R2)*P_r/SMC(flu.h_r) - E_r * ((AC(rad.R2) - R0_r)/R0_r));

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

end