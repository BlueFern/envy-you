function [NE,AC,SMC,EC] = all_fluxes (t,state)
%% load the constants for the fluxes and pointers:
    all_indices();
    all_constants();
    global stretch_ch only_Koenig NVU
%% Calculate the fluxes for the Astrocyte (AC)

% Below all the additional equations are calculated and stores in AC, SMC
% and EC
NE=[];

AC(flu.R_s)    = R_tot - state(ind.R_k);                               % m

AC(flu.N_Cl_s) = state(ind.N_Na_s) + state(ind.N_K_s) - state(ind.N_HCO3_s);   % uMm

AC(flu.Na_k  ) = state(ind.N_Na_k)/state(ind.R_k);  % uM
AC(flu.K_k   ) = state(ind.N_K_k)/state(ind.R_k);   % uM
AC(flu.HCO3_k) = state(ind.N_HCO3_k)/state(ind.R_k);% uM
AC(flu.Cl_k  ) = state(ind.N_Cl_k)/state(ind.R_k);  % uM
AC(flu.Na_s  ) = state(ind.N_Na_s)/AC(flu.R_s);     % uM
AC(flu.K_s   ) = state(ind.N_K_s) /AC(flu.R_s);     % uM
AC(flu.HCO3_s) = state(ind.N_HCO3_s)/AC(flu.R_s);   % uM
AC(flu.Cl_s  ) = AC(flu.N_Cl_s)/AC(flu.R_s);        % uM

AC(flu.ck)     = state(ind.ck);
AC(flu.sk)     = state(ind.sk);
AC(flu.hk)     = state(ind.h_k_star);
AC(flu.ik)     = state(ind.ik);
AC(flu.eetk)   = state(ind.eetk);

AC(flu.E_Na_k ) = (R_gas * Temp) / (z_Na * Farad) * log(AC(flu.Na_s)/AC(flu.Na_k)); % V
AC(flu.E_K_k )  = (R_gas * Temp) / (z_K  * Farad) * log(AC(flu.K_s )/AC(flu.K_k )); % V
AC(flu.E_Cl_k ) = (R_gas * Temp) / (z_Cl * Farad) * log(AC(flu.Cl_s)/AC(flu.Cl_k)); % V
AC(flu.E_NBC_k )= (R_gas * Temp) / (z_NBC* Farad) * ...
              log((AC(flu.Na_s)*AC(flu.HCO3_s)^2)/(AC(flu.Na_k)*AC(flu.HCO3_k)^2));     % V 
AC(flu.E_BK_k)  = (R_gas * Temp) / (z_K  * Farad) * log((state(ind.K_p_star)*K_0)/AC(flu.K_k));% V


AC(flu.J_NaK_k  ) = J_NaK_max * Hill(AC(flu.Na_k), K_Na_k, 1.5) * ...
                Hill(AC(flu.K_s),K_K_s,1);              % uMm s-1 

AC(flu.v_k )   = ( g_Na_k  * AC(flu.E_Na_k )...
             + g_K_k   * AC(flu.E_K_k  )...
             + g_Cl_k  * AC(flu.E_Cl_k )...
             + g_NBC_k * AC(flu.E_NBC_k)...
             - AC(flu.J_NaK_k)*Farad/unitcon...
             + g_BK_k *state(ind.w_k) * AC(flu.E_BK_k)          )...
            /(g_Na_k + g_K_k + g_Cl_k + g_NBC_k + g_BK_k*state(ind.w_k));  % V

       
AC(flu.J_KCC1_k ) = getRef(t,'fluxft')*...
                (R_gas * Temp * g_KCC1_k) / (Farad^2) * log((AC(flu.K_s)...
                *AC(flu.Cl_s))/(AC(flu.K_k)*AC(flu.Cl_k)))*unitcon;                  %uMm s-1

AC(flu.J_NBC_k  ) = g_NBC_k / Farad * (AC(flu.v_k) - AC(flu.E_NBC_k))*unitcon;       %uMm s-1
AC(flu.J_NKCC1_k) = getRef(t,'fluxft')*...
                (g_NKCC1_k * R_gas * Temp) / (Farad^2) ...
                * log((AC(flu.K_s) * AC(flu.Na_s) * AC(flu.Cl_s)^2)...
                     /(AC(flu.K_k) * AC(flu.Na_k) * AC(flu.Cl_k)^2))*unitcon;        %uMm s-1            
AC(flu.J_Na_k  ) = g_Na_k / Farad * (AC(flu.v_k) - AC(flu.E_Na_k))*unitcon;          %uMm s-1
AC(flu.J_K_k   ) = g_K_k  / Farad * (AC(flu.v_k) - AC(flu.E_K_k ))*unitcon;          %uMm s-1
AC(flu.J_BK_k)   = g_BK_k / Farad * state(ind.w_k)*(AC(flu.v_k)-AC(flu.E_BK_k))*unitcon; %uMm s-1

     %[-]
%Hannah/NVU
  AC(flu.vh_3)     = vh_5/2*tanh((AC(flu.ck)-Ca_3)/Ca_4);
 
 if NVU == 1
     AC(flu.w_inf)    = 0.5*(1+tanh((AC(flu.v_k)+v_6)/(v_4)));  %NVU
     AC(flu.phi_w)    = psi_w*cosh((AC(flu.v_k)+v_6)/(2*v_4));      %NVU
 elseif NVU == 2
     AC(flu.w_inf) =0.5*(1+tanh(((AC(flu.v_k)+(eet_shift*AC(flu.eetk))-AC(flu.vh_3)))/(vh_4))); %NVU+EET&CA2+
     AC(flu.phi_w)    = psi_h*cosh((AC(flu.v_k)-AC(flu.vh_3))/(2*vh_4));  %NVU+Ca2+
 elseif NVU == 3
     AC(flu.w_inf)    = 0.5*(1+tanh((AC(flu.v_k)+eet_shift*AC(flu.eetk)+v_6)/vh_4)); %NVU+EET
     AC(flu.phi_w)    = psi_w*cosh((AC(flu.v_k)+v_6)/(2*v_4));      %NVU
 else
     AC(flu.w_inf) =0.5*(1+tanh(((AC(flu.v_k))-AC(flu.vh_3))/(vh_4)));  %NVU+Ca2+
     AC(flu.phi_w)    = psi_h*cosh((AC(flu.v_k)-AC(flu.vh_3))/(2*vh_4));  %NVU+Ca2+
 end


%astrocyte fluxes by Hannah

AC(flu.J_ERleak)=P_L*(1-AC(flu.ck)/AC(flu.sk)); %calcium leak flux from ER to the cytosol
AC(flu.J_pump)=V_max*AC(flu.ck)^2/(AC(flu.ck)^2+k_pump^2); %ATP dependant calcium flux from cytoplasm to ER
AC(flu.J_ip3)=J_max*((AC(flu.ik)/(AC(flu.ik)+K_I))*(AC(flu.ck)/(AC(flu.ck)+K_act))*AC(flu.hk))^3*(1-AC(flu.ck)/AC(flu.sk)); %calcium flux from ER to cytosolic by IP3 receptors!
AC(flu.G_pr)=(getRef(t,'rho')+sig)/(K_G+getRef(t,'rho')+sig);
if NVU ==1
    AC(flu.B_cyt)= 0.0244; %NVU
else
    AC(flu.B_cyt)=(1+BK_end+(K_ex*B_ex)/(K_ex+AC(flu.ck))^2)^-1;  %Loes and Evert determined extra equation
end

%% SMC

SMC(flu.M)                   = 1 - state(ind.Mp_star) - state(ind.AM_star) - state(ind.AMp_star);                         
SMC(flu.E_K_i)              = (R_gas * Temp) / (z_K  * Farad)*unitcon*log((state(ind.K_p_star)*K_0)/(state(ind.K_i_star)*K_0));
% SMC(flu.h_r)                 =  -state(ind.R) + sqrt(state(ind.R)^2 + 2*rb_r*h0_r + h0_r^2);
SMC(flu.h_r)                 = 0.1* state(ind.R);

SMC(flu.v_coup_i)            = - g_hat * ( (state(ind.v_i_star)*v_0) - (state(ind.v_j_star)*v_0) );   
SMC(flu.Ca_coup_i)           = - p_hat * ( (state(ind.Ca_i_star)*Ca_0) - (state(ind.Ca_j_star)*Ca_0) );
SMC(flu.IP3_coup_i)          = - p_hatIP3 * ( (state(ind.I_i_star)*I_0) - (state(ind.I_j_star)*I_0) );
SMC(flu.rho_i)               = 1;%( K_d + state(ind.Ca_i ))^2 / ( ( K_d + (state(ind.Ca_i_star)*Ca_0) )^2 + ( K_d * B_T ) );
SMC(flu.J_IP3_i)             = Fmax_i * ( (state(ind.I_i_star)*I_0)^2 ) / ( Kr_i^2 + (state(ind.I_i_star)*I_0)^2 );
SMC(flu.J_SRuptake_i)        = B_i * ( (state(ind.Ca_i_star)*Ca_0)^2 ) / ( (state(ind.Ca_i_star)*Ca_0)^2 + cb_i^2 );
SMC(flu.J_CICR_i)            = C_i * ( (state(ind.s_i_star)*Ca_0)^2 ) / ( sc_i^2 + (state(ind.s_i_star)*Ca_0)^2 ) * ( (state(ind.Ca_i_star)*Ca_0)^4 ) / ( cc_i^4 + (state(ind.Ca_i_star)*Ca_0)^4 );
SMC(flu.J_extrusion_i)       = D_i * (state(ind.Ca_i_star)*Ca_0) * ( 1 + ( (state(ind.v_i_star)*v_0) - vd_i ) / ( Rd_i ) );
SMC(flu.J_leak_i)            = L_i * (state(ind.s_i_star)*Ca_0);
SMC(flu.J_VOCC_i)            = G_Ca * ( (state(ind.v_i_star)*v_0) - v_Ca1) / ( 1 + exp( - ( (state(ind.v_i_star)*v_0) - v_Ca2 ) / ( R_Ca ) ) );
SMC(flu.J_NaCa_i)            = G_NaCa * (state(ind.Ca_i_star)*Ca_0)* ( (state(ind.v_i_star)*v_0) - v_NaCa ) / ( (state(ind.Ca_i_star)*Ca_0) + c_NaCa );
SMC(flu.J_NaK_i)             = F_NaK;
SMC(flu.J_Cl_i)              = G_Cl * ( (state(ind.v_i_star)*v_0) - v_Cl );
SMC(flu.J_K_i)              = G_K * state(ind.w_i) * ( (state(ind.v_i_star)*v_0) - vK_i );
SMC(flu.Kactivation_i)      = ( (state(ind.Ca_i_star)*Ca_0) + c_w )^2 / ( ((state(ind.Ca_i_star)*Ca_0) + c_w)^2 + bet*exp(-((state(ind.v_i_star)*v_0) - v_Ca3)/R_K) );
SMC(flu.J_degrad_i)         = k_i * (state(ind.I_i_star)*I_0);

if strcmp(stretch_ch,'ON') == 1
 SMC(flu.J_stretch_i) = G_stretch/(1+exp(-alpha1*(P_str*state(ind.R)/SMC(flu.h_r) - sig0))) * ((state(ind.v_i_star)*v_0) - Esac);
elseif strcmp(stretch_ch,'OFF') == 1
   SMC(flu.J_stretch_i)     = 0;
end

if strcmp(only_Koenig,'OFF') == 1
   SMC(flu.v_KIR_i)    = z_1 * (state(ind.K_p_star)*K_0)/unitcon + z_2;                                               % mV
   SMC(flu.G_KIR_i)    = exp( z_5 * (state(ind.v_i_star)*v_0) + z_3 * (state(ind.K_p_star)*K_0)/unitcon + z_4 ); %exp( z_5 * (state(ind.v_i_star)*v_0) + z_3 * (state(ind.K_p_star)*K_0)/unitcon + z_4 );                     % pS pF-1 =s-1
   SMC(flu.J_KIR_i)    = F_il/gam * SMC(flu.G_KIR_i)*((state(ind.v_i_star)*v_0)-SMC(flu.v_KIR_i));                                % mV s-1
elseif strcmp(only_Koenig,'ON') == 1
   SMC(flu.v_KIR_i)    = 0;
   SMC(flu.G_KIR_i)    = 0;
   SMC(flu.J_KIR_i)    = 0;
end

SMC(flu.K1_c)       = gam_cross*(state(ind.Ca_i_star)*Ca_0)^3;
SMC(flu.K6_c)       = SMC(flu.K1_c);

%% EC

EC(flu.v_coup_j)             = - g_hat * ( (state(ind.v_j_star)*v_0) - (state(ind.v_i_star)*v_0) );  
EC(flu.Ca_coup_j)            = - p_hat * ( (state(ind.Ca_j_star)*Ca_0) - (state(ind.Ca_i_star)*Ca_0) );
EC(flu.IP3_coup_j)           = - p_hatIP3 * ( (state(ind.I_j_star)*I_0) - (state(ind.I_i_star)*I_0) );
EC(flu.rho_j)               = 1;%( K_d + (state(ind.Ca_j_star)*Ca_0) )^2 / ( ( K_d + (state(ind.Ca_j_star)*Ca_0) )^2 + ( K_d * B_T ) );
EC(flu.J_0_j)               = J0_j;
EC(flu.J_IP3_j)             = Fmax_j * ( (state(ind.I_j_star)*I_0)^2 ) / ( Kr_j^2 + (state(ind.I_j_star)*I_0)^2 );
EC(flu.J_ERuptake_j)        = B_j * ( (state(ind.Ca_j_star)*Ca_0)^2 ) / ( (state(ind.Ca_j_star)*Ca_0)^2 + cb_j^2 );
EC(flu.J_CICR_j)            = C_j *  ( (state(ind.s_j_star)*Ca_0)^2 ) / ( sc_j^2 + (state(ind.s_j_star)*Ca_0)^2 ) *  ( (state(ind.Ca_j_star)*Ca_0)^4 ) / ( cc_j^4 + (state(ind.Ca_j_star)*Ca_0)^4 );
EC(flu.J_extrusion_j)       = D_j * (state(ind.Ca_j_star)*Ca_0); 
EC(flu.J_leak_j)            = L_j * (state(ind.s_j_star)*Ca_0);
EC(flu.J_cation_j)          = G_cat * ( E_Ca - (state(ind.v_j_star)*v_0) )* 0.5 * ( 1 + tanh(( log10( (state(ind.Ca_j_star)*Ca_0) ) - m3cat )/( m4cat)) );
EC(flu.J_BKCa_j) 			= 0.4/2 * ( 1 + tanh( ( (  log10((state(ind.Ca_j_star)*Ca_0)) - c) * ( (state(ind.v_j_star)*v_0) - b ) - a1 ) / ( m3b*( (state(ind.v_j_star)*v_0) + a2 * ( log10( state(ind.Ca_j )) - c ) - b )^2 + m4b ) ) );
EC(flu.J_SKCa_j) 			= 0.6/2 * ( 1 + tanh( ( log10((state(ind.Ca_j_star)*Ca_0)) - m3s ) / ( m4s ) ) );
EC(flu.J_K_j)               = G_tot * ( (state(ind.v_j_star)*v_0) - vK_j ) * ( EC(flu.J_BKCa_j) + EC(flu.J_SKCa_j) );
EC(flu.J_R_j)               = G_R * ( (state(ind.v_j_star)*v_0) - v_rest);
EC(flu.J_degrad_j)          = k_j * (state(ind.I_j_star)*I_0);

if strcmp(stretch_ch,'ON') == 1
   EC(flu.J_stretch_j)      = G_stretch/(1+exp(-alpha1*(P_str*state(ind.R)/SMC(flu.h_r) - sig0))) * ((state(ind.v_j_star)*v_0) - Esac);
elseif strcmp(stretch_ch,'OFF') == 1
   EC(flu.J_stretch_j)      =  0;
end
end


%    A function that corrects for negative concentrations and sets them to 1e-180
%    Note that, when the system is running correctly this function is not
%    used.


function out = negCheck(input,wx)
    if (input / wx) > 0
        out = input/wx;
    else
        out = 1e-180;
    end
end











