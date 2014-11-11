function [NE,AC,SMC,EC] = all_fluxes (t,state)
%% load the constants for the fluxes and pointers:
    all_indices();
    all_constants();
    global stretch_ch only_Koenig NVU
%% Calculate the fluxes for the Astrocyte (AC)

% Below all the additional equations are calculated and stores in AC, SMC
% and EC
NE=[];

% n=1;
% if n*20<t<n*20+5;
%     AC(rad.R2)=30*1e-6;
%     n=n+1;
% else 
%     AC(rad.R2)=0;
% end

AC(rad.R2)=(20+50*getRef(t,'rho'))*1e-6;
AC(flu.R_s)    = R_tot - state(ind.R_k);                               % m    
AC(flu.N_Cl_s) = state(ind.N_Na_s) + state(ind.N_K_s) - state(ind.N_HCO3_s);   % uMm

% AC(flu.Na_k  ) = negCheck(state(ind.N_Na_k)  ,6.1e-8);  % uM
% AC(flu.K_k   ) = negCheck(state(ind.N_K_k)   ,6.1e-8);  % uM
% AC(flu.HCO3_k) = negCheck(state(ind.N_HCO3_k),6.1e-8);  % uM
% AC(flu.Cl_k  ) = negCheck(state(ind.N_Cl_k)  ,6.1e-8);  % uM
% AC(flu.Na_s  ) = negCheck(state(ind.N_Na_s)  ,2.69e-8);     % uM
% AC(flu.K_s   ) = negCheck(state(ind.N_K_s)   ,2.69e-8);     % uM
% AC(flu.HCO3_s) = negCheck(state(ind.N_HCO3_s),2.69e-8);     % uM
% AC(flu.Cl_s  ) = negCheck(AC(flu.N_Cl_s)     ,2.69e-8);     % uM

AC(flu.Na_k  ) = negCheck(state(ind.N_Na_k)  ,state(ind.R_k));  % uM
AC(flu.K_k   ) = negCheck(state(ind.N_K_k)   ,state(ind.R_k));  % uM
AC(flu.HCO3_k) = negCheck(state(ind.N_HCO3_k),state(ind.R_k));  % uM
AC(flu.Cl_k  ) = negCheck(state(ind.N_Cl_k)  ,state(ind.R_k));  % uM
AC(flu.Na_s  ) = negCheck(state(ind.N_Na_s)  ,AC(flu.R_s));     % uM
AC(flu.K_s   ) = negCheck(state(ind.N_K_s)   ,AC(flu.R_s));     % uM
AC(flu.HCO3_s) = negCheck(state(ind.N_HCO3_s),AC(flu.R_s));     % uM
AC(flu.Cl_s  ) = negCheck(AC(flu.N_Cl_s)     ,AC(flu.R_s));     % uM

AC(flu.ck)     = state(ind.ck);
AC(flu.sk)     = state(ind.sk);
AC(flu.hk)     = state(ind.hk);
AC(flu.ik)     = state(ind.ik);
AC(flu.eetk)   = state(ind.eetk);

AC(flu.Ca_p)   = state(ind.Ca_p);

AC(flu.E_Na_k ) = (R_gas * Temp) / (z_Na * Farad) * log(AC(flu.Na_s)/AC(flu.Na_k)); % V
AC(flu.E_K_k )  = (R_gas * Temp) / (z_K  * Farad) * log(AC(flu.K_s )/AC(flu.K_k )); % V
AC(flu.E_Cl_k ) = (R_gas * Temp) / (z_Cl * Farad) * log(AC(flu.Cl_s)/AC(flu.Cl_k)); % V
AC(flu.E_NBC_k )= (R_gas * Temp) / (z_NBC* Farad) * ...
              log((AC(flu.Na_s)*AC(flu.HCO3_s)^2)/(AC(flu.Na_k)*AC(flu.HCO3_k)^2));     % V 
AC(flu.E_BK_k)  = (R_gas * Temp) / (z_K  * Farad) * log(state(ind.K_p)/AC(flu.K_k));% V


AC(flu.J_NaK_k  ) = J_NaK_max * Hill(AC(flu.Na_k), K_Na_k, 1.5) * ...
                Hill(AC(flu.K_s),K_K_s,1);              % uMm s-1 
%Flux for BM
%AC(flu.E_TRPV_k) = 0.5;
AC(flu.E_TRPV_k) = (R_gas * Temp) / (z_Ca * Farad) * log(AC(flu.Ca_p)/AC(flu.ck));

AC(flu.v_k )   = ( g_Na_k  * AC(flu.E_Na_k )...
             + g_K_k   * AC(flu.E_K_k  )...
             + g_Cl_k  * AC(flu.E_Cl_k )...
             + g_NBC_k * AC(flu.E_NBC_k)...
             - AC(flu.J_NaK_k)*Farad/unitcon...
             + g_BK_k *state(ind.w_k) * AC(flu.E_BK_k)        ...
             + g_trpv*state(ind.z_k)*AC(flu.E_TRPV_k) )... 
             /(g_Na_k + g_K_k + g_Cl_k + g_NBC_k + g_BK_k*state(ind.w_k)+g_trpv*state(ind.z_k));  % V

         %    %  
         %  
       
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
AC(flu.G_pr)=(0*getRef(t,'rho')+sig)/(K_G+0*getRef(t,'rho')+sig);
%AC(flu.G_pr)=(getRef(t,'rho')+sig)/(K_G+getRef(t,'rho')+sig);
if NVU ==1
    AC(flu.B_cyt)= 0.0244; %NVU
else
    AC(flu.B_cyt)=(1+BK_end+(K_ex*B_ex)/(K_ex+AC(flu.ck))^2)^-1;  %Loes and Evert determined extra equation
end

%Fluxes for BM
AC(flu.t_Ca)=t_trpv/AC(flu.Ca_p);
AC(flu.H_Ca)= AC(flu.ck)/gam_cai+AC(flu.Ca_p)/gam_cae;
%AC(flu.z_inf) = (1/(1+exp(-((state(ind.R)-R0)/R0-eps_h)/k)))*(1/(1+AC(flu.H_Ca)))*(AC(flu.H_Ca)+tanh((AC(flu.v_k)-v_1i)/v_2i));
AC(flu.z_inf) = (1/(1+exp(-((AC(rad.R2)-R0)/R0-eps_h)/k)))*(1/(1+AC(flu.H_Ca)))*(AC(flu.H_Ca)+tanh((AC(flu.v_k)-v_1i)/v_2i));
AC(flu.J_TRPV_k)= (g_trpv*state(ind.z_k)*(AC(flu.v_k)-v_trpv)/(C_ast*gam));

%% SMC

SMC(flu.M)                   = 1 - state(ind.Mp) - state(ind.AM) - state(ind.AMp);                         
SMC(flu.E_K_i)              = (R_gas * Temp) / (z_K  * Farad)*unitcon*log(state(ind.K_p)/state(ind.K_i));
% SMC(flu.h_r)                 =  -state(ind.R) + sqrt(state(ind.R)^2 + 2*rb_r*h0_r + h0_r^2);
%SMC(flu.h_r)                 = 0.1* state(ind.R);
SMC(flu.h_r)                 = 0.1* AC(rad.R2);

SMC(flu.v_coup_i)            = - g_hat * ( state(ind.v_i) - state(ind.v_j) );   
SMC(flu.Ca_coup_i)           = - p_hat * ( state(ind.Ca_i) - state(ind.Ca_j) );
SMC(flu.IP3_coup_i)          = - p_hatIP3 * ( state(ind.I_i) - state(ind.I_j) );
SMC(flu.rho_i)              = 1;%( K_d + state(ind.Ca_i ))^2 / ( ( K_d + state(ind.Ca_i) )^2 + ( K_d * B_T ) );
SMC(flu.J_IP3_i)            = Fmax_i * ( state(ind.I_i)^2 ) / ( Kr_i^2 + state(ind.I_i)^2 );
SMC(flu.J_SRuptake_i)       = B_i * ( state(ind.Ca_i)^2 ) / ( state(ind.Ca_i)^2 + cb_i^2 );
SMC(flu.J_CICR_i)           = C_i *  ( state(ind.s_i)^2 ) / ( sc_i^2 + state(ind.s_i)^2 ) *  ( state(ind.Ca_i)^4 ) / ( cc_i^4 + state(ind.Ca_i)^4 );
SMC(flu.J_extrusion_i)      = D_i * state(ind.Ca_i) * ( 1 + ( state(ind.v_i) - vd_i ) / ( Rd_i ) );
SMC(flu.J_leak_i)           = L_i * state(ind.s_i);
SMC(flu.J_VOCC_i)           = G_Ca * ( state(ind.v_i) - v_Ca1) / ( 1 + exp( - ( state(ind.v_i) - v_Ca2 ) / ( R_Ca ) ) );
SMC(flu.J_NaCa_i)           = G_NaCa * state(ind.Ca_i)* ( state(ind.v_i) - v_NaCa ) / ( state(ind.Ca_i) + c_NaCa );
SMC(flu.J_NaK_i)            = F_NaK;
SMC(flu.J_Cl_i)             = G_Cl * ( state(ind.v_i) - v_Cl );
SMC(flu.J_K_i)              = G_K * state(ind.w_i) * ( state(ind.v_i) - vK_i );
SMC(flu.Kactivation_i)      = ( state(ind.Ca_i) + c_w )^2 / ( (state(ind.Ca_i) + c_w)^2 + bet*exp(-(state(ind.v_i) - v_Ca3)/R_K) );
SMC(flu.J_degrad_i)         = k_i * state(ind.I_i);

if strcmp(stretch_ch,'ON') == 1
   %SMC(flu.J_stretch_i)     = G_stretch/(1+exp(-alpha1*(P_str*state(ind.R)/SMC(flu.h_r) - sig0))) * (state(ind.v_i) - Esac);
   SMC(flu.J_stretch_i)     = G_stretch/(1+exp(-alpha1*(P_str*AC(rad.R2)/SMC(flu.h_r) - sig0))) * (state(ind.v_i) - Esac);
elseif strcmp(stretch_ch,'OFF') == 1
   SMC(flu.J_stretch_i)     = 0;
end

if strcmp(only_Koenig,'OFF') == 1
   SMC(flu.v_KIR_i)    = z_1 * state(ind.K_p)/unitcon + z_2;                                               % mV
   SMC(flu.G_KIR_i)    = exp( z_5 * state(ind.v_i) + z_3 * state(ind.K_p)/unitcon + z_4 ); %exp( z_5 * state(ind.v_i) + z_3 * state(ind.K_p)/unitcon + z_4 );                     % pS pF-1 =s-1
   SMC(flu.J_KIR_i)    = (F_il/gam * SMC(flu.G_KIR_i)*(state(ind.v_i)-SMC(flu.v_KIR_i)));                                % mV s-1
% elseif strcmp(only_Koenig,'ON') == 1
%    SMC(flu.v_KIR_i)    = 57*log10(state(ind.K_p))-130;
%    SMC(flu.G_KIR_i)    = 145e-9/A_ef_k*sqrt(state(ind.K_p));
%    SMC(flu.J_KIR_i)    = (F_il/gam * SMC(flu.G_KIR_i)*(state(ind.v_i)-SMC(flu.v_KIR_i)));
% end
elseif strcmp(only_Koenig,'ON') == 1
   SMC(flu.v_KIR_i)    = 0;
   SMC(flu.G_KIR_i)    = 0;
   SMC(flu.J_KIR_i)    = 0;
end

SMC(flu.K1_c)       = gam_cross*(state(ind.Ca_i))^3;
SMC(flu.K6_c)       = SMC(flu.K1_c);

%Fluxes for BM
% SMC(flu.m_inf) = 0.5*(1+tanh((state(ind.v_i)-v_1i)/v_2i));
% SMC(flu.J_Ca) =  (g_Ca*SMC(flu.m_inf)*(state(ind.v_i)-v_Ca)/(C_smc*gam));


%% EC

EC(flu.v_coup_j)             = - g_hat * ( state(ind.v_j) - state(ind.v_i) );  
EC(flu.Ca_coup_j)            = - p_hat * ( state(ind.Ca_j) - state(ind.Ca_i) );
EC(flu.IP3_coup_j)           = - p_hatIP3 * ( state(ind.I_j) - state(ind.I_i) );
EC(flu.rho_j)               = 1;%( K_d + state(ind.Ca_j) )^2 / ( ( K_d + state(ind.Ca_j) )^2 + ( K_d * B_T ) );
EC(flu.J_0_j)               = J0_j;
EC(flu.J_IP3_j)             = Fmax_j * ( state(ind.I_j)^2 ) / ( Kr_j^2 + state(ind.I_j)^2 );
EC(flu.J_ERuptake_j)        = B_j * ( state(ind.Ca_j)^2 ) / ( state(ind.Ca_j)^2 + cb_j^2 );
EC(flu.J_CICR_j)            = C_j *  ( state(ind.s_j)^2 ) / ( sc_j^2 + state(ind.s_j)^2 ) *  ( state(ind.Ca_j)^4 ) / ( cc_j^4 + state(ind.Ca_j)^4 );
EC(flu.J_extrusion_j)       = D_j * state(ind.Ca_j); 
EC(flu.J_leak_j)            = L_j * state(ind.s_j);
EC(flu.J_cation_j)          = G_cat * ( E_Ca - state(ind.v_j) )* 0.5 * ( 1 + tanh(( log10( state(ind.Ca_j) ) - m3cat )/( m4cat)) );
EC(flu.J_BKCa_j) 			= 0.4/2 * ( 1 + tanh( ( (  log10(state(ind.Ca_j)) - c) * ( state(ind.v_j) - b ) - a1 ) / ( m3b*( state(ind.v_j) + a2 * ( log10( state(ind.Ca_j )) - c ) - b )^2 + m4b ) ) );
EC(flu.J_SKCa_j) 			= 0.6/2 * ( 1 + tanh( ( log10(state(ind.Ca_j)) - m3s ) / ( m4s ) ) );
EC(flu.J_K_j)               = G_tot * ( state(ind.v_j) - vK_j ) * ( EC(flu.J_BKCa_j) + EC(flu.J_SKCa_j) );
EC(flu.J_R_j)               = G_R * ( state(ind.v_j) - v_rest);
EC(flu.J_degrad_j)          = k_j * state(ind.I_j);

if strcmp(stretch_ch,'ON') == 1
  %EC(flu.J_stretch_j)      = G_stretch/(1+exp(-alpha1*(P_str*state(ind.R)/SMC(flu.h_r) - sig0))) * (state(ind.v_j) - Esac);
  EC(flu.J_stretch_j)      = G_stretch/(1+exp(-alpha1*(P_str*AC(rad.R2)/SMC(flu.h_r) - sig0))) * (state(ind.v_j) - Esac);
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











