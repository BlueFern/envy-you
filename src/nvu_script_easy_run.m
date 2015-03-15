%% Demonstration script for 00-NVU model
% This script runs the NVU Model under pre-defined base conditions. The
% purpose of this script is to give a quick overview of the outcomes of the
% code for first time viewers. 

% Please note this is not an exhaustive instruction list. Please refer to
% nvu_script for more detail.

%% NVU Model Overview
% The NVU consists of a number of submodules, implemented as MATLAB
% classes, presently an astrocyte, a lumped SMC/EC model, and a model of
% the wall mechanics. 
%
% The parameters of each of these submodules are specified when the 
% modules are constructed. Note: if theses are to be changed please refer 
% to the nvu_script for more detail.
clear; clc; close all
%% Options 
% First the options need to be set for the ODE solver (currently |ode15s|) 
odeopts = odeset('RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 1, 'Vectorized', 1);

nv = NVU(Astrocyte(), ...
    WallMechanics(), ...
    SMCEC(), ...
    'odeopts', odeopts);
%% Run a basic simulation
nv.simulate()

%% Plots
% Lists of quantities that can be retrieved from the NVU model after
% simulation are given in the documentation pages of the individual model
% components.

figure(1) % Plot the Radius
plot(nv.T, 1e6 * nv.out('R'))
title('Radius (R)');   xlabel('time (s)'); ylabel('(\mum)'); grid on

figure(2) % Plot variables from the Astrocyte 
suptitle('Astrocyte')
subplot(3,2,1)
plot(nv.T, 0.001*nv.out('N_K_k')./(nv.out('R_k')))
title('Potassium in AC (K_k)');   xlabel('Time [s]'); ylabel('[K^+]   (mM)') ; grid on
subplot(3,2,2)
plot(nv.T, 0.001*nv.out('N_Na_k')./(nv.out('R_k')))
title('Sodium in AC (Na_k)');   xlabel('Time [s]'); ylabel('[Na^+]   (mM)'); grid on
subplot(3,2,3)
plot(nv.T, 0.001*nv.out('N_HCO3_k')./(nv.out('R_k')))
title('HCO_3 in AC (HCO3_k)');   xlabel('Time [s]'); ylabel('[HCO_3]   (mM)'); grid on
subplot(3,2,4)
plot(nv.T, 0.001*nv.out('N_Cl_k')./(nv.out('R_k')))
title('Chlorine in AC (Cl_k)');   xlabel('Time [s]'); ylabel('[Cl]   (mM)'); grid on
subplot(3,2,5)
plot(nv.T, nv.out('w_k'))
title('Open probability of the BK Channel in AC (N_K_k)');   xlabel('Time [s]'); ylabel('(-)'); grid on
subplot(3,2,6)
plot(nv.T, (nv.out('J_BK_k'))'./(nv.out('R_k')))
title('K^+ flux through the BK channel in AC (J\_BK\_k/R\_k)');   xlabel('Time [s]'); ylabel('K^+ flux [\muM/s]'); grid on

figure(3) % Plot variables from the Perivascular Space and Synaptic Cleft
subplot(1,2,1)
plot(nv.T, 0.001*nv.out('K_p'))
title('Potassium in the Perivascular Space (K_p)');   xlabel('Time [s]'); ylabel('[K^+]  (mM)'); grid on
subplot(1,2,2)
plot(nv.T, 0.001*nv.out('K_s'))
title('[K^+] in synaptic cleft: K_s'); xlabel('Time [s]'); ylabel('[K^+]_s [mM]'); grid on

figure(4)
suptitle('Smooth Muscle Cell')
subplot(2,2,1)
plot(nv.T, nv.out('Ca_i'))
title('[Ca^{2+}] in smooth muscle cell: Ca_i'); xlabel('Time [s]'); ylabel('[Ca^{2+}]_i [\muM]'); grid on
subplot(2,2,2)
plot(nv.T, nv.out('v_i'))
xlabel('Time [s]'); ylabel('v_i [mV]'); title('Membrane voltage smooth muscle cell: v_i'); grid on
subplot(2,2,3)
plot(nv.T, nv.out('J_VOCC_i'))
title('Ca^{2+} flux through the VOCC channel: J\_VOCC\_i'); xlabel('Time [s]'); ylabel('Ca^{2+} flux [\muM/s]'); grid on
subplot(2,2,4)
plot(nv.T,  nv.out('J_KIR_i'))
title('K^+ flux through the KIR channel: J\_KIR\_i'); xlabel('Time [s]'); ylabel('K^+ flux [\muM m/s]'); grid on

figure(5)
suptitle('Endothlial Cell')
subplot(2,2,1)
plot(nv.T, nv.out('Ca_j'))
title('[Ca^{2+}] in Endothlial cell: Ca_j'); xlabel('Time [s]'); ylabel('[Ca^{2+}]_i [\muM]'); grid on
subplot(2,2,2)
plot(nv.T, nv.out('v_j'))
xlabel('Time [s]'); ylabel('v_i [mV]'); title('Membrane voltage Endothlial cell: v_j'); grid on
subplot(2,2,[3:4])
plot(nv.T, nv.out('J_IP3_j'))
title('IP3 release in EC: J\_IP_3\_j'); xlabel('Time [s]'); ylabel(' [\muM/s]'); grid on

