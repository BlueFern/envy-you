%% Demonstration script for new NVU model reverting back to version 1.0
%      Author: Michelle Goodman
%      Date: 3/2/2015

%      For reference 
%           "Old" refers to version 1.0 orgional prior to Ca and speedfix 
%           "New" refers to version 1.1 with Ca in astrocyte and speedfix edits

% First copy from old code running types
close all; clc; clear
odeopts = odeset('RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 1); 
% odeopts = odeset('Vectorized', 1); This was the latest

nv = NVU(Astrocyte(), ...
    WallMechanics(), ...
    SMCEC('J_PLC', 0.4), ...
    'odeopts', odeopts);


%% Changing parameters, initial conditions, simulation time to be the same as the old
nv.smcec.params.J_PLC = 0.4; % new is 0.4;  % (muM s-1) EC agonist concentration  
%nv.astrocyte.u0(nv.astrocyte.index.K_p) = 3e3;
nv.T = linspace(0, 500, 1245);
nv.simulate();
% Plot, e.g. Ca_i
figure(1)
plot(nv.T, nv.out('Ca_i'))
xlabel('time (s)')
ylabel('[Ca^{2+}] (\muM)')

figure(2)
subplot(3, 2, 1)
plot(nv.T, nv.out('M'))
xlabel('Time')
ylabel('Fraction [-]')
title('[M]')

subplot(3, 2, 2)
plot(nv.T, nv.out('Mp'))
xlabel('Time')
ylabel('Fraction [-]')
title('[Mp]')

subplot(3, 2, 3)
plot(nv.T, nv.out('AMp'))
xlabel('Time')
ylabel('Fraction [-]')
title('[AMp]')

subplot(3, 2, 4)
plot(nv.T, nv.out('AM'))
xlabel('Time')
ylabel('Fraction [-]')
title('[AM]')

subplot(3, 2, 5)
plot(nv.T, nv.out('F_r'))
xlabel('Time')
ylabel('Fraction [-]')
title('F_r')
%% IMPORTANT QUANTITY
subplot(3, 2, 6)
plot(nv.T, 1e6 * nv.out('R'))
xlabel('Time')
ylabel('\mu m')
title('Radius')
cd('C:\Users\mlg77\Local Documents\NVU Documentation')
%desired_save  = [nv.T, 1e6 * nv.out('R'), ];
desired_save  = [nv.T,1e6 * nv.out('R'), nv.out('Ca_i') , nv.out('J_VOCC_i')', nv.out('v_i'), nv.out('J_KIR_i')', 0.001*nv.out('K_p'), (nv.out('J_BK_k'))'./(nv.out('R_k')), 0.001*nv.out('K_s')' ];
desired_save(1,10) = nv.smcec.params.J_PLC;
save('Radius My attempt at old, mu m', 'desired_save')
cd('C:\Users\mlg77\Local Documents\Git\envy-you-fork\src')

%%


%% plot all state variables:
for i = 1:15
    figure(3)
    set(gcf,'name','Astrocyte')
    subplot(3,5,i)
    plot(nv.T, nv.out(nv.astrocyte.varnames{i}))
    h1(i) = gca();
    ylabel(nv.astrocyte.varnames{i})
    hold on
end
    linkaxes(h1, 'x');


for i = 1:10
    figure(4)
    set(gcf,'name','ECSMC & WallMechanics')
    subplot(3,5,i)
    plot(nv.T, nv.out(nv.smcec.varnames{i}))
    h2(i) = gca();
    ylabel(nv.smcec.varnames{i})
    hold on
end

    
for i = 1:4
    figure(4)
    subplot(3,5,10+i)
    plot(nv.T, nv.out(nv.wall.varnames{i}))
    h2(10+i) = gca();
    ylabel(nv.wall.varnames{i})
    hold on
end
    linkaxes(h2, 'x');
    hold on