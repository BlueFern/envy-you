%% Main script to plot FIG7 for NO paper
% open probability of BK channel in SMC 
% need to set dt(ind.v_i) & cGMP=0 in DEsyst.m
% need to put vivi in initCond.m 
% posibility to switch variable Ca_i on/off

clean
tic

% global variables
global CASE J_PLC startpulse lengthpulse C_Hillmann stretch_ch only_Koenig NVU Glu_start Glu_end wss_start wss_end c_w_switch t_wss_switch vivi nNOS_switch eNOS_switch
global cGMP_activated cGMP_switch Ca_switch %(1: Ca variable, 0: Ca locked at 1uM) 
%% NO pathway
global m %(cGMP coupling (0 - lowest influence to 2 - highest influence))
m = 2;
lalaa = 1;

c_w_switch = 1;
t_wss_switch = 1;

NO_switch = 1;
nNOS_switch = NO_switch;
eNOS_switch = NO_switch;

Ca_switch = 1;

for cGMP_activated = [0,1]
for vivi = -70:5:-30

%% Parameters to adjust the model:
t_start = 0;
t_end = 300;
startpulse  = 500;  % (s) 
lengthpulse = 2000;  % (s) 
Glu_start   = startpulse;
Glu_end     = startpulse + lengthpulse;
wss_start   = 100000; 
wss_end     = 120000;
CASE        = 2;    % (see all_constants.m for details)
J_PLC 		= 0.18;  % 0.18(steady) %0.4(fluctuating) (muM s-1) EC agonist concentration  
C_Hillmann  = 1;    % scaling factor for the Hai&Murphy rate constants (see all_constants.m for details)
stretch_ch  = 'ON'; % choose 'ON'/'OFF' to activate/deactivate stretch-activated channels in EC and SMC
only_Koenig = 'OFF';% choose 'ON'/'OFF' to simulate only the Koenigsberger model (other sub-models will still be considered, but the KIR channel is set to 0)
NVU         = 1;     % 1=NVU 1.0 , 2=NVU 1.1, 3=NVU 1.0 + EET, 4= NVU 1.0 + Ca2+
%% load the constants for the fluxes and pointers:
all_indices();
all_constants();
%% load the initial conditions of the system:
state0 = InitCond();
%% Ensure single filenames for the writing of data in other files
global csvfilename
csvfilename = 'Data_simulation.csv';
try
delete(csvfilename) % remove file, if present from older simulation.
end
%% Solve the proces from initial position tot Steady State:
options = odeset('OutputFcn',@odeprogWD,'Events',@odeabort,'Stats','on');%,'RelTol', 1e-03, 'AbsTol', 1e-06, 'MaxStep', 1); 
[t,state] = ode15s(@DEsyst,[t_start t_end],state0,options);

%% Write output and info to file/cmd
output.info.completiontime = toc;
fprintf('ODE solution time: %.3f seconds\n', output.info.completiontime)

%% Plot statement:
% % plot_all()
% 
DATA = csvread(csvfilename);

n = zeros(1,9);
a = zeros(1,31);
s = zeros(1,37);
e = zeros(1,19);
f= zeros(1,48); 
dfdt= zeros(1,48);
t= zeros(3,1);
input= zeros(1,3);

neoff   = 0;
acoff   = neoff + length(n);
smcoff  = acoff + length(a);
ecoff   = smcoff + length(s);
stoff   = ecoff  + length(e);
dfdtoff = stoff   + length(f);
tijdoff = dfdtoff+ length(dfdt);
inputoff= tijdoff+ 1;


time = DATA(:,length(DATA(1,:))-5);


vi(lalaa)   = vivi;
Po(lalaa)   = DATA(end,smcoff+flu.Kactivation_i);
vca3(lalaa) = DATA(end,smcoff+flu.v_Ca3);
cai(lalaa)  = state(end,ind.Ca_i) ;
wi(lalaa)  = state(end,ind.w_i) ;

lalaa = lalaa + 1

figure(112), plot(time, state(:,ind.cGMP))
hold all
figure(113), plot(time, state(:,ind.Ca_i))
hold all
figure(114), plot(time, state(:,ind.w_i))
hold all
figure(115), plot(time, state(:,ind.v_i))
hold all
end

figure(1), plot(vi,wi); xlabel('vi'), ylabel('wi'), title([cGMP_switch, Ca_switch]), legend('Basal state', 'Activated state'), hold all
% figure, plot(vi,Po); xlabel('vi'), ylabel('Kact'), title([cGMP_switch, Ca_switch])
figure(2), plot(cai,wi); xlabel('cai'), ylabel('wi'), title([cGMP_switch, Ca_switch]), legend('Basal state', 'Activated state'), hold all
% figure, plot(cai,Po); xlabel('cai'), ylabel('Kact'), title([cGMP_switch, Ca_switch])

clear vi Po vca3 cai wi
lalaa = 1;
end










% % figure; plot(time,DATA(:,smcoff+flu.M)+state(:,ind.AMp)+state(:,ind.AM)+state(:,ind.Mp)); hold on;
% % plot(time,DATA(:,smcoff+flu.M),'r'); plot(time,state(:,ind.Mp),'g'); plot(time,state(:,ind.AMp),'b');plot(time,state(:,ind.AM),'k');
% % legend('Total Myosin','[M]','[Mp]','[AMp]','[AM]')
% % title('New')
% 
% % to plot a single flux, type in plot(time,DATA(:,flu.(name))     
% % to plot a single state variable, type in plot(time,DATA(:,ind.(name))
% %(don't forget to put the offset!! e.g. smcoff+flu.1_c)
