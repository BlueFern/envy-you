%% Main script to plot FIG2 for NO paper
% Influence of length of stimulus on vasodilation
% Option to switch on and off c_w, t_wss, NO_pathway! -> legend entries

clean
tic

% global variables
global CASE J_PLC startpulse lengthpulse C_Hillmann stretch_ch only_Koenig NVU Glu_start Glu_end wss_start wss_end c_w_switch t_wss_switch vivi NO_switch
global nNOS_switch eNOS_switch 

%% NO pathway
global m %(cGMP coupling (0 - lowest influence to 2 - highest influence))
m = 2;
lalaa = 1;
% 
% cai = [];
% logcai = [];
% Po = [];
% vca3 = [];
% rcgmp1 = [];
% vii = [];
% timee = [];
% for c_w_switch = [1]
% 
% for t_wss_switch = [1]
% 
% for NO_switch = [1]
% for switch1 = [1]
    
% c_w_switch = switch1;
% t_wss_switch = switch1;
% NO_switch = switch1;
%     

startpulse  = 300;  % (s) 
for lengthpulse = [20,500]  % (s) ------> give different legths here!


for c_w_switch = [1];
for t_wss_switch = [1];
for NO_switch = [1];
nNOS_switch = NO_switch;
eNOS_switch = NO_switch;
    
c_w_switch
t_wss_switch
NO_switch
% for vivi = -100:10:200

%% Parameters to adjust the model:
t_start = 0;
t_end = 1200;
Glu_start   = startpulse;
Glu_end     = startpulse+lengthpulse;
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
options = odeset('OutputFcn',@odeprogWD,'Events',@odeabort,'Stats','on'); %,'RelTol', 1e-03, 'AbsTol', 1e-06, 'MaxStep', 1); 
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
f= zeros(1,42); 
dfdt= zeros(1,42);
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

figure(4), plot(time, state(:,ind.R))
legend('1','2','3','4','5','6','7','8')
hold all

figure(5)
subplot(2,6,1)
plot(time, state(:,ind.nNOS_act))
xlabel('time in s')
ylabel('[nNOS_{act}]_n in \muM')
legend('1','2','3','4','5','6','7','8')
hold all

subplot(2,6,2)
plot(time, state(:,ind.eNOS_act))
xlabel('time in s')
ylabel('[eNOS_{act}]_n in \muM')
legend('1','2','3','4','5','6','7','8')
hold all

subplot(2,6,3)
plot(time, state(:,ind.cGMP))
xlabel('time in s')
ylabel('[cGMP] in \muM')
legend('1','2','3','4','5','6','7','8')
hold all

subplot(2,6,4)
plot(time,DATA(:,smcoff+flu.K2_c))
xlabel('time in s')
ylabel('K_2 and K_5 (dim.less)')
legend('1','2','3','4','5','6','7','8')
hold all

% subplot(2,6,5)
% plot(time, state(:,ind.w_i))
% xlabel('time in s')
% ylabel('w_i')
% legend('1','2','3','4','5','6','7','8')
% hold all

subplot(2,6,5)
plot(time, state(:,ind.w_i))
xlabel('time in s')
ylabel('w_i')
legend('1','2','3','4','5','6','7','8')
hold all

subplot(2,6,6)
plot(time, state(:,ind.R)*1e6)
ylabel('Radius in um')
xlabel('time in s')
legend('1','2','3','4','5','6','7','8')
hold all

subplot(2,6,7)
plot(time, state(:,ind.NO_n),time, state(:,ind.NO_k),time, state(:,ind.NO_j),time, state(:,ind.NO_i))
xlabel('time in s')
ylabel('[NO] in \muM')
legend('1[NO]_n','1[NO]_a','1[NO]_j','1[NO]_i','2[NO]_n','2[NO]_a','2[NO]_j','2[NO]_i','3[NO]_n','3[NO]_a','3[NO]_j','3[NO]_i','4[NO]_n','4[NO]_a','4[NO]_j','4[NO]_i','Location','NorthWest')
hold all

subplot(2,6,8)
plot(time, DATA(:,smcoff+flu.R_cGMP2))
ylabel('R\_cGMP2')
xlabel('time in s')
legend('1','2','3','4','5','6','7','8')
hold all

subplot(2,6,9)
plot(time, DATA(:,smcoff+flu.Act_eNOS_Ca))
ylabel('Act\_eNOS\_Ca')
xlabel('time in s')
legend('1','2','3','4','5','6','7','8')
hold all

subplot(2,6,10)
plot(time, DATA(:,smcoff+flu.Act_eNOS_wss))
ylabel('Act\_eNOS\_wss')
xlabel('time in s')
legend('1','2','3','4','5','6','7','8')
hold all

subplot(2,6,11)
plot(time, DATA(:,smcoff+flu.P_NOj_eNOS))
ylabel('P\_NOj\_eNOS')
xlabel('time in s')
legend('1','2','3','4','5','6','7','8')
hold all

subplot(2,6,12)
plot(time, DATA(:,neoff+flu.tau_w))
ylabel('tau\_wss')
xlabel('time in s')
legend('1','2','3','4','5','6','7','8')
hold all

%% FIG2
figure(7)
set(gcf, 'Position', [400 300 700 300]);
plot(time, state(:,ind.R)*1e6,'LineWidth',1,'LineStyle','--');
ylabel('Radius (\mum)')
xlabel('Time (s)');
ylim([20 32]);
xlim([0 1200]);
legend('100 s stimulation','500 s stimulation')
% legend('1','2','3','4','5','6','7','8')
hold all


end
end
end
end

% 
% %% save figures & parameters
% %save_all()
% 
% 
% % to create .tikz figures:
% % matlab2tikz('test.tikz', 'height', '\figureheight', 'width', '\figurewidth');
% 
% % 
% % % figure
% % % subplot(4,1,1)
% % % plot(time,DATA(:,smcoff+flu.Kactivation_i),'-x')
% % % subplot(4,1,2)
% % % plot(time,DATA(:,smcoff+flu.v_Ca3),'-x')
% % % subplot(4,1,3)
% % % plot(time,state(:,ind.Ca_i),'-x')
% % % subplot(4,1,4)
% % % plot(time,DATA(:,smcoff+flu.R_cGMP1),'-x')
% % 
% % foo0 = state(:,ind.Ca_i);
% % foo1 = DATA(:,smcoff+flu.Kactivation_i);
% % foo2 = DATA(:,smcoff+flu.v_Ca3);
% % foo3 = DATA(:,smcoff+flu.R_cGMP1);
% % foo4 = state(:,ind.v_i);
% % 
% % cai{lalaa} = mean(foo0(200:end));
% % % logcai{lalaa} = log10(cai);
% % Po{lalaa} = mean(foo1(200:end));
% % vca3{lalaa} = mean(foo2(200:end));
% % rcgmp1{lalaa} = mean(foo3(200:end));
% % vii{lalaa} = mean(foo4(200:end));
% % timee{lalaa} = time;
% % 
% % clear foo0 foo1 foo2 foo3 foo4
% % 
% % lalaa = lalaa + 1
% % end
% 
% 




% vi(lalaa)   = vivi;
% Po(lalaa)   = DATA(end,smcoff+flu.Kactivation_i);
% vca3(lalaa) = DATA(end,smcoff+flu.v_Ca3);
% cai(lalaa)  = state(end,ind.Ca_i) ;
% wi(lalaa)  = state(end,ind.w_i) ;
% 
% lalaa = lalaa+1
% end
% 
% figure, plot(vi,wi); xlabel('vi'), ylabel('wi')
% figure, plot(vi,Po); xlabel('vi'), ylabel('Kact')
% figure, plot(cai,wi); xlabel('cai'), ylabel('wi')
% figure, plot(cai,Po); xlabel('cai'), ylabel('Kact')
% 
% figure(112), plot(time, state(:,ind.cGMP))
% hold all
% figure(113), plot(time, state(:,ind.R))
% hold all
% end











% % figure; plot(time,DATA(:,smcoff+flu.M)+state(:,ind.AMp)+state(:,ind.AM)+state(:,ind.Mp)); hold on;
% % plot(time,DATA(:,smcoff+flu.M),'r'); plot(time,state(:,ind.Mp),'g'); plot(time,state(:,ind.AMp),'b');plot(time,state(:,ind.AM),'k');
% % legend('Total Myosin','[M]','[Mp]','[AMp]','[AM]')
% % title('New')
% 
% % to plot a single flux, type in plot(time,DATA(:,flu.(name))     
% % to plot a single state variable, type in plot(time,DATA(:,ind.(name))
% %(don't forget to put the offset!! e.g. smcoff+flu.1_c)
