% clean
tic

% global variables
global CASE J_PLC startpulse lengthpulse C_Hillmann stretch_ch only_Koenig NVU Glu_start Glu_end wss_start wss_end c_w_switch t_wss_switch vivi nNOS_switch eNOS_switch

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

global caicai ccGMP

% caicai = 1; %uM  
vivi = 40; %mV
ccGMP = 11; %uM
for caicai = logspace(-9, -3, 30)
    
c_w_switch = 1;
nNOS_switch = 1;
eNOS_switch = 1;
t_wss_switch = 1;


%% Parameters to adjust the model:
t_start = 0;
t_end = 200;
startpulse  = 400;  % (s) 
lengthpulse = 200;  % (s) 
Glu_start   = 400;
Glu_end     = 600;
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
f= zeros(1,46); 
dfdt= zeros(1,46);
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


%% FIG3
figure(7)
set(gcf, 'Position', [400 300 700 300]);
plot(time, state(:,ind.R)*1e6,'LineWidth',1);
ylabel('Radius (\mum)')
xlabel('Time (s)');
ylim([17 31]);
hold all


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



vi(lalaa)   = vivi;
Po(lalaa)   = DATA(end,smcoff+flu.Kactivation_i);
vca3(lalaa) = DATA(end,smcoff+flu.v_Ca3);
cai(lalaa)  = state(end,ind.Ca_i) ;
wi(lalaa)  = state(end,ind.w_i) ;

lalaa = lalaa+1

end

figure(1), plot(vi,wi); xlabel('vi'), ylabel('wi')
hold all
figure(2), plot(vi,Po); xlabel('vi'), ylabel('Kact')
hold all
figure(3), plot(cai,wi); xlabel('cai'), ylabel('wi')
hold all
figure(4), plot(cai,Po); xlabel('cai'), ylabel('Kact')
hold all










% % figure; plot(time,DATA(:,smcoff+flu.M)+state(:,ind.AMp)+state(:,ind.AM)+state(:,ind.Mp)); hold on;
% % plot(time,DATA(:,smcoff+flu.M),'r'); plot(time,state(:,ind.Mp),'g'); plot(time,state(:,ind.AMp),'b');plot(time,state(:,ind.AM),'k');
% % legend('Total Myosin','[M]','[Mp]','[AMp]','[AM]')
% % title('New')
% 
% % to plot a single flux, type in plot(time,DATA(:,flu.(name))     
% % to plot a single state variable, type in plot(time,DATA(:,ind.(name))
% %(don't forget to put the offset!! e.g. smcoff+flu.1_c)
