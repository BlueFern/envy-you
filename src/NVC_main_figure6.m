%% Main script to plot FIG6 for NO paper
% Distribution of NO concentration among cell types

clean
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

for c_w_switch = [1];
for t_wss_switch = [1];
nNOS_switch = 1;
eNOS_switch = 1;

c_w_switch
t_wss_switch

% for vivi = -100:10:200

%% Parameters to adjust the model:
t_start = 0;
t_end = 1000;
startpulse  = 200;  % (s) 
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

%% FIG6
figure(7)
set(gcf, 'Position', [400 300 700 300]);
plot(time, state(:,ind.NO_n),time, state(:,ind.NO_k),time, state(:,ind.NO_j),time, state(:,ind.NO_i),'LineWidth',1);
xlabel('Time (s)');
ylabel('[NO] (\muM)')
legend('NE','AC','EC','SMC','Location','NorthWest')
hold all
% ylim([20 32]);
hold all

% figure(13)
% set(gcf, 'Position', [400 300 700 300]);
% plot(time, state(:,ind.NOn_max),time, state(:,ind.NOk_max),time, state(:,ind.NOj_max),time, state(:,ind.NOi_max),'LineWidth',1);
% xlabel('Time (s)');
% ylabel('[NO] (\muM)')
% legend('NE\_max','AC\_max','EC\_max','SMC\_max','Location','NorthWest')
% hold all
% % ylim([20 32]);
% hold all

figure(8)
set(gcf, 'Position', [400 300 700 300]);
y = [state(end,ind.NO_n) state(end,ind.NO_k) state(end,ind.NO_i) state(end,ind.NO_j)];
bar(y,0.4)
ylabel('[NO] (\muM)')
hold all

% figure(12)
% set(gcf, 'Position', [400 300 700 300]);
% y = [state(end,ind.NOn_max) state(end,ind.NOk_max) state(end,ind.NOi_max) state(end,ind.NOj_max)];
% bar(y,0.4)
% ylabel('[NO] (\muM)')
% hold all


% figure(9)
% set(gcf, 'Position', [400 300 700 300]);
% y = [1.0615 0.2003; 0.6621 0.1438; 0.2627 0.08725; 0.2538 0.08598];
% bar(y) %,0.4)
% ylabel('[NO] (\muM)')
% legend('[NO]_n', '[NO]_j')
% hold all


% figure(10)
% set(gcf, 'Position', [400 300 700 300]);
% NOn = 1.0615;
% NOj = 0.2003; %0.08598;
% set(gcf, 'Position', [400 300 700 300]);
% y = [1.0615/NOn 0.2003/NOj; 0.6621/NOn 0.1438/NOj; 0.2627/NOn 0.08725/NOj; 0.2538/NOn 0.08598/NOj];
% bar(y) %,0.4)
% ylabel('[NO] (\muM)')
% legend('[NO]_n', '[NO]_j')

% figure(11)
% set(gcf, 'Position', [400 300 700 300]);
% y = [1.065 1.0615 0.2003; 0.6656 0.6621 0.1438; 0.2663 0.2627 0.08725; 0.2573 0.2538 0.08598];
% bar(y) %,0.4)
% Labels = {'NE', 'AC', 'SMC', 'EC'};
% set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
% ylabel('[NO] (\muM)')
% legend('maximum neuronal and endothelial NO release', 'only neuronal NO release', 'only endothelial NO release')



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
