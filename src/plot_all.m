DATA = csvread(csvfilename);
close all
all_indices();
all_constants();

linestylle = ':';

n = zeros(1,9);
a = zeros(1,31);
s = zeros(1,34);
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

% plots report

min(DATA(:,acoff+flu.J_BK_k)./(DATA(:,stoff+ind.R_k)))
   
total=4;


figure(1)
hold all

subplot(total,1,1)
    plot( time,DATA(:,inputoff+1)*k_C);
title('Input signal from the neuron into the synaptic cleft')
xlabel('Time [s]')
ylabel('K^+ flux [\muM m/s]')

% subplot(total,1,1)
%     plot( time,DATA(:,acoff+flu.K_s));
% title('K^+ in SC')
% xlabel('Time [s]')
% ylabel('K_s [\muM]')

subplot(total,1,2)
plot(time,DATA(:,acoff+flu.v_k)*1000)
title('Membrane Potential of the astrocyte')
xlabel('Time [s]')
ylabel('v_k [mV]')

subplot(total,1,3)
[AX2,H12,H23] = plotyy(time,DATA(:,acoff+flu.J_BK_k)./(DATA(:,stoff+ind.R_k)*VR_pa)...
    ,time,DATA(:,smcoff+flu.J_KIR_i)/(VR_ps)...
    );
xlabel('Time [s]')
set(get(AX2(1),'Ylabel'),'String','J_{BK} [\muM/s]') 
set(get(AX2(2),'Ylabel'),'String','J_{KIR} [\muM/s]')
%legend('BK-channel','KIR-channel')
title('The contribution of the BK- and KIR-channel to K_p')

subplot(total,1,4)
plot(time,DATA(:,stoff+ind.K_p)...
     )
title('Potassium concentration in perivascular space')
xlabel('Time [s]')
ylabel('K_p  [\muM]')

%%

SMCtitle{1}='v_ coup_i';
SMCtitle{2}='Ca_ coup_i';
SMCtitle{3}='IP3_ coup_i';
SMCtitle{4}='rho_i';
SMCtitle{5}='J_ IP3_i';
SMCtitle{6}='J_ SRuptake_i';
SMCtitle{7}='J_ CICR_i';
SMCtitle{8}='J_ extrusion_i';
SMCtitle{9}='J_ leak_i';
SMCtitle{10}='J_ VOCC_i';
SMCtitle{11}='J _ NaCa_i';
SMCtitle{12}='J_ NaK_i';
SMCtitle{13}='J_ Cl_i';
SMCtitle{14}='J_ K_i';
SMCtitle{15}='Kactivation_i';
SMCtitle{16}='J_ degrad_i';

ECtitle{1}='v_ coup_j';
ECtitle{2}='Ca_ coup_j';
ECtitle{3}='IP3_ coup_j';
ECtitle{4}='rho_j';
ECtitle{5}='J_ 0_j';
ECtitle{6}='J_ IP3_j';
ECtitle{7}='J_ ERuptake_j';
ECtitle{8}='J_ CICR_j';
ECtitle{9}='J_ extrusion_j';
ECtitle{10}='J_ leak_j';
ECtitle{11}='J_ cation_j';
ECtitle{12}='J_ BKCa_j';
ECtitle{13}='J_ SKCa_j';
ECtitle{14}='J_ K_j';
ECtitle{15}='J_ R_j';
ECtitle{16}='J_ degrad_j';

figure(2)
set(gcf,'Name','SMC fluxes')
set(gcf,'Position', [24 62 1616 904],...
        'PaperPosition', [0.634517 6.34517 20.3046 15.2284]...
        );
for j = [1:3 5:16]
subplot(4,4,j)
plot(time,DATA(:,smcoff+j))
% h101(j+2) = gca();
xlabel('time in s')
ylabel(SMCtitle{j})
title(SMCtitle{j})
hold all
end

subplot(4,4,4)
plot(time,DATA(:,smcoff+flu.J_stretch_i))
xlabel('time in s')
ylabel('J stretch_i')
title('J stretch_i')
hold all

figure(3)
set(gcf,'Name','EC fluxes')
set(gcf,'Position', [24 62 1616 904],...
        'PaperPosition', [0.634517 6.34517 20.3046 15.2284]...
        );
for k = [1:3 5:16]
subplot(4,4,k)
plot(time,DATA(:,ecoff+k))
% h101(j+2) = gca();
xlabel('time in s')
ylabel(ECtitle{k}) 
title(ECtitle{k})
hold all
end

subplot(4,4,4)
plot(time,DATA(:,ecoff+flu.J_stretch_j))
ylabel('J stretch_j')
title('J stretch_j')
hold all

%% DFDT plot

figure(4)
set(gcf,'Name','DFDT')
set(gcf,'Position', [24 62 1616 904],...
        'PaperPosition', [0.634517 6.34517 20.3046 15.2284]...
        );
subplot(3,3,1)
plot(time, state(:,ind.Ca_i) )
h6(1) = gca();
xlabel('time in s')
ylabel('Ca_i in uM')
title('SMC [Ca^{2+}]')
hold all

subplot(3,3,2)
plot(time, state(:,ind.Ca_j) )
h6(2) = gca();
xlabel('time in s')
ylabel('Ca_j in uM')
title('EC [Ca^{2+}]')
hold all

subplot(3,3,3)
plot(time, state(:,ind.I_i) )
h6(8) = gca();
xlabel('time in s')
ylabel('I_i in uM')
title('SMC [IP_{3}]')
hold all

subplot(3,3,4)
plot(time, state(:,ind.I_j) )
h6(9) = gca();
xlabel('time in s')
ylabel('I_j in uM')
title('EC [IP_{3}]')
hold all

subplot(3,3,5)
plot(time, state(:,ind.s_i) )
h6(3) = gca();
xlabel('time in s')
ylabel('s_i in uM')
title('SR [Ca^{2+}]')
hold all

subplot(3,3,6)
plot(time, state(:,ind.s_j) )
h6(4) = gca();
xlabel('time in s')
ylabel('s_j in uM')
title('ER [Ca^{2+}]')
hold all

subplot(3,3,7)
plot(time, state(:,ind.v_i) )
h6(5) = gca();
xlabel('time in s')
ylabel('v_i in mV')
title('SMC membrane voltage')
hold all

subplot(3,3,8)
plot(time, state(:,ind.v_j) )
h6(6) = gca();
xlabel('time in s')
ylabel('v_j in mV')
title('EC membrane voltage')
hold all

subplot(3,3,9)
plot(time, state(:,ind.w_i) )
h6(7) = gca();
xlabel('time in s')
ylabel('w_i [-]')
title('open probability K^+ channel')
hold all


linkaxes(h6, 'x');

%% plot myosin crossbridge model

MCtitle{1}='[M]';
MCtitle{2}='[Mp]';
MCtitle{3}='[AMp]';
MCtitle{4}='[AM]';

figure(5)
set(gcf,'Name','Myosin model and radius')

subplot(3,2,1)
plot(time,DATA(:,smcoff+flu.M))
xlabel('Time')
ylabel('Fraction [-]')
title(MCtitle{1})
hold all

for j = 2:4
    subplot(3,2,j)
    plot(time,state(:,ind.Mp+j-2))
    xlabel('Time')
    ylabel('Fraction [-]')
    title(MCtitle{j})
    hold all
end

subplot(3,2,5)
plot(time,state(:,ind.AMp)+state(:,ind.AM))
xlabel('Time')
ylabel('Fraction [-]')
title(' F_r')
hold all

subplot(3,2,6)
plot(time,1e6*state(:,ind.R))
xlabel('Time')
ylabel('Radius in um')
title('Radius')
hold all


% close all;
figure(6)
% subplot(5,1,1)
% plot( time,DATA(:,inputoff+1)*k_C);
% title('Input signal from the neuron into the synaptic cleft')
% xlabel('Time [s]')
% ylabel('K^+ flux [\muM m/s]')

subplot(4,2,1)
plot( time,0.001*DATA(:,acoff+flu.K_s));
title('[K^+] in synaptic cleft')
xlabel('Time [s]')
ylabel('[K^+]_s [mM]')
hold all

subplot(4,2,2)
plot(time,(DATA(:,acoff+flu.J_BK_k)./(DATA(:,stoff+ind.R_k))))
title('K^+ flux through the BK channel')
xlabel('Time [s]')
ylabel('K^+ flux [\muM/s]')
hold all

subplot(4,2,3)
plot(time,0.001*DATA(:,stoff+ind.K_p))
title('[K^+] in perivascular space')
xlabel('Time [s]')
ylabel('[K^+]_s [mM]')
hold all

subplot(4,2,4)
plot(time,DATA(:,smcoff+flu.J_KIR_i))
title('K^+ flux through the KIR channel')
xlabel('Time [s]')
ylabel('K^+ flux [\muM m/s]')
hold all

subplot(4,2,5)
plot(time, state(:,ind.v_i) )
xlabel('Time [s]')
ylabel('v_i [mV]')
title('Membrane voltage smooth muscle cell')
hold all

subplot(4,2,6)
plot(time,DATA(:,smcoff+flu.J_VOCC_i))
title('Ca^{2+} flux through the VOCC channel')
xlabel('Time [s]')
ylabel('Ca^{2+} flux [\muM/s]')
%axis([0 250 -0.15 0])
hold all

subplot(4,2,7)
plot(time, state(:,ind.Ca_i) )
title('[Ca^{2+}] in smooth muscle cell')
xlabel('Time [s]')
ylabel('[Ca^{2+}]_i [\muM]')
hold all


% subplot(5,1,4)
% plot(time,state(:,ind.AMp)+state(:,ind.AM))
% xlabel('Time')
% ylabel('Fraction [-]')
% title(' F_r')

subplot(4,2,8)
plot(time,1e6*state(:,ind.R))
title('Radius')
xlabel('Time [s]')
ylabel('Radius [\mum]')
%axis([0 250 15 27])
hold all
% axis([0 250 15 27])




figure(7)
subplot(2,1,1)
plot(time, state(:,ind.Ca_i) )
xlabel('time in s')
ylabel('Ca_i in uM')
title('SMC [Ca^{2+}]')
hold all

subplot(2,1,2)
plot(time,1e6*state(:,ind.R))
title('Radius')
xlabel('Time [s]')
ylabel('Radius [\mum]')
hold all


%% Potassium plots

% figure;
% subplot(2,2,1)
% plot(time,state(:,ind.K_i))
% title('K^+ in SMC')
% subplot(2,2,2)
% plot(time,DATA(:,smcoff+flu.J_KIR_i))
% title('KIR channel outflux')
% subplot(2,2,3)
% plot(time,-DATA(:,smcoff+flu.J_K_i))
% title('J K_i outflux')
% subplot(2,2,4)
% plot(time,DATA(:,smcoff+flu.J_NaK_i))
% title('J NaK_i influx')




%% presentation
figure(99)
set(gcf, 'Position', [1 1 573 475]);
axes1 = axes('Position',[0.150779896013865 0.129554655870445 0.79896013864818 0.821862348178138],...
    'FontSize',16,...
    'FontName','Calibri',...
    'Color',[0.933333333333333 0.933333333333333 0.933333333333333]);
ylim(axes1,[15 30]);
box(axes1,'on');
hold(axes1,'all');
xlabel('Time (s)','FontSize',24,'FontName','Calibri');
ylabel('Vessel Radius (\mum)','FontSize',24,'FontName','Calibri');
plot(time,1e6*state(:,ind.R),'LineWidth',2);
% legend('Neuron','Endothelial Cell','Smooth Muscle Cell')

%% NO pathway

figure(8)
set(gcf,'Name','NO pathway')
set(gcf,'Position', [24 62 1616 904],...
        'PaperPosition', [0.634517 6.34517 20.3046 15.2284]...
        );
   
subplot(4,2,1)
    plot( time,DATA(:,inputoff+5)./1e3);
xlabel('Time [s]')
ylabel('[Glu] in mM s^{-1}')
% [AX3,H56,H67] = plotyy( time,DATA(:,inputoff+4),time,DATA(:,inputoff+5)./1e3);
% set(get(AX3(1),'Ylabel'),'String','wss in Pa') 
% set(get(AX3(2),'Ylabel'),'String','[Glu] in mM s^{-1}') 

subplot(4,2,2)
    plot( time,DATA(:,inputoff+4));
xlabel('Time [s]')
ylabel('wss in Pa')

subplot(4,2,3)
plot(time, state(:,ind.Ca_n) )
xlabel('time in s')
ylabel('[Ca^{2+}]_n in \muM')

subplot(4,2,6)
plot(time, state(:,ind.eNOS_act))
xlabel('time in s')
ylabel('[eNOS_{act}] in \muM')

subplot(4,2,5)
plot(time, state(:,ind.nNOS_act))
xlabel('time in s')
ylabel('[nNOS_{act}] in \muM')

subplot(4,2,7)
plot(time, state(:,ind.NOn),time, state(:,ind.NOa),time, state(:,ind.NOj),time, state(:,ind.NOi))
xlabel('time in s')
ylabel('[NO] in \muM')
legend('[NO]_n','[NO]_a','[NO]_j','[NO]_i','Location','NorthWest')

subplot(4,2,8)
plot(time, state(:,ind.cGMP))
xlabel('time in s')
ylabel('[cGMP] in \muM')


figure(9)
set(gcf,'Name','Report1 - NO pathway')
set(gcf,'Position', [24 62 734 904],...
        'PaperPosition', [0.634517 6.34517 20.3046 15.2284]...
        );
   
subplot(4,2,1)
% plot( time,DATA(:,inputoff+5)./1e3);
% ylabel('[Glu] in mM s^{-1}')
[AX4,H78,H89] = plotyy( time,DATA(:,inputoff+4),time,DATA(:,inputoff+3)./1e3);
set(get(AX4(1),'Ylabel'),'String','wss in Pa') 
set(get(AX4(2),'Ylabel'),'String','[Glu] in mM s^{-1}') 
% set(get(AX4(1),'YLim'),[0 4]) 


subplot(4,2,3)
plot(time, state(:,ind.Ca_n) )
%xlabel('time in s')
ylabel('[Ca^{2+}]_n in \muM')

subplot(4,2,5)
plot(time, state(:,ind.nNOS_act))
%xlabel('time in s')
ylabel('[nNOS_{act}]_n in \muM')

subplot(4,2,7)
plot(time, state(:,ind.NOn),time, state(:,ind.NOa),time, state(:,ind.NOj),time, state(:,ind.NOi))
xlabel('time in s')
ylabel('[NO] in \muM')
legend('[NO]_n','[NO]_a','[NO]_j','[NO]_i','Location','NorthWest')

subplot(4,2,2)
plot(time, state(:,ind.cGMP))
%xlabel('time in s')
ylabel('[cGMP] in \muM')

subplot(4,2,4)
plot(time,DATA(:,smcoff+flu.Kactivation_i))
ylabel('K_{activation_i} (dim.less)')

subplot(4,2,6)
plot(time,DATA(:,smcoff+flu.K2_c))
%xlabel('time in s')
ylabel('K_2 and K_5 (dim.less)')

subplot(4,2,8)
plot(time,1e6*state(:,ind.R))
xlabel('time in s')
ylabel('Radius in \mum')



figure(10)
set(gcf,'Name','Report2 - NO pathway')
set(gcf,'Position', [24 62 734 904],...
        'PaperPosition', [0.634517 6.34517 20.3046 15.2284]...
        );
   
subplot(4,2,1)
% plot( time,DATA(:,inputoff+4));
% ylabel('wss in Pa')
[AX3,H56,H67] = plotyy( time,DATA(:,inputoff+4),time,DATA(:,inputoff+5)./1e3);
set(get(AX3(1),'Ylabel'),'String','wss in Pa') 
set(get(AX3(2),'Ylabel'),'String','[Glu] in mM s^{-1}') 

subplot(4,2,3)
plot(time, state(:,ind.eNOS_act))
%xlabel('time in s')
ylabel('[eNOS_{act}]_j in \muM')

subplot(4,2,5)
plot(time, state(:,ind.NOn),time, state(:,ind.NOa),time, state(:,ind.NOj),time, state(:,ind.NOi))
%xlabel('time in s')
ylabel('[NO] in \muM')
legend('[NO]_n','[NO]_a','[NO]_j','[NO]_i','Location','NorthWest')

subplot(4,2,7)
plot(time, state(:,ind.cGMP))
xlabel('time in s')
ylabel('[cGMP] in \muM')

subplot(4,2,2)
plot(time,DATA(:,smcoff+flu.Kactivation_i))
%xlabel('time in s')
ylabel('K_{activation_i} (dim.less)')

subplot(4,2,4)
plot(time,DATA(:,smcoff+flu.K2_c))
%xlabel('time in s')
ylabel('K_2 and K_5 (dim.less)')

subplot(4,2,6)
plot(time,1e6*state(:,ind.R))
xlabel('time in s')
ylabel('Radius in \mum')


figure(11)
set(gcf,'Name','Report3 - NO pathway')
set(gcf,'Position', [24 62 734 904],...
        'PaperPosition', [0.634517 6.34517 20.3046 15.2284]...
        );
   
subplot(4,2,1)
[AX1,H34,H45] = plotyy( time,DATA(:,inputoff+4),time,DATA(:,inputoff+5)./1e3);
set(get(AX1(1),'Ylabel'),'String','wss in Pa') 
set(get(AX1(2),'Ylabel'),'String','[Glu] in mM s^{-1}') 

subplot(4,2,3)
plot(time, state(:,ind.eNOS_act))
%xlabel('time in s')
ylabel('[eNOS_{act}]_j in \muM')

subplot(4,2,5)
plot(time, state(:,ind.nNOS_act))
%xlabel('time in s')
ylabel('[nNOS_{act}]_n in \muM')

subplot(4,2,7)
plot(time, state(:,ind.NOn),time, state(:,ind.NOj),time, state(:,ind.NOi))
xlabel('time in s')
ylabel('[NO] in \muM')
legend('[NO]_n','[NO]_j','[NO]_i','Location','NorthWest')

subplot(4,2,2)
plot(time, state(:,ind.cGMP))
% xlabel('time in s')
ylabel('[cGMP] in \muM')

subplot(4,2,4)
plot(time,DATA(:,smcoff+flu.Kactivation_i))
%xlabel('time in s')
ylabel('K_{activation_i} (dim.less)')

subplot(4,2,6)
plot(time,DATA(:,smcoff+flu.K2_c))
%xlabel('time in s')
ylabel('K_2 and K_5 (dim.less)')

subplot(4,2,8)
plot(time,1e6*state(:,ind.R))
xlabel('time in s')
ylabel('Radius in \mum')






figure(12)
set(gcf,'Name','Report4 - NO pathway')
set(gcf,'Position', [24 62 734 904],...
        'PaperPosition', [0.634517 6.34517 20.3046 15.2284]...
        );
   
subplot(4,2,1)
plot(time,DATA(:,smcoff+flu.Kactivation_i))
ylabel('K_{activation_i} (dim.less)')


subplot(4,2,3)
plot(time, state(:,ind.w_i))
ylabel('w_i')

subplot(4,2,5)
plot(time,DATA(:,smcoff+flu.J_K_i))
ylabel('J_{K_i}')

subplot(4,2,7)
plot(time, state(:,ind.v_i))
ylabel('v_i')

subplot(4,2,2)
plot(time,DATA(:,smcoff+flu.J_VOCC_i))
ylabel('J_{VOCC_i}')

subplot(4,2,4)
plot(time,state(:,ind.Ca_i))
ylabel('Ca_i in \mu M')

subplot(4,2,6)
plot(time,DATA(:,smcoff+flu.K1_c))
ylabel('K_1 and K_6 (dim.less)')

subplot(4,2,8)
plot(time,1e6*state(:,ind.R))
xlabel('time in s')
ylabel('Radius in \mum')

