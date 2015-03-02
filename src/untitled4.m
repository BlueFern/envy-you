
figure(44)
subplot(2,6,1)
plot(time,DATA(:,smcoff+flu.J_K_i))
xlabel('time in s')
ylabel('J\_K_i')
hold all

subplot(2,6,2)
plot(time,DATA(:,smcoff+flu.J_VOCC_i))
xlabel('time in s')
ylabel('J\_VOCC_i')
hold all

subplot(2,6,3)
plot(time, state(:,ind.Ca_i))
xlabel('time in s')
ylabel('[Ca]_i (\muM)')
hold all

subplot(2,6,4)
plot(time, state(:,ind.v_i))
xlabel('time in s')
ylabel('v_i')
hold all

subplot(2,6,5)
plot(time, state(:,ind.w_i))
xlabel('time in s')
ylabel('w_i')
hold all

subplot(2,6,6)
plot(time,DATA(:,smcoff+flu.J_NaK_i))
ylabel('J\_NaK')
xlabel('time in s')
hold all

subplot(2,6,7)
plot(time,DATA(:,smcoff+flu.J_Cl_i))
xlabel('time in s')
ylabel('J\_Cl_i')
hold all

subplot(2,6,8)
plot(time,DATA(:,smcoff+flu.J_NaCa_i))
ylabel('J\_NaCa_i')
xlabel('time in s')
hold all

subplot(2,6,9)
plot(time,DATA(:,smcoff+flu.J_stretch_i))
ylabel('J\_stretch_i')
xlabel('time in s')
hold all

subplot(2,6,10)
plot(time,DATA(:,smcoff+flu.J_KIR_i))
ylabel('J\_KIR_i')
xlabel('time in s')
hold all

subplot(2,6,11)
plot(time,DATA(:,smcoff+flu.v_coup_i))
ylabel('V\_coupl')
xlabel('time in s')
hold all

% subplot(2,6,12)
% plot(time,DATA(:,smcoff+flu.J_VOCC_i))
% ylabel('tau\_wss')
% xlabel('time in s')
hold all