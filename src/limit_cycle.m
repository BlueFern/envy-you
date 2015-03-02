figure, 
subplot(1,3,1)
plot(time(200:466),state((200:466),ind.v_i))
xlabel('time')
ylabel('v\_i')

subplot(1,3,2)
plot(time(200:466),state((200:466),ind.w_i))
xlabel('time')
ylabel('w\_i')

subplot(1,3,3)
plot(state((200:466),ind.v_i),state((200:466),ind.w_i))
xlabel('v\_i')
ylabel('w\_i')