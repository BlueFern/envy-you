time = 0:0.001:200;


flux_min = 0;
%flux_max = 1;  
t_up   = 20;
t_down = 80;
lengthpulse = t_down - t_up;
lengtht1 = 10;
t0 = t_up;	
t1 = t0 + lengtht1;
%ramp = 0.001;	
%x_centre = -0.0008;
%y_centre = -0.0008;
for i = 1:length(time)
    flux_time(i) = 0.5 * tanh((time(i)-t0)/0.0005) - 0.5 * tanh((time(i)-t1-lengthpulse)/0.0005);
    %flux_space = (flux_max-flux_min) * exp(- ((pow((x-x_centre),2)+pow((y-y_centre),2)) / (2 * pow(ramp,2)))); 
    %flux_out = flux_min + flux_time * flux_space;
    flux_out(i) = flux_min + flux_time(i);
end

figure, plot(time, flux_out)
%% 
K_input_min = 0;
K_input_max = 2.5;
ramp = 0.001; % 0.004; 
%x_centre = -0.0008; % 0;
%y_centre = -0.0008; % 0;
F_input = 2.5;
t2 = t0 + lengthpulse;
t3 = t1 + lengthpulse;
alpha = 2;
beta = 5;
deltat= 10; 	
gab = factorialK(alpha + beta - 1)
ga = factorialK(alpha - 1)
gb = factorialK(beta - 1)
%K_space = (K_input_max-K_input_min) * exp(- ((pow((x-x_centre),2)+pow((y-y_centre),2)) / (2 * pow(ramp,2)))); 
for i = 1:length(time)
    if (time(i) >= t0 && time(i) <= t1) 
        K_time = F_input * gab / (ga * gb) * ((1-(time(i)-t0) / deltat))^(beta - 1) * (((time(i) - t0) / deltat))^(alpha-1);  
    elseif (time(i) >= t2 && time(i) <= t3) 
        K_time = - F_input;
    else 
        K_time = 0;
    end
    K_out(i) = K_input_min + K_time;
end
figure, plot(time, K_out)

%%

for i = 1:length(time)
    out(i) = getRef(time(i),'ft');
end

figure, plot(time, out)


