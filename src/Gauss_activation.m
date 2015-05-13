%% Gauss activation
close all
% coordinates
% xlin = linspace(-0.0128,0.0128,33);         %n_bif = 13
% ylin = linspace(-0.0128,0.0128,33);         
xlin = linspace(-0.0016,0.0016,33);           %n_bif = 7
ylin = linspace(-0.0016,0.0016,33);    

A = 3;
ramp = 0.003; %0.002;
x_centre = 0;%0.0008;%0.0064;
y_centre = 0;%0.0008;%0.0064;
PLC_min = 0.18;
PLC_max = 0.4;

[x,y] = meshgrid(xlin,ylin);
% out = min(1,A .* exp(- (((x-x_centre).^2+(y-y_centre).^2) ./ (2 * ramp.^2))));
% out = -20.3125 .* x + 0.44; % gradient for nbif = 13 (Jplc)
%out = -162.5 .* x + 0.44; % gradient for nbif = 7 (Jplc) 0.18-0.4
% out = 937.5 .* x + 1.5; % K_input 0-3
%out = PLC_min + (PLC_max-PLC_min) * (312.5 .* x + 0.5);
%%out = -100 .* x -100 .* y + 0.44; %diagonal
%out = (312.5 .* y + 0.5); % gradient 0-1
out = 625 * y + 1.2; % 0-1


figure
surf(x,y,out)
xlabel('x')
ylabel('y')
zlabel('out') 


%%

K_input_min = 0;
K_input_max = 2.5;
ramp = 0.003; % 0.004; 
ampl = 3;
x_centre = 0;%-0.0008; % 0;
y_centre = 0;%-0.0008; % 0;
t_up   = 200;
t_down = 800;
lengthpulse = t_down - t_up;	
lengtht1 = 10;
F_input = 2.5;
t0 = t_up;
t1 = t0 + lengtht1;
t2 = t0 + lengthpulse;
t3 = t1 + lengthpulse;
alpha = 2;
beta = 5;
deltat= 10; 	
gab = factorial(alpha + beta - 1);
ga = factorial(alpha - 1);
gb = factorial(beta - 1);
K_space = fmin(1.0,ampl*((K_input_max-K_input_min) * exp(- ((pow((x-x_centre),2)+pow((y-y_centre),2)) / (2 * pow(ramp,2)))))); 



