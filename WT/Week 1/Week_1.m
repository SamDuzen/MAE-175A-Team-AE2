close all; clear all; clc;

%% raw data input
speed = [10,20,30,40,50,60,70,80]; %percent wind tunnel speed
pressure = [0.082,0.381,0.929,1.685,2.659,3.833,5.289,6.81]; %pressure reading (inches H2O)
SD_airflow = [8.92,8.32,8.66,7.8,7.64,7.93,6.35,5.34];
SD_pressure = [0.036,0.037,0.037,0.037,0.039,0.042,0.048,0.047];

P_amb = 101000; %Pa
T_amb = 21+273.15; %Kelvin
R_air = 287; %J/kg*K

%% calculate density
density = P_amb/(R_air*T_amb);

%% airspeed calculation

%pressure conversion to N/m
p_si = 248.84.*pressure;

%calculation of airspeed
U = sqrt((2*p_si)/density);

%linear fit - percent as a function of airspeed
    %useful for finding what percent to set wind tunnel for specific
    %airspeed as done in weeks 2 and 3
p2 = 1.251; %generated from curve fitting toolbox
p1 = 1.471; %generated from curve fitting toolbox
percent = @(x) p2 + p1*(x);

%linear fit - speed as a function of percent
    %useful for plotting
p4 = -0.8471; %generated from curve fitting toolbox
p3 = 0.6796; %generated from curve fitting toolbox
airspeed = @(x) p4 + p3*(x);

%vector creation
m = linspace(10,80,100);
n = airspeed(m);

speed_estimate = airspeed(speed);

%% Error Analysis

%Airspeed Error
dP = 248.84.*SD_pressure;
dU = (((density)./(2*p_si)).*(dP.^2)).^(1/2);

%% figure

figure()
hold on
plot(m,n,'b','LineWidth',1)
plot(speed,U,'ok','MarkerFaceColor','k','MarkerSize',4)
errorbar(speed,speed_estimate,dU,'m')
xlim([5,85])
ylim([0,60])
grid on
xlabel('Wind Tunnel Speed (%)')
ylabel('Airspeed (m/s)')
legend('Linear Fit - [Airspeed] = 0.6796*[Wind Tunnel Speed] - 0.8471','Data Points','Error Bars','Location','best','FontSize',8)





