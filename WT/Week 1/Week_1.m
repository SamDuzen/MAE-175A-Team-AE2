close all; clear all; clc;

%raw data input
speed = [10,20,30,40,50,60,70,80]; %percent wind tunnel speed
pressure = [0.082,0.381,0.929,1.685,2.659,3.833,5.289,6.81]; %pressure reading (inches H2O)

%pressure conversion to N/m
p_si = 248.84.*pressure;

%calculation of airspeed
U = sqrt((2*p_si)/(1.196384077));

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


%error
ls = 0.6726;
li = -1.2;
us = 0.6866;
ui = -0.4946;

lower = ls*speed + li;
upper = us*speed + ui;

error = upper-lower;

speed_estimate = airspeed(speed);



figure()
hold on
plot(m,n,'b','LineWidth',1)
plot(speed,U,'ok','MarkerFaceColor','k','MarkerSize',4)
errorbar(speed,speed_estimate,error,'m')
xlim([5,85])
ylim([0,60])
grid on
xlabel('Wind Tunnel Speed (%)')
ylabel('Airspeed (m/s)')
legend('Linear Fit - [Airspeed] = 0.6796*[Wind Tunnel Speed] - 0.8471','Data Points','Location','best','FontSize',8)