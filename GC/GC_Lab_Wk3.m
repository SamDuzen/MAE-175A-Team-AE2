clear all; close all; clc;
%% Read Data

%Read Open Loop Data
    %Trial 5 Open Loop
    fileID = fopen('Wk3_OL_5.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{1} = cell2mat(data); % Convert to numeric matrix
    Sample{1} = data{1}(:,1);  
    Time{1} = data{1}(:,2);     
    Encoder4Pos{1} = data{1}(:,3);  
    ControlEffort2{1} = data{1}(:,4);

%Plot Open Loop Data
figure()
hold on
plot(Time{1},Encoder4Pos{1})
xline(3,'--k')
grid on
xlabel('Time [s]')
ylabel('Encoder 4 Position [counts]')
title('3000ms Dwell Time')

%Velocity Calculation
i=1;
k=1;
n = 1;
while i < 130
    VelTime(k) = (Time{1}(i+n)+Time{1}(i))/2;
    OL_Vel(k) = (Encoder4Pos{1}(i+n)-Encoder4Pos{1}(i))/(Time{1}(i+n)-Time{1}(i));
    k = k+1;
    i = i+n;
end

figure()
hold on
plot(VelTime,OL_Vel)
xline(3,'--k')
xlabel('Time [s]')
ylabel('Encoder 4 Velocity [Counts/s]')
title('3000ms Dwell Time')
grid on


%% Sensitivity Analysis

%Extract data from maelab.m figures
fig1 = openfig('Wk3_Reference.fig');
    axObjs = fig1.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{1} = dataObjs(1).XData; %time
    y{1} = dataObjs(1).YData; %encoder position
fig2 = openfig('Wk3_Wn_09.fig');
    axObjs = fig2.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{2} = dataObjs(1).XData;
    y{2} = dataObjs(1).YData;
fig3 = openfig('Wk3_Wn_11.fig');
    axObjs = fig3.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{3} = dataObjs(1).XData;
    y{3} = dataObjs(1).YData;
fig4 = openfig('Wk3_Beta_09.fig');
    axObjs = fig4.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{4} = dataObjs(1).XData;
    y{4} = dataObjs(1).YData;
fig5 = openfig('Wk3_Beta_11.fig');
    axObjs = fig5.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{5} = dataObjs(1).XData;
    y{5} = dataObjs(1).YData;
fig6 = openfig('Wk3_K_09.fig');
    axObjs = fig6.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{6} = dataObjs(1).XData;
    y{6} = dataObjs(1).YData;
fig7 = openfig('Wk3_K_11.fig');
    axObjs = fig7.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{7} = dataObjs(1).XData;
    y{7} = dataObjs(1).YData;

    %TEST
fig8 = openfig('test2.fig');
    axObjs = fig8.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{8} = dataObjs(1).XData; %time
    y{8} = dataObjs(1).YData; %encoder position

%Velocity Calculation
i=1;
k=1;
n = 5;
while i < 1340
    VelTime(k) = (x{1}(i+n)+x{1}(i))/2;
    Vel{1}(k) = (y{1}(i+n)-y{1}(i))/(x{1}(i+n)-x{1}(i));
    Vel{2}(k) = (y{2}(i+n)-y{2}(i))/(x{2}(i+n)-x{2}(i));
    Vel{3}(k) = (y{3}(i+n)-y{3}(i))/(x{3}(i+n)-x{3}(i));
    Vel{4}(k) = (y{4}(i+n)-y{4}(i))/(x{4}(i+n)-x{4}(i));
    Vel{5}(k) = (y{5}(i+n)-y{5}(i))/(x{5}(i+n)-x{5}(i));
    Vel{6}(k) = (y{6}(i+n)-y{6}(i))/(x{6}(i+n)-x{6}(i));
    Vel{7}(k) = (y{7}(i+n)-y{7}(i))/(x{7}(i+n)-x{7}(i));
    k = k+1;
    i = i+n;
end

figure()
hold on
plot(VelTime,Vel{1})
plot(VelTime,Vel{2})
plot(VelTime,Vel{3})
plot(VelTime,Vel{4})
plot(VelTime,Vel{5})
plot(VelTime,Vel{6})
plot(VelTime,Vel{7})
xline(3)
grid on
xlabel('Time [s]')
ylabel('Encoder 4 Velocity [Counts/s]')
legend('Reference','-10% W_n','+10% W_n','-10% Beta','+10% Beta','-10% K', '+10% K','Location','best')
title('Simulated Closed Loop Step Response (3V)')


%Target Values
for i = 1:7
    Target(i) = mean(Vel{i}(1:136)); %Mean of values before the dip
end

%Overshoot
for i = 1:7
    OS_Sim(i) = 100*(max(Vel{i})-Target(i))/Target(i);
end


%Plot
figure()
hold on
plot(x{1},y{1})
plot(x{2},y{2})
plot(x{3},y{3})
plot(x{4},y{4})
plot(x{5},y{5})
plot(x{6},y{6})
plot(x{7},y{7})
% yline(Target*1.02,'--r')
% yline(Target*1.25,'--k')
% yline(Target*0.98,'--r')
xline(3,'--k')
xlabel('Time [s]')
ylabel('Encoder 4 Position [Counts]')
legend('Reference','-10% W_n','+10% W_n','-10% Beta','+10% Beta','-10% K', '+10% K','Location','best')
title('Simulated Closed Loop Step Response (3V)')
grid on
