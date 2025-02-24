clear all; close all; clc;
%% Read Data

%Read Open Loop Data
    %Trial 1 (3V)
    fileID = fopen('OL_1.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{1} = cell2mat(data); % Convert to numeric matrix
    Sample{1} = data{1}(:,1);  
    Time{1} = data{1}(:,2);     
    Encoder2Pos{1} = data{1}(:,3);
    Encoder4Pos{1} = data{1}(:,4);
    ControlEffort2{1} = data{1}(:,5);

    %Trial 2 (3V)
    fileID = fopen('OL_2.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{2} = cell2mat(data); % Convert to numeric matrix
    Sample{2} = data{2}(:,1);  
    Time{2} = data{2}(:,2);     
    Encoder2Pos{2} = data{2}(:,3);
    Encoder4Pos{2} = data{2}(:,4);
    ControlEffort2{2} = data{2}(:,5);

    %Trial 3 (3V)
    fileID = fopen('OL_3.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{3} = cell2mat(data); % Convert to numeric matrix
    Sample{3} = data{3}(:,1);  
    Time{3} = data{3}(:,2);     
    Encoder2Pos{3} = data{3}(:,3);
    Encoder4Pos{3} = data{3}(:,4);
    ControlEffort2{3} = data{3}(:,5);

    %Trial 4 (3V)
    fileID = fopen('OL_4.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{4} = cell2mat(data); % Convert to numeric matrix
    Sample{4} = data{4}(:,1);  
    Time{4} = data{4}(:,2);     
    Encoder2Pos{4} = data{4}(:,3);
    Encoder4Pos{4} = data{4}(:,4);
    ControlEffort2{4} = data{4}(:,5);

    %Trial 5 (3V)
    fileID = fopen('OL_5.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{5} = cell2mat(data); % Convert to numeric matrix
    Sample{5} = data{5}(:,1);  
    Time{5} = data{5}(:,2);     
    Encoder2Pos{5} = data{5}(:,3);
    Encoder4Pos{5} = data{5}(:,4);
    ControlEffort2{5} = data{5}(:,5);

%% Plot Open Loop

%Plot Open Loop Data
figure()
hold on
plot(Time{1},Encoder2Pos{1})
plot(Time{2},Encoder2Pos{2})
plot(Time{3},Encoder2Pos{3})
plot(Time{4},Encoder2Pos{4})
plot(Time{5},Encoder2Pos{5})
xline(3,'--k')
grid on
xlabel('Time [s]')
ylabel('Encoder 2 Position [counts]')
title('3000ms Dwell Time')
legend('Trial 1','Trial 2', 'Trial 3','Trial 4', 'Trial 5','Location','best')

% %Average Line (Error Analysis) - Position (All 5 Trials)
% for i = 1:numel(Time{1})
%     %Sample Mean
%         SM(i) = (1/5)*(Encoder2Pos{1}(i) + Encoder2Pos{2}(i) + Encoder2Pos{3}(i) + Encoder2Pos{4}(i) + Encoder2Pos{5}(i)); %Sample mean at each point
%     %Sample Standard Deviation
%         SSD(i) = sqrt(((Encoder2Pos{1}(i)-SM(i))^2+(Encoder2Pos{2}(i)-SM(i))^2+(Encoder2Pos{3}(i)-SM(i))^2+(Encoder2Pos{4}(i)-SM(i))^2+(Encoder2Pos{5}(i)-SM(i))^2)/5);
%     %Standard Deviation Bounds
%         Upper(i) = SM(i) + SSD(i);
%         Lower(i) = SM(i) - SSD(i);
% end
% 
% %Position Error Plot
% figure()
% hold on
% plot(Time{1},SM)
% plot(Time{1},Upper,'--r')
% plot(Time{1},Lower,'--r')
% xline(3,'--k')
% grid on
% xlabel('Time [s]')
% ylabel('Encoder 2 Position [Counts]')
% title('3000ms Dwell Time')
% legend('Average','Stand. Dev. Bounds')

%Average Line (Error Analysis) - Position (Excluding Trial 1)
for i = 1:numel(Time{1})
    %Sample Mean
        SM(i) = (1/4)*(Encoder2Pos{2}(i) + Encoder2Pos{3}(i) + Encoder2Pos{4}(i) + Encoder2Pos{5}(i)); %Sample mean at each point
    %Sample Standard Deviation
        SSD(i) = sqrt(((Encoder2Pos{2}(i)-SM(i))^2+(Encoder2Pos{3}(i)-SM(i))^2+(Encoder2Pos{4}(i)-SM(i))^2+(Encoder2Pos{5}(i)-SM(i))^2)/4);
    %Standard Deviation Bounds
        Upper(i) = SM(i) + SSD(i);
        Lower(i) = SM(i) - SSD(i);
end

%Position Error Plot
figure()
hold on
plot(Time{1},SM)
plot(Time{1},Upper,'--r')
plot(Time{1},Lower,'--r')
xline(3,'--k')
grid on
xlabel('Time [s]')
ylabel('Encoder 2 Position [Counts]')
title('3000ms Dwell Time')
legend('Average','Stand. Dev. Bounds')

%% Open Loop Model Validation
fig = openfig('Week2Validation.fig');
grid on
xlabel('Time [s]')
ylabel('Encoder 2 Position [Counts]')

%% Error

%Local Maxima
for i = 1:5
    LM{i} = islocalmax(Encoder2Pos{i}); %Locate local maxima
    Peaks{i} = Encoder2Pos{i}(LM{i}); %Find y-axis value of maxima
    Peak_Time{i} = Time{i}(LM{i}); %Find x-axis value of maxima
end

%y-infinity
for i = 1:5
    Reference(i) = Encoder2Pos{i}(68); %Find steady state value right before dip
end

%Omega_d // BOm_n
n = 2;
for i = 1:5
    Omega_d(i) = 2*pi*(n/(Peak_Time{i}(n+1)-Peak_Time{i}(1)));
    BOm_n(i) = (1/(Peak_Time{i}(n+1)-Peak_Time{i}(1)))*log((Peaks{i}(1)-Reference(i))/(Peaks{i}(n+1)-Reference(i)));
end

%Reporting Parameters
Voltage = 3; %Step input
for i = 1:5
    Omega_n(i) = sqrt(Omega_d(i)^2 + BOm_n(i)^2);
    Beta(i) = BOm_n(i)/Omega_n(i);
    K(i) = Reference(i)/Voltage;
end

%Omega_n Error
Mean_Omega_n = sum(Omega_n)/numel(Omega_n);
SD_Omega_n = sqrt(((Omega_n(1)-Mean_Omega_n)^2 + (Omega_n(2)-Mean_Omega_n)^2 + (Omega_n(3)-Mean_Omega_n)^2 + (Omega_n(4)-Mean_Omega_n)^2 + (Omega_n(5)-Mean_Omega_n)^2)/(numel(Omega_n)));
Conf_95_Bound_Omega_n = 1.96*SD_Omega_n;
Conf_95_Lower_Omega_n = Mean_Omega_n - Conf_95_Bound_Omega_n;
Conf_95_Upper_Omega_n = Mean_Omega_n + Conf_95_Bound_Omega_n;

%Beta Error
Mean_Beta = sum(Beta)/numel(Beta);
SD_Beta = sqrt(((Beta(1)-Mean_Beta)^2 + (Beta(2)-Mean_Beta)^2 + (Beta(3)-Mean_Beta)^2 + (Beta(4)-Mean_Beta)^2 + (Beta(5)-Mean_Beta)^2)/(numel(Beta)));
Conf_95_Bound_Beta = 1.96*SD_Beta;
Conf_95_Lower_Beta = Mean_Beta - Conf_95_Bound_Beta;
Conf_95_Upper_Beta = Mean_Beta + Conf_95_Bound_Beta;

%K Error
Mean_K = sum(K)/numel(K);
SD_K = sqrt(((K(1)-Mean_K)^2 + (K(2)-Mean_K)^2 + (K(3)-Mean_K)^2 + (K(4)-Mean_K)^2 + (K(5)-Mean_K)^2)/(numel(K)));
Conf_95_Bound_K = 1.96*SD_K;
Conf_95_Lower_K = Mean_K - Conf_95_Bound_K;
Conf_95_Upper_K = Mean_K + Conf_95_Bound_K;

%Transfer Function
G22 = tf(Mean_K*Mean_Omega_n^2,[1, 2*Mean_Beta*Mean_Omega_n, Mean_Omega_n^2]);


%% Closed Loop
clear all;

%Open Control Figure
fig = openfig('Wk2_CL.fig');

%Read Open Loop Data
    %Trial 1 (1.5V)
    fileID = fopen('KP017_KD0008.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{1} = cell2mat(data); % Convert to numeric matrix
    Sample{1} = data{1}(:,1);  
    Time{1} = data{1}(:,2);     
    Encoder2Pos{1} = data{1}(:,3);
    Encoder4Pos{1} = data{1}(:,4);
    ControlEffort2{1} = data{1}(:,5);

    %Trial 2 (1.5V)
    fileID = fopen('KP017_KD0008_2.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{2} = cell2mat(data); % Convert to numeric matrix
    Sample{2} = data{2}(:,1);  
    Time{2} = data{2}(:,2);     
    Encoder2Pos{2} = data{2}(:,3);
    Encoder4Pos{2} = data{2}(:,4);
    ControlEffort2{2} = data{2}(:,5);

    %Trial 3 (1.5V)
    fileID = fopen('KP017_KD0008_3.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{3} = cell2mat(data); % Convert to numeric matrix
    Sample{3} = data{3}(:,1);  
    Time{3} = data{3}(:,2);     
    Encoder2Pos{3} = data{3}(:,3);
    Encoder4Pos{3} = data{3}(:,4);
    ControlEffort2{3} = data{3}(:,5);

    %Trial 4 (1.5V)
    fileID = fopen('KP017_KD0008_4.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{4} = cell2mat(data); % Convert to numeric matrix
    Sample{4} = data{4}(:,1);  
    Time{4} = data{4}(:,2);     
    Encoder2Pos{4} = data{4}(:,3);
    Encoder4Pos{4} = data{4}(:,4);
    ControlEffort2{4} = data{4}(:,5);

    %Trial 5 (1.5V)
    fileID = fopen('KP017_KD0008_5.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{5} = cell2mat(data); % Convert to numeric matrix
    Sample{5} = data{5}(:,1);  
    Time{5} = data{5}(:,2);     
    Encoder2Pos{5} = data{5}(:,3);
    Encoder4Pos{5} = data{5}(:,4);
    ControlEffort2{5} = data{5}(:,5);

%Overshoot
for i = 1:5
    Reference(i) = Encoder2Pos{i}(68);
    Overshoot(i) = 100*(max(Encoder2Pos{i}-Reference(i)))/Reference(i);
end

OS_Mean = mean(Overshoot)
OS_SD = std(Overshoot)


%Plot Closed Loop Data
figure()
hold on
plot(Time{1},Encoder2Pos{1})
plot(Time{2},Encoder2Pos{2})
plot(Time{3},Encoder2Pos{3})
plot(Time{4},Encoder2Pos{4})
plot(Time{5},Encoder2Pos{5})
xline(3,'--k')
grid on
xlabel('Time [s]')
ylabel('Encoder 2 Position [counts]')
title('3000ms Dwell Time')
legend('Trial 1','Trial 2', 'Trial 3','Trial 4', 'Trial 5','Location','best')

%Average of Trials
for i = 1:numel(Time{1})
    %Sample Mean
        SM(i) = (1/5)*(Encoder2Pos{1}(i) + Encoder2Pos{2}(i) + Encoder2Pos{3}(i) + Encoder2Pos{4}(i) + Encoder2Pos{5}(i)); %Sample mean at each point
    %Sample Standard Deviation
        SSD(i) = sqrt(((Encoder2Pos{1}(i)-SM(i))^2+(Encoder2Pos{2}(i)-SM(i))^2+(Encoder2Pos{3}(i)-SM(i))^2+(Encoder2Pos{4}(i)-SM(i))^2+(Encoder2Pos{5}(i)-SM(i))^2)/5);
    %Standard Deviation Bounds
        Upper(i) = SM(i) + SSD(i);
        Lower(i) = SM(i) - SSD(i);
end

%Plot the averaged trials
figure()
hold on
plot(Time{1},SM)
plot(Time{1},Upper,'--r')
plot(Time{1},Lower,'--r')
xline(3,'--k')
grid on
xlabel('Time [s]')
ylabel('Encoder 2 Position [Counts]')
title('3000ms Dwell Time')
legend('Average','Stand. Dev. Bounds')

%Average Plot Overshoot
    Reference(6) = SM(68);
    Overshoot(6) = 100*(max(SM)-Reference(6))/Reference(6);

%% Sensitivity Analysis

%% Sensitivity Analysis

%Extract data from maelab.m figures
fig1 = openfig('Wk2_CL_Reference.fig');
    axObjs = fig1.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{1} = dataObjs(1).XData; %time
    y{1} = dataObjs(1).YData; %encoder position
fig2 = openfig('Wk2_CL_Wn_09.fig');
    axObjs = fig2.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{2} = dataObjs(1).XData;
    y{2} = dataObjs(1).YData;
fig3 = openfig('Wk2_CL_Wn_11.fig');
    axObjs = fig3.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{3} = dataObjs(1).XData;
    y{3} = dataObjs(1).YData;
fig4 = openfig('Wk2_CL_Beta_09.fig');
    axObjs = fig4.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{4} = dataObjs(1).XData;
    y{4} = dataObjs(1).YData;
fig5 = openfig('Wk2_CL_Beta_11.fig');
    axObjs = fig5.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{5} = dataObjs(1).XData;
    y{5} = dataObjs(1).YData;
fig6 = openfig('Wk2_CL_K_09.fig');
    axObjs = fig6.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{6} = dataObjs(1).XData;
    y{6} = dataObjs(1).YData;
fig7 = openfig('Wk2_CL_K_11.fig');
    axObjs = fig7.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{7} = dataObjs(1).XData;
    y{7} = dataObjs(1).YData;

%Overshoot and Settling Time
Target = 3500; %target value for each case (settled value)
Threshold = 0; %Initialize
for j = 1:7
    %Overshoot
    Target(j) = y{j}(675); %At about 2.98 seconds, so should be the settled value
    OS_Sim(j) = 100*(max(y{j})-Target(j))/Target(j);

    %Settling Time - 2%
    for i = 1:numel(y{j})
    check = y{j}(i)-Target(j);
    if 100*abs(check)/Target(j) <= 2
        if 100*(y{j}(i+5)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+10)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+25)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+35)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+45)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+50)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+75)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+100)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+125)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+150)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+200)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+225)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+300)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+350)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+400)-Target(j))/Target(j) <= 2
            Threshold = Threshold + 1;
        end
        if Threshold == 15
            break
        end
    end
    Threshold = 0;
end

ST_CL(j) = x{j}(i);
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
ylabel('Encoder 2 Position [Counts]')
legend('Reference','-10% W_n','+10% W_n','-10% Beta','+10% Beta','-10% K', '+10% K')
title('Simulated Closed Loop Step Response (1.5V)')
grid on