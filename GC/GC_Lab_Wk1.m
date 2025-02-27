clear all; close all; clc;
%% Week 1

%Read Open Loop Data
    %2V Open Loop
    fileID = fopen('GC_2000c_4000ms.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{1} = cell2mat(data); % Convert to numeric matrix
    Sample{1} = data{1}(:,1);  
    Time{1} = data{1}(:,2);     
    Encoder3Pos{1} = data{1}(:,3);  
    ControlEffort1{1} = data{1}(:,4);

    %1.9V Open Loop
    fileID = fopen('GC_1900c_4000ms.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{2} = cell2mat(data); % Convert to numeric matrix
    Sample{2} = data{2}(:,1);  
    Time{2} = data{2}(:,2);     
    Encoder3Pos{2} = data{2}(:,3);  
    ControlEffort1{2} = data{2}(:,4);

    %1.8V Open Loop
    fileID = fopen('GC_1800c_4000ms.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{3} = cell2mat(data); % Convert to numeric matrix
    Sample{3} = data{3}(:,1);  
    Time{3} = data{3}(:,2);     
    Encoder3Pos{3} = data{3}(:,3);  
    ControlEffort1{3} = data{3}(:,4);

    %1.7V Open Loop
    fileID = fopen('GC_1700c_4000ms.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{4} = cell2mat(data); % Convert to numeric matrix
    Sample{4} = data{4}(:,1);  
    Time{4} = data{4}(:,2);     
    Encoder3Pos{4} = data{4}(:,3);  
    ControlEffort1{4} = data{4}(:,4); 

    %1.6V Open Loop
    fileID = fopen('GC_1600c_4000ms.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{5} = cell2mat(data); % Convert to numeric matrix
    Sample{5} = data{5}(:,1);  
    Time{5} = data{5}(:,2);     
    Encoder3Pos{5} = data{5}(:,3);  
    ControlEffort1{5} = data{5}(:,4);

%Plot Open Loop Data - Position
figure()
hold on
plot(Time{1},Encoder3Pos{1})
plot(Time{2},Encoder3Pos{2})
plot(Time{3},Encoder3Pos{3})
plot(Time{4},Encoder3Pos{4})
plot(Time{5},Encoder3Pos{5})
xline(4,'--k')
grid on
xlabel('Time [s]')
ylabel('Encoder 3 Position [counts]')
title('4000ms Dwell Time')
legend('Trial 1 (2 Volts)','Trial 2 (1.9 Volts)', 'Trial 3 (1.8 Volts)','Trial 4 (1.7 Volts)', 'Trial 5 (1.6 Volts)','Location','best')

%Average Line (Error Analysis) - Position
for i = 1:numel(Time{1})
    %Sample Mean
        SM(i) = (1/5)*(Encoder3Pos{1}(i) + Encoder3Pos{2}(i) + Encoder3Pos{3}(i) + Encoder3Pos{4}(i) + Encoder3Pos{5}(i)); %Sample mean at each point
    %Sample Standard Deviation
        SSD(i) = sqrt(((Encoder3Pos{1}(i)-SM(i))^2+(Encoder3Pos{2}(i)-SM(i))^2+(Encoder3Pos{3}(i)-SM(i))^2+(Encoder3Pos{4}(i)-SM(i))^2+(Encoder3Pos{5}(i)-SM(i))^2)/5);
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
xline(4,'--k')
grid on
xlabel('Time [s]')
ylabel('Encoder 3 Position [Counts]')
title('4000ms Dwell Time')
legend('Average','Stand. Dev. Bounds')

%Plot Open Loop Data - Control Effort
figure()
hold on
plot(Time{1},ControlEffort1{1})
plot(Time{2},ControlEffort1{2})
plot(Time{3},ControlEffort1{3})
plot(Time{4},ControlEffort1{4})
plot(Time{5},ControlEffort1{5})
xline(4,'--k')
grid on
xlabel('Time [s]')
ylabel('Control Effort 1 [counts]')
title('4000ms Dwell Time')
legend('Trial 1 (2 Volts)','Trial 2 (1.9 Volts)', 'Trial 3 (1.8 Volts)','Trial 4 (1.7 Volts)', 'Trial 5 (1.6 Volts)','Location','best')

%Average Line (Error Analysis) - Position
for i = 1:numel(Time{1})
    %Sample Mean
        SMC(i) = (1/5)*(ControlEffort1{1}(i) + ControlEffort1{2}(i) + ControlEffort1{3}(i) + ControlEffort1{4}(i) + ControlEffort1{5}(i)); %Sample mean at each point
    %Sample Standard Deviation
        SSDC(i) = sqrt(((ControlEffort1{1}(i)-SMC(i))^2+(ControlEffort1{2}(i)-SMC(i))^2+(ControlEffort1{3}(i)-SMC(i))^2+(ControlEffort1{4}(i)-SMC(i))^2+(ControlEffort1{5}(i)-SMC(i))^2)/5);
    %Standard Deviation Bounds
        UpperC(i) = SMC(i) + SSDC(i);
        LowerC(i) = SMC(i) - SSDC(i);
end

%Position Error Plot - Control
figure()
hold on
plot(Time{1},SMC)
plot(Time{1},UpperC,'--r')
plot(Time{1},LowerC,'--r')
xline(4,'--k')
grid on
xlabel('Time [s]')
ylabel('Control Effort 1 [Counts]')
title('4000ms Dwell Time')
legend('Average','Stand. Dev. Bounds')

%Velocity Calculation
i=1;
k=1;
n = 10;
while i < 900
%for i = 1:numel(SM)-1
    VelTime(k) = (Time{1}(i+n)+Time{1}(i))/2;
    SMV(k) = (SM(i+n)-SM(i))/(Time{1}(i+n)-Time{1}(i));
    Tr1_Vel(k) = (Encoder3Pos{1}(i+n)-Encoder3Pos{1}(i))/(Time{1}(i+n)-Time{1}(i));
    Tr2_Vel(k) = (Encoder3Pos{2}(i+n)-Encoder3Pos{2}(i))/(Time{1}(i+n)-Time{1}(i));
    Tr3_Vel(k) = (Encoder3Pos{3}(i+n)-Encoder3Pos{3}(i))/(Time{1}(i+n)-Time{1}(i));
    Tr4_Vel(k) = (Encoder3Pos{4}(i+n)-Encoder3Pos{4}(i))/(Time{1}(i+n)-Time{1}(i));
    Tr5_Vel(k) = (Encoder3Pos{5}(i+n)-Encoder3Pos{5}(i))/(Time{1}(i+n)-Time{1}(i));
    k = k+1;
    i = i+n;
end

%Average Line (Error Analysis)
for i = 1:numel(SMV)
    %Sample Mean
        SMV(i) = (1/5)*(Tr1_Vel(i) + Tr2_Vel(i) + Tr3_Vel(i) + Tr4_Vel(i) + Tr5_Vel(i)); %Sample mean at each point
    %Sample Standard Deviation
        SSDV(i) = sqrt(((Tr1_Vel(i)-SMV(i))^2+(Tr2_Vel(i)-SMV(i))^2+(Tr3_Vel(i)-SMV(i))^2+(Tr4_Vel(i)-SMV(i))^2+(Tr5_Vel(i)-SMV(i))^2)/5);
    %Standard Deviation Bounds
        UpperV(i) = SMV(i) + SSDV(i);
        LowerV(i) = SMV(i) - SSDV(i);
end

%Velocity Plot
figure()
hold on
plot(VelTime,SMV)
plot(VelTime,UpperV,'--r')
plot(VelTime,LowerV,'--r')
grid on
xlabel('Time [s]')
ylabel('Encoder 3 Velocity [Counts/s]')
legend('Average','Stand. Dev. Bounds','Location','best')

%Acceleration Calculation
k=1;
j = 1;
m = 1;
while j < i-m
    AccTime(k) = (VelTime(j+m)+VelTime(j))/2;
    SMA(k) = (SMV(j+m)-SMV(j))/(VelTime(j+m)-VelTime(j));
    Tr1_Acc(k) = (Tr1_Vel(j+m)-Tr1_Vel(j))/(Time{1}(j+m)-VelTime(j));
    Tr2_Acc(k) = (Tr2_Vel(j+m)-Tr2_Vel(j))/(VelTime(j+m)-VelTime(j));
    Tr3_Acc(k) = (Tr3_Vel(j+m)-Tr3_Vel(j))/(VelTime(j+m)-VelTime(j));
    Tr4_Acc(k) = (Tr4_Vel(j+m)-Tr4_Vel(j))/(VelTime(j+m)-VelTime(j));
    Tr5_Acc(k) = (Tr5_Vel(j+m)-Tr5_Vel(j))/(VelTime(j+m)-VelTime(j));
    k = k+1;
    j = j+m;
end

%Average Line (Error Analysis)
for i = 1:numel(SMA)
    %Sample Mean
        SMA(i) = (1/5)*(Tr1_Acc(i) + Tr2_Acc(i) + Tr3_Acc(i) + Tr4_Acc(i) + Tr5_Acc(i)); %Sample mean at each point
    %Sample Standard Deviation
        SSDA(i) = sqrt(((Tr1_Acc(i)-SMA(i))^2+(Tr2_Acc(i)-SMA(i))^2+(Tr3_Acc(i)-SMA(i))^2+(Tr4_Acc(i)-SMA(i))^2+(Tr5_Acc(i)-SMA(i))^2)/5);
    %Standard Deviation Bounds
        UpperA(i) = SMA(i) + SSDA(i);
        LowerA(i) = SMA(i) - SSDA(i);
end

%Acceleration Plot
figure()
hold on
plot(AccTime,SMA)
plot(AccTime,UpperA,'--r')
plot(AccTime,LowerA,'--r')
grid on
xlabel('Time [s]')
ylabel('Encoder 3 Acceleration [Counts/s^2]')
legend('Average','Stand. Dev. Bounds','Location','best')

%Average Parameter Calculations
(SMV(numel(SMV)/2-5)-SMV(1))/(VelTime(numel(SMV)/2-5)-VelTime(1));

%% Open Loop Error (Parameters)

%Velocity Slope at t=0
    dy_dt(1) = (Tr1_Vel(10)-Tr1_Vel(1))/(VelTime(10)-VelTime(1));
    dy_dt(2) = (Tr2_Vel(10)-Tr2_Vel(1))/(VelTime(10)-VelTime(1));
    dy_dt(3) = (Tr3_Vel(10)-Tr3_Vel(1))/(VelTime(10)-VelTime(1));
    dy_dt(4) = (Tr4_Vel(10)-Tr4_Vel(1))/(VelTime(10)-VelTime(1));
    dy_dt(5) = (Tr5_Vel(10)-Tr5_Vel(1))/(VelTime(10)-VelTime(1));
dy_dt_mean = sum(dy_dt)/numel(dy_dt);
dy_dt_SD = sqrt(((dy_dt(1)-dy_dt_mean)^2 + (dy_dt(2)-dy_dt_mean)^2 + (dy_dt(3)-dy_dt_mean)^2 + (dy_dt(4)-dy_dt_mean)^2 + (dy_dt(5)-dy_dt_mean)^2)/(numel(dy_dt)));
dy_dt_Conf_95_Bound = 1.96*dy_dt_SD;
dy_dt_Conf_95_Lower = dy_dt_mean - dy_dt_Conf_95_Bound;
dy_dt_Conf_95_Upper = dy_dt_mean + dy_dt_Conf_95_Bound;

%Gain K
    K(1) = dy_dt(1)*(1/2);
    K(2) = dy_dt(2)*(1/1.9);
    K(3) = dy_dt(3)*(1/1.8);
    K(4) = dy_dt(4)*(1/1.7);
    K(5) = dy_dt(5)*(1/1.6);
K_mean = sum(K)/numel(K);
K_SD = sqrt(((K(1)-K_mean)^2 + (K(2)-K_mean)^2 + (K(3)-K_mean)^2 + (K(4)-K_mean)^2 + (K(5)-K_mean)^2)/(numel(K)));
K_Conf_95_Bound = 1.96*K_SD;
K_Conf_95_Lower = K_mean - K_Conf_95_Bound;
K_Conf_95_Upper = K_mean + K_Conf_95_Bound;

%Beta
    Beta = [0.1,0.02,0.05,0.02,0.125];
Beta_mean = mean(Beta);
Beta_SD = sqrt(((Beta(1)-Beta_mean)^2 + (Beta(2)-Beta_mean)^2 + (Beta(3)-Beta_mean)^2 + (Beta(4)-Beta_mean)^2 + (Beta(5)-Beta_mean)^2)/(numel(Beta)))
Beta_Conf_95_Bound = 1.96*Beta_SD
Beta_Conf_95_Lower = Beta_mean - Beta_Conf_95_Bound
Beta_Conf_95_Upper = Beta_mean + Beta_Conf_95_Bound

%% Closed Loop

%File Read
    %3.5V Closed Loop (With Controller)
    fileID = fopen('KPN08_KDN01.txt', 'r'); % Open file for reading
    fgetl(fileID); % Skip the first line (header)
    data = textscan(fileID, '%f %f %f %f', 'Delimiter', {' ', '\t', ';','[',']'}, 'MultipleDelimsAsOne', true);
    fclose(fileID); % Close the file
    data{6} = cell2mat(data); % Convert to numeric matrix
    Sample{2} = data{6}(:,1);  
    Time{6} = data{6}(:,2);     
    Encoder3Pos{6} = data{6}(:,3);  
    ControlEffort1{6} = data{6}(:,4);

%Determine Overshoot Bounds (47-90)
Settled = mean(Encoder3Pos{6}(47:90));
OS_Bound = 1.25*Settled;

%Determine Overshoot
Overshoot = 100*(max(Encoder3Pos{6})-Settled)/Settled; %percent
fprintf('The overshoot is %d percent\n',Overshoot)

%Settling Time - 2%
Threshold = 0;
for i = 1:numel(Encoder3Pos{6})
    check = Encoder3Pos{6}(i)-Settled;
    if 100*abs(check)/Settled <= 2
        if 100*(Encoder3Pos{6}(i+1)-Settled)/Settled <= 2
            Threshold = Threshold + 1;
        end
        if 100*(Encoder3Pos{6}(i+2)-Settled)/Settled <= 2
            Threshold = Threshold + 1;
        end
        if 100*(Encoder3Pos{6}(i+3)-Settled)/Settled <= 2
            Threshold = Threshold + 1;
        end
        if Threshold == 3
            break
        end
    end
    Threshold = 0;
end

Settling_Time = Time{6}(i);
fprintf('The settling time is %d seconds\n',Settling_Time)

%Plot Closed Loop Data - Position
figure()
hold on
plot(Time{6},Encoder3Pos{6})
xline(4,'--k')
yline(OS_Bound,'--r')
grid on
xlabel('Time [s]')
ylabel('Encoder 3 Position [counts]')
title('4000ms Dwell Time')
legend('Measured Response', 'Control Input Cutoff', '25% Overshoot Bound','Location','best')

%Plot Closed Loop Data
figure()
hold on
plot(Time{6},ControlEffort1{6})
xline(4,'--k')
grid on
xlabel('Time [s]')
ylabel('Control Effort 1 [counts]')
title('4000ms Dwell Time')
legend('Closed Loop Control Effort','Location','best')

%Figure
openfig('Week1PDSIM.fig')

%% Sensitivity Analysis

%Extract data from maelab.m figures
fig1 = openfig('K09_CL.fig');
    axObjs = fig1.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{1} = dataObjs(1).XData;
    y{1} = dataObjs(1).YData;
fig2 = openfig('K11_CL.fig');
    axObjs = fig2.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{2} = dataObjs(1).XData;
    y{2} = dataObjs(1).YData;
fig3 = openfig('B09_CL.fig');
    axObjs = fig3.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{3} = dataObjs(1).XData;
    y{3} = dataObjs(1).YData;
fig4 = openfig('B11_CL.fig');
    axObjs = fig4.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{4} = dataObjs(1).XData;
    y{4} = dataObjs(1).YData;
fig5 = openfig('CL_Reference.fig');
    axObjs = fig5.Children;
    topAx = axObjs(3);
    dataObjs = topAx.Children; % Get plot data
    x{5} = dataObjs(1).XData;
    y{5} = dataObjs(1).YData;

%Overshoot and Settling Time
Target = 3500; %target value for each case (settled value)
Threshold = 0; %Initialize
for j = 1:5
    %Overshoot
    OS_Sim(j) = 100*(max(y{j})-Target)/Target;

    %Settling Time - 2%
    for i = 1:numel(y{j})
    check = y{j}(i)-Target;
    if 100*abs(check)/Target <= 2
        if 100*(y{j}(i+10)-Target)/Target <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+75)-Target)/Target <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+150)-Target)/Target <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+225)-Target)/Target <= 2
            Threshold = Threshold + 1;
        end
        if 100*(y{j}(i+300)-Target)/Target <= 2
            Threshold = Threshold + 1;
        end
        if Threshold == 5
            break
        end
    end
    Threshold = 0;
end

ST_CL(j) = x{j}(i);
end


figure()
hold on
plot(x{1},y{1})
plot(x{2},y{2})
plot(x{3},y{3})
plot(x{4},y{4})
plot(x{5},y{5})
yline(Target*1.02,'--r')
yline(Target*1.25,'--k')
yline(Target*0.98,'--r')
xlabel('Time [s]')
ylabel('Encoder 3 Position [Counts]')
legend('-10% K','+10% K','-10% Beta','+10% Beta','Reference','2% Settle Bounds','25% Overshoot Limit','FontSize',10)
title('Simulated Closed Loop Step Response (3.5V)')
grid on
