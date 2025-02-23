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

%Average Line (Error Analysis) - Position
for i = 1:numel(Time{1})
    %Sample Mean
        SM(i) = (1/5)*(Encoder2Pos{1}(i) + Encoder2Pos{2}(i) + Encoder2Pos{3}(i) + Encoder2Pos{4}(i) + Encoder2Pos{5}(i)); %Sample mean at each point
    %Sample Standard Deviation
        SSD(i) = sqrt(((Encoder2Pos{1}(i)-SM(i))^2+(Encoder2Pos{2}(i)-SM(i))^2+(Encoder2Pos{3}(i)-SM(i))^2+(Encoder2Pos{4}(i)-SM(i))^2+(Encoder2Pos{5}(i)-SM(i))^2)/5);
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

%Average Line (Error Analysis) - Position
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

%% Closed Loop
openfig('Wk2_CL.fig')