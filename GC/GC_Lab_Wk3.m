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
