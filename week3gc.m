% Week 3: Kp=-0.05 Kd=0 Ki=-0.075
% Load the data from the text file

data_3 = readmatrix('KPN005_KIN0075.txt');

% Extract data points
time = data_3(:,1);
% encoder_4 = data_3(:,2);
control_2 = data_3(:,3);

% Plot both over time
figure;
% plot(time, encoder_4, 'b', 'LineWidth', 1.5);
hold on;
plot(time, control_2, 'r', 'LineWidth', 1.5);
hold off;

% Labels and title
xlabel('Time');
ylabel('Values');
title('Control 2 over Time');
legend('Control 2');
grid on;