% Three speeds (U): 20, 35, 50 [m/s]
% Five AoA (alp): 0, 4, 8, 12, 16 [deg]

% Equations
% Cp = (rho - rho_inf) / q_inf

% Measured Constants
rho = 1.1963; % [kg/m^3] air density

% Assumed Constants
mu = 18.13e-06; % [Pa-s]

% Channel Locations Matrix
loc = [0 7.5 10 20 30 40 50 60 70 80]/100*3.5; % Total chord is 3.5"
loc_u = loc;
loc_l = loc(1:end-1);
chord_spread = linspace(0,3.5,1000);

% Chord Length Positions
chord = linspace(0,3.5,1000);

% Velocity Matrix
U_matrix = [20 35 50];

% Plot Axis Limit [Turned off]
axis_limit = [0 3.5 -4*10^-3 4*10^-3];

% Rake Probe Locations
rake_spread = linspace(0,1.7,18)-1.7/2; % [in];

%% Notes

% Unfortnuately, you have to manually import each column vector.

% Each column vector is denoted by AOA and then speed.
% First index corresponds to channel 0 (Leading Edge).
% Indices 2 through 10 correspond to channels 1 through 9 (Upper).
% Indices 11 through 18 correspond to channels 10 through 17 (Lower).

% For channel positions, refer to the appendix of the laboratory handout.

%% 0 Degrees AOA
load('zero_deg_0.mat'), load('zero_deg_20.mat'), load('zero_deg_35.mat'), load('zero_deg_50.mat')

% Calibration (Zeroing)
zero_deg_20 = (zero_deg_20 - zero_deg_0)*248.84;
zero_deg_35 = (zero_deg_35 - zero_deg_0)*248.84;
zero_deg_50 = (zero_deg_50 - zero_deg_0)*248.84;


% Generating Upper & Lower Pressure Matrices
zero_deg_20_u = zero_deg_20(2:10);
zero_deg_20_l = zero_deg_20(11:18);

zero_deg_35_u = zero_deg_35(2:10);
zero_deg_35_l = zero_deg_35(11:18);

zero_deg_50_u = zero_deg_50(2:10);
zero_deg_50_l = zero_deg_50(11:18);

% Average Pressure Matrix
zero_deg_20_avg = (zero_deg_20_u(1:end-1) - zero_deg_20_l)/2;
zero_deg_20_avg = [zero_deg_20(1)' zero_deg_20_avg(1:8)' zero_deg_20(10)' zero_deg_20_avg(9:end)'];

zero_deg_35_avg = (zero_deg_35_u(1:end-1) - zero_deg_35_l)/2;
zero_deg_35_avg = [zero_deg_35(1)' zero_deg_35_avg(1:8)' zero_deg_35(10)' zero_deg_35_avg(9:end)'];

zero_deg_50_avg = (zero_deg_50_u(1:end-1) - zero_deg_50_l)/2;
zero_deg_50_avg = [zero_deg_50(1)' zero_deg_50_avg(1:8)' zero_deg_50(10)' zero_deg_50_avg(9:end)'];

% Calculating and Graphing Cp
Cp_zero_deg_20 = zero_deg_20_avg / zero_deg_20(end);
Cp_zero_deg_35 = zero_deg_35_avg / zero_deg_35(end);
Cp_zero_deg_50 = zero_deg_50_avg / zero_deg_50(end);

figure(1)
plot(loc,Cp_zero_deg_20,'-r')
hold on
plot(loc,Cp_zero_deg_35,'-b')
plot(loc,Cp_zero_deg_50,'-g')
hold off
grid('on')
%axis(axis_limit)
legend('20 m/s', '35 m/s', '50 m/s')
title('Cp Along Chord at 0 AOA and Various Speeds')

%% 4 Degrees AOA
load('four_deg_0.mat'), load('four_deg_20.mat'), load('four_deg_35.mat'), load('four_deg_50.mat')

% Calibration (Zeroing)
four_deg_20 = (four_deg_20 - four_deg_0)*248.84;
four_deg_35 = (four_deg_35 - four_deg_0)*248.84;
four_deg_50 = (four_deg_50 - four_deg_0)*248.84;


% Generating Upper & Lower Pressure Matrices
four_deg_20_u = four_deg_20(2:10);
four_deg_20_l = four_deg_20(11:18);

four_deg_35_u = four_deg_35(2:10);
four_deg_35_l = four_deg_35(11:18);

four_deg_50_u = four_deg_50(2:10);
four_deg_50_l = four_deg_50(11:18);

% Average Pressure Matrix
four_deg_20_avg = (four_deg_20_u(1:end-1) - four_deg_20_l)/2;
four_deg_20_avg = [four_deg_20(1)' four_deg_20_avg(1:8)' four_deg_20(10)' four_deg_20_avg(9:end)'];

four_deg_35_avg = (four_deg_35_u(1:end-1) - four_deg_35_l)/2;
four_deg_35_avg = [four_deg_35(1)' four_deg_35_avg(1:8)' four_deg_35(10)' four_deg_35_avg(9:end)'];

four_deg_50_avg = (four_deg_50_u(1:end-1) - four_deg_50_l)/2;
four_deg_50_avg = [four_deg_50(1)' four_deg_50_avg(1:8)' four_deg_50(10)' four_deg_50_avg(9:end)'];

% Calculating and Graphing Cp
Cp_four_deg_20 = four_deg_20_avg / four_deg_20(end);
Cp_four_deg_35 = four_deg_35_avg / four_deg_35(end);
Cp_four_deg_50 = four_deg_50_avg / four_deg_50(end);

figure(2)
plot(loc,Cp_four_deg_20,'-r')
hold on
plot(loc,Cp_four_deg_35,'-b')
plot(loc,Cp_four_deg_50,'-g')
hold off
grid('on')
%axis(axis_limit)
legend('20 m/s', '35 m/s', '50 m/s')
title('Cp Along Chord at 4 AOA and Various Speeds')

%% 8 Degrees AOA
load('eight_deg_0.mat'), load('eight_deg_20.mat'), load('eight_deg_35.mat'), load('eight_deg_50.mat')

% Calibration (Zeroing)
eight_deg_20 = (eight_deg_20 - eight_deg_0)*248.84;
eight_deg_35 = (eight_deg_35 - eight_deg_0)*248.84;
eight_deg_50 = (eight_deg_50 - eight_deg_0)*248.84;


% Generating Upper & Lower Pressure Matrices
eight_deg_20_u = eight_deg_20(2:10);
eight_deg_20_l = eight_deg_20(11:18);


eight_deg_35_u = eight_deg_35(2:10);
eight_deg_35_l = eight_deg_35(11:18);

eight_deg_50_u = eight_deg_50(2:10);
eight_deg_50_l = eight_deg_50(11:18);

% Average Pressure Matrix
eight_deg_20_avg = (eight_deg_20_u(1:end-1) - eight_deg_20_l)/2;
eight_deg_20_avg = [eight_deg_20(1)' eight_deg_20_avg(1:8)' eight_deg_20(10)' eight_deg_20_avg(9:end)'];

eight_deg_35_avg = (eight_deg_35_u(1:end-1) - eight_deg_35_l)/2;
eight_deg_35_avg = [eight_deg_35(1)' eight_deg_35_avg(1:8)' eight_deg_35(10)' eight_deg_35_avg(9:end)'];

eight_deg_50_avg = (eight_deg_50_u(1:end-1) - eight_deg_50_l)/2;
eight_deg_50_avg = [eight_deg_50(1)' eight_deg_50_avg(1:8)' eight_deg_50(10)' eight_deg_50_avg(9:end)'];

% Calculating and Graphing Cp
Cp_eight_deg_20 = eight_deg_20_avg / eight_deg_20(end);
Cp_eight_deg_35 = eight_deg_35_avg / eight_deg_35(end);
Cp_eight_deg_50 = eight_deg_50_avg / eight_deg_50(end);

figure(3)
plot(loc,Cp_eight_deg_20,'-r')
hold on
plot(loc,Cp_eight_deg_35,'-b')
plot(loc,Cp_eight_deg_50,'-g')
hold off
grid('on')
%axis(axis_limit)
legend('20 m/s', '35 m/s', '50 m/s', Location='SouthEast')
title('C_p Along Chord at 8 AOA and Various Speeds')

%% 12 Degrees AOA
load('twelve_deg_0.mat'), load('twelve_deg_20.mat'), load('twelve_deg_35.mat'), load('twelve_deg_50.mat')

% Calibration (Zeroing)
twelve_deg_20 = (twelve_deg_20 - twelve_deg_0)*248.84;
twelve_deg_35 = (twelve_deg_35 - twelve_deg_0)*248.84;
twelve_deg_50 = (twelve_deg_50 - twelve_deg_0)*248.84;


% Generating Upper & Lower Pressure Matrices
twelve_deg_20_u = twelve_deg_20(2:10);
twelve_deg_20_l = twelve_deg_20(11:18);

twelve_deg_35_u = twelve_deg_35(2:10);
twelve_deg_35_l = twelve_deg_35(11:18);

twelve_deg_50_u = twelve_deg_50(2:10);
twelve_deg_50_l = twelve_deg_50(11:18);

% Average Pressure Matrix
twelve_deg_20_avg = (twelve_deg_20_u(1:end-1) - twelve_deg_20_l)/2;
twelve_deg_20_avg = [twelve_deg_20(1)' twelve_deg_20_avg(1:8)' twelve_deg_20(10)' twelve_deg_20_avg(9:end)'];

twelve_deg_35_avg = (twelve_deg_35_u(1:end-1) - twelve_deg_35_l)/2;
twelve_deg_35_avg = [twelve_deg_35(1)' twelve_deg_35_avg(1:8)' twelve_deg_35(10)' twelve_deg_35_avg(9:end)'];

twelve_deg_50_avg = (twelve_deg_50_u(1:end-1) - twelve_deg_50_l)/2;
twelve_deg_50_avg = [twelve_deg_50(1)' twelve_deg_50_avg(1:8)' twelve_deg_50(10)' twelve_deg_50_avg(9:end)'];

% Calculating and Graphing Cp
Cp_twelve_deg_20 = twelve_deg_20_avg / twelve_deg_20(end);
Cp_twelve_deg_35 = twelve_deg_35_avg / twelve_deg_35(end);
Cp_twelve_deg_50 = twelve_deg_50_avg / twelve_deg_50(end);

figure(4)
plot(loc,Cp_twelve_deg_20,'-r')
hold on
plot(loc,Cp_twelve_deg_35,'-b')
plot(loc,Cp_twelve_deg_50,'-g')
hold off
grid('on')
%axis(axis_limit)
legend('20 m/s', '35 m/s', '50 m/s', Location='SouthEast')
title('C_p Along Chord at 12 AOA and Various Speeds')

%% 16 Degrees AOA
load('sixteen_deg_0.mat'), load('sixteen_deg_20.mat'), load('sixteen_deg_35.mat'), load('sixteen_deg_50.mat')

% Calibration (Zeroing)
sixteen_deg_20 = (sixteen_deg_20 - sixteen_deg_0)*248.84;
sixteen_deg_35 = (sixteen_deg_35 - sixteen_deg_0)*248.84;
sixteen_deg_50 = (sixteen_deg_50 - sixteen_deg_0)*248.84;


% Generating Upper & Lower Pressure Matrices
sixteen_deg_20_u = sixteen_deg_20(2:10);
sixteen_deg_20_l = sixteen_deg_20(11:18);

sixteen_deg_35_u = sixteen_deg_35(2:10);
sixteen_deg_35_l = sixteen_deg_35(11:18);

sixteen_deg_50_u = sixteen_deg_50(2:10);
sixteen_deg_50_l = sixteen_deg_50(11:18);

% Average Pressure Matrix
sixteen_deg_20_avg = (sixteen_deg_20_u(1:end-1) - sixteen_deg_20_l)/2;
sixteen_deg_20_avg = [sixteen_deg_20(1)' sixteen_deg_20_avg(1:8)' sixteen_deg_20(10)' sixteen_deg_20_avg(9:end)'];

sixteen_deg_35_avg = (sixteen_deg_35_u(1:end-1) - sixteen_deg_35_l)/2;
sixteen_deg_35_avg = [sixteen_deg_35(1)' sixteen_deg_35_avg(1:8)' sixteen_deg_35(10)' sixteen_deg_35_avg(9:end)'];

sixteen_deg_50_avg = (sixteen_deg_50_u(1:end-1) - sixteen_deg_50_l)/2;
sixteen_deg_50_avg = [sixteen_deg_50(1)' sixteen_deg_50_avg(1:8)' sixteen_deg_50(10)' sixteen_deg_50_avg(9:end)'];

% Calculating and Graphing Cp
dynamic_pressure = 0.5 * rho * U_matrix.^2;
Cp_sixteen_deg_20 = sixteen_deg_20_avg / sixteen_deg_20(end);
Cp_sixteen_deg_35 = sixteen_deg_35_avg / sixteen_deg_35(end);
Cp_sixteen_deg_50 = sixteen_deg_50_avg / sixteen_deg_50(end);

figure(5)
plot(loc,Cp_sixteen_deg_20,'-r')
hold on
plot(loc,Cp_sixteen_deg_35,'-b')
plot(loc,Cp_sixteen_deg_50,'-g')
hold off
grid('on')
%axis(axis_limit)
legend('20 m/s', '35 m/s', '50 m/s', Location='SouthEast')
title('C_p Along Chord at 16 AOA and Various Speeds')

%% Calculating Cp_u, Cp_l, Cn, and their fitted values

% Zero Deg AOA
Cp_u_zero_deg_20 = [zero_deg_20(1) zero_deg_20_u'] / zero_deg_20(end);
Cp_l_zero_deg_20 = [zero_deg_20(1) zero_deg_20_l'] / zero_deg_20(end);
cpu020_pval = polyfit(loc_u,Cp_u_zero_deg_20,2);
cpl020_pval = polyfit(loc_l,Cp_l_zero_deg_20,2);
cn020_pval = cpl020_pval - cpu020_pval;
Cp_u_zero_20_fit = polyval(cpu020_pval, chord_spread);
Cp_l_zero_20_fit = polyval(cpl020_pval, chord_spread);
%Cn_zero_20_fit = calculate_cn_values(cn020_pval,chord,chord_spread);
Cn_zero_20_fit = calc_integral(cn020_pval,3.5)/3.5;

Cp_u_zero_deg_35 = [zero_deg_35(1) zero_deg_35_u'] / zero_deg_35(end);
Cp_l_zero_deg_35 = [zero_deg_35(1) zero_deg_35_l'] / zero_deg_35(end);
cpu035_pval = polyfit(loc_u,Cp_u_zero_deg_35,2);
cpl035_pval = polyfit(loc_l,Cp_l_zero_deg_35,2);
cn035_pval = cpl035_pval - cpu035_pval;
Cp_u_zero_35_fit = polyval(cpu035_pval, chord_spread);
Cp_l_zero_35_fit = polyval(cpl035_pval, chord_spread);
%Cn_zero_35_fit = calculate_cn_values(cn035_pval,chord,chord_spread);
Cn_zero_35_fit = calc_integral(cn035_pval,3.5)/3.5;

Cp_u_zero_deg_50 = [zero_deg_50(1) zero_deg_50_u'] / zero_deg_50(end);
Cp_l_zero_deg_50 = [zero_deg_50(1) zero_deg_50_l'] / zero_deg_50(end);
cpu050_pval = polyfit(loc_u,Cp_u_zero_deg_50,2);
cpl050_pval = polyfit(loc_l,Cp_l_zero_deg_50,2);
cn050_pval = cpl050_pval - cpu050_pval;
Cp_u_zero_50_fit = polyval(cpu050_pval, chord_spread);
Cp_l_zero_50_fit = polyval(cpl050_pval, chord_spread);
%Cn_zero_50_fit = calculate_cn_values(cn050_pval,chord,chord_spread);
Cn_zero_50_fit = calc_integral(cn050_pval,3.5)/3.5;

% Four Deg AOA
Cp_u_four_deg_20 = [four_deg_20(1) four_deg_20_u'] / four_deg_20(end);
Cp_l_four_deg_20 = [four_deg_20(1) four_deg_20_l'] / four_deg_20(end);
cpu420_pval = polyfit(loc_u,Cp_u_four_deg_20,2);
cpl420_pval = polyfit(loc_l,Cp_l_four_deg_20,2);
cn420_pval = cpl420_pval - cpu420_pval;
Cp_u_four_20_fit = polyval(cpu420_pval, chord_spread);
Cp_l_four_20_fit = polyval(cpl420_pval, chord_spread);
%Cn_four_20_fit = calculate_cn_values(cn420_pval,chord,chord_spread);
Cn_four_20_fit = calc_integral(cn420_pval,3.5)/3.5;

Cp_u_four_deg_35 = [four_deg_35(1) four_deg_35_u'] / four_deg_35(end);
Cp_l_four_deg_35 = [four_deg_35(1) four_deg_35_l'] / four_deg_35(end);
cpu435_pval = polyfit(loc_u,Cp_u_four_deg_35,2);
cpl435_pval = polyfit(loc_l,Cp_l_four_deg_35,2);
cn435_pval = cpl435_pval - cpu435_pval;
Cp_u_four_35_fit = polyval(cpu435_pval, chord_spread);
Cp_l_four_35_fit = polyval(cpl435_pval, chord_spread);
%Cn_four_35_fit = calculate_cn_values(cn435_pval,chord,chord_spread);
Cn_four_35_fit = calc_integral(cn435_pval,3.5)/3.5;

Cp_u_four_deg_50 = [four_deg_50(1) four_deg_50_u'] / four_deg_50(end);
Cp_l_four_deg_50 = [four_deg_50(1) four_deg_50_l'] / four_deg_50(end);
cpu450_pval = polyfit(loc_u,Cp_u_four_deg_50,2);
cpl450_pval = polyfit(loc_l,Cp_l_four_deg_50,2);
cn450_pval = cpl450_pval - cpu450_pval;
Cp_u_four_50_fit = polyval(cpu450_pval, chord_spread);
Cp_l_four_50_fit = polyval(cpl450_pval, chord_spread);
%Cn_four_50_fit = calculate_cn_values(cn450_pval,chord,chord_spread);
Cn_four_50_fit = calc_integral(cn450_pval,3.5)/3.5;

% Eight Deg AOA
Cp_u_eight_deg_20 = [eight_deg_20(1) eight_deg_20_u'] / eight_deg_20(end);
Cp_l_eight_deg_20 = [eight_deg_20(1) eight_deg_20_l'] / eight_deg_20(end);
cpu820_pval = polyfit(loc_u,Cp_u_eight_deg_20,2);
cpl820_pval = polyfit(loc_l,Cp_l_eight_deg_20,2);
cn820_pval = cpl820_pval - cpu820_pval;
Cp_u_eight_20_fit = polyval(cpu820_pval, chord_spread);
Cp_l_eight_20_fit = polyval(cpl820_pval, chord_spread);
%Cn_eight_20_fit = calculate_cn_values(cn820_pval,chord,chord_spread);
Cn_eight_20_fit = calc_integral(cn820_pval,3.5)/3.5;

Cp_u_eight_deg_35 = [eight_deg_35(1) eight_deg_35_u'] / eight_deg_35(end);
Cp_l_eight_deg_35 = [eight_deg_35(1) eight_deg_35_l'] / eight_deg_35(end);
cpu835_pval = polyfit(loc_u,Cp_u_eight_deg_35,2);
cpl835_pval = polyfit(loc_l,Cp_l_eight_deg_35,2);
cn835_pval = cpl835_pval - cpu835_pval;
Cp_u_eight_35_fit = polyval(cpu835_pval, chord_spread);
Cp_l_eight_35_fit = polyval(cpl835_pval, chord_spread);
%Cn_eight_35_fit = calculate_cn_values(cn835_pval,chord,chord_spread);
Cn_eight_35_fit = calc_integral(cn835_pval,3.5)/3.5;

Cp_u_eight_deg_50 = [eight_deg_50(1) eight_deg_50_u'] / eight_deg_50(end);
Cp_l_eight_deg_50 = [eight_deg_50(1) eight_deg_50_l'] / eight_deg_50(end);
cpu850_pval = polyfit(loc_u,Cp_u_eight_deg_50,2);
cpl850_pval = polyfit(loc_l,Cp_l_eight_deg_50,2);
cn850_pval = cpl850_pval - cpu850_pval;
Cp_u_eight_50_fit = polyval(cpu850_pval, chord_spread);
Cp_l_eight_50_fit = polyval(cpl850_pval, chord_spread);
%Cn_eight_50_fit = calculate_cn_values(cn850_pval,chord,chord_spread);
Cn_eight_50_fit = calc_integral(cn850_pval,3.5)/3.5;

% Twelve Deg AOA
Cp_u_twelve_deg_20 = [twelve_deg_20(1) twelve_deg_20_u'] / twelve_deg_20(end);
Cp_l_twelve_deg_20 = [twelve_deg_20(1) twelve_deg_20_l'] / twelve_deg_20(end);
cpu1220_pval = polyfit(loc_u,Cp_u_twelve_deg_20,2);
cpl1220_pval = polyfit(loc_l,Cp_l_twelve_deg_20,2);
cn1220_pval = cpl1220_pval - cpu1220_pval;
Cp_u_twelve_20_fit = polyval(cpu1220_pval, chord_spread);
Cp_l_twelve_20_fit = polyval(cpl1220_pval, chord_spread);
%Cn_twelve_20_fit = calculate_cn_values(cn1220_pval,chord,chord_spread);
Cn_twelve_20_fit = calc_integral(cn1220_pval,3.5)/3.5;

Cp_u_twelve_deg_35 = [twelve_deg_35(1) twelve_deg_35_u'] / twelve_deg_35(end);
Cp_l_twelve_deg_35 = [twelve_deg_35(1) twelve_deg_35_l'] / twelve_deg_35(end);
cpu1235_pval = polyfit(loc_u,Cp_u_twelve_deg_35,2);
cpl1235_pval = polyfit(loc_l,Cp_l_twelve_deg_35,2);
cn1235_pval = cpl1235_pval - cpu1235_pval;
Cp_u_twelve_35_fit = polyval(cpu1235_pval, chord_spread);
Cp_l_twelve_35_fit = polyval(cpl1235_pval, chord_spread);
%Cn_twelve_35_fit = calculate_cn_values(cn1235_pval,chord,chord_spread);
Cn_twelve_35_fit = calc_integral(cn1235_pval,3.5)/3.5;

Cp_u_twelve_deg_50 = [twelve_deg_50(1) twelve_deg_50_u'] / twelve_deg_50(end);
Cp_l_twelve_deg_50 = [twelve_deg_50(1) twelve_deg_50_l'] / twelve_deg_50(end);
cpu1250_pval = polyfit(loc_u,Cp_u_twelve_deg_50,2);
cpl1250_pval = polyfit(loc_l,Cp_l_twelve_deg_50,2);
cn1250_pval = cpl1250_pval - cpu1250_pval;
Cp_u_twelve_50_fit = polyval(cpu1250_pval, chord_spread);
Cp_l_twelve_50_fit = polyval(cpl1250_pval, chord_spread);
%Cn_twelve_50_fit = calculate_cn_values(cn1250_pval,chord,chord_spread);
Cn_twelve_50_fit = calc_integral(cn1250_pval,3.5)/3.5;

% Sixteen Deg AOA
Cp_u_sixteen_deg_20 = [sixteen_deg_20(1) sixteen_deg_20_u'] / sixteen_deg_20(end);
Cp_l_sixteen_deg_20 = [sixteen_deg_20(1) sixteen_deg_20_l'] / sixteen_deg_20(end);
cpu1620_pval = polyfit(loc_u,Cp_u_sixteen_deg_20,2);
cpl1620_pval = polyfit(loc_l,Cp_l_sixteen_deg_20,2);
cn1620_pval = cpl1620_pval - cpu1620_pval;
Cp_u_sixteen_20_fit = polyval(cpu1620_pval, chord_spread);
Cp_l_sixteen_20_fit = polyval(cpl1620_pval, chord_spread);
%Cn_sixteen_20_fit = calculate_cn_values(cn1620_pval,chord,chord_spread);
Cn_sixteen_20_fit = calc_integral(cn1620_pval,3.5)/3.5;

Cp_u_sixteen_deg_35 = [sixteen_deg_35(1) sixteen_deg_35_u'] / sixteen_deg_35(end);
Cp_l_sixteen_deg_35 = [sixteen_deg_35(1) sixteen_deg_35_l'] / sixteen_deg_35(end);
cpu1635_pval = polyfit(loc_u,Cp_u_sixteen_deg_35,2);
cpl1635_pval = polyfit(loc_l,Cp_l_sixteen_deg_35,2);
cn1635_pval = cpl1635_pval - cpu1635_pval;
Cp_u_sixteen_35_fit = polyval(cpu1635_pval, chord_spread);
Cp_l_sixteen_35_fit = polyval(cpl1635_pval, chord_spread);
%Cn_sixteen_35_fit = calculate_cn_values(cn1635_pval,chord,chord_spread);
Cn_sixteen_35_fit = calc_integral(cn1635_pval,3.5)/3.5;

Cp_u_sixteen_deg_50 = [sixteen_deg_50(1) sixteen_deg_50_u'] / sixteen_deg_50(end);
Cp_l_sixteen_deg_50 = [sixteen_deg_50(1) sixteen_deg_50_l'] / sixteen_deg_50(end);
cpu1650_pval = polyfit(loc_u,Cp_u_sixteen_deg_50,2);
cpl1650_pval = polyfit(loc_l,Cp_l_sixteen_deg_50,2);
cn1650_pval = cpl1650_pval - cpu1650_pval;
Cp_u_sixteen_50_fit = polyval(cpu1650_pval, chord_spread);
Cp_l_sixteen_50_fit = polyval(cpl1650_pval, chord_spread);
%Cn_sixteen_50_fit = calculate_cn_values(cn1650_pval,chord,chord_spread);
Cn_sixteen_50_fit = calc_integral(cn1650_pval,3.5)/3.5;

%% Calculating Cl and Cd

% Zero Deg AOA
Cl_zero_deg_20 = Cn_zero_20_fit*cosd(0);
Cd_zero_deg_20 = Cn_zero_20_fit*sind(0);
Cl_zero_deg_35 = Cn_zero_35_fit*cosd(0);
Cd_zero_deg_35 = Cn_zero_35_fit*sind(0);
Cl_zero_deg_50 = Cn_zero_50_fit*cosd(0);
Cd_zero_deg_50 = Cn_zero_50_fit*sind(0);

% Four Deg AOA
Cl_four_deg_20 = Cn_four_20_fit*cosd(4);
Cd_four_deg_20 = Cn_four_20_fit*sind(4);
Cl_four_deg_35 = Cn_four_35_fit*cosd(4);
Cd_four_deg_35 = Cn_four_35_fit*sind(4);
Cl_four_deg_50 = Cn_four_50_fit*cosd(4);
Cd_four_deg_50 = Cn_four_50_fit*sind(4);

% Eight Deg AOA
Cl_eight_deg_20 = Cn_eight_20_fit*cosd(8);
Cd_eight_deg_20 = Cn_eight_20_fit*sind(8);
Cl_eight_deg_35 = Cn_eight_35_fit*cosd(8);
Cd_eight_deg_35 = Cn_eight_35_fit*sind(8);
Cl_eight_deg_50 = Cn_eight_50_fit*cosd(8);
Cd_eight_deg_50 = Cn_eight_50_fit*sind(8);

% Twelve Deg AOA
Cl_twelve_deg_20 = Cn_twelve_20_fit*cosd(12);
Cd_twelve_deg_20 = Cn_twelve_20_fit*sind(12);
Cl_twelve_deg_35 = Cn_twelve_35_fit*cosd(12);
Cd_twelve_deg_35 = Cn_twelve_35_fit*sind(12);
Cl_twelve_deg_50 = Cn_twelve_50_fit*cosd(12);
Cd_twelve_deg_50 = Cn_twelve_50_fit*sind(12);

% Sixteen Deg AOA
Cl_sixteen_deg_20 = Cn_sixteen_20_fit*cosd(16);
Cd_sixteen_deg_20 = Cn_sixteen_20_fit*sind(16);
Cl_sixteen_deg_35 = Cn_sixteen_35_fit*cosd(16);
Cd_sixteen_deg_35 = Cn_sixteen_35_fit*sind(16);
Cl_sixteen_deg_50 = Cn_sixteen_50_fit*cosd(16);
Cd_sixteen_deg_50 = Cn_sixteen_50_fit*sind(16);

%% Organizing Values into Cl and Cd Matrices by their Speeds
Cl_20 = [Cl_zero_deg_20 Cl_four_deg_20 Cl_eight_deg_20 Cl_twelve_deg_20 Cl_sixteen_deg_20];
Cl_35 = [Cl_zero_deg_35 Cl_four_deg_35 Cl_eight_deg_35 Cl_twelve_deg_35 Cl_sixteen_deg_35];
Cl_50 = [Cl_zero_deg_50 Cl_four_deg_50 Cl_eight_deg_50 Cl_twelve_deg_50 Cl_sixteen_deg_50];

Cd_20 = [Cd_zero_deg_20 Cd_four_deg_20 Cd_eight_deg_20 Cd_twelve_deg_20 Cd_sixteen_deg_20];
Cd_35 = [Cd_zero_deg_35 Cd_four_deg_35 Cd_eight_deg_35 Cd_twelve_deg_35 Cd_sixteen_deg_35];
Cd_50 = [Cd_zero_deg_50 Cd_four_deg_50 Cd_eight_deg_50 Cd_twelve_deg_50 Cd_sixteen_deg_50];

%% Plotting Cp_u and Cp_l
loc_u_perc = [0 7.5 10 20 30 40 50 60 70 80];
loc_l_perc = [0 7.5 10 20 30 40 50 60 70];

figure(6)
hold on
plot(loc_u_perc, Cp_u_zero_deg_20,'-r')
plot(loc_u_perc, Cp_u_four_deg_20,'-m')
plot(loc_u_perc, Cp_u_eight_deg_20,'-b')
plot(loc_u_perc, Cp_u_twelve_deg_20,'-g')
plot(loc_u_perc, Cp_u_sixteen_deg_20,'-c')
hold off
title('C_{pu} Along Chord at 20 m/s and Various AOA')
xlabel('Percent Along Chord (%)')
ylabel('C_{pu}')
grid('on')
legend('Zero deg', 'Four deg', 'Eight deg', 'Twelve deg', 'Sixteen deg', Location='SouthEast')

figure(7)
hold on
plot(loc_u_perc, Cp_u_zero_deg_35,'-r')
plot(loc_u_perc, Cp_u_four_deg_35,'-m')
plot(loc_u_perc, Cp_u_eight_deg_35,'-b')
plot(loc_u_perc, Cp_u_twelve_deg_35,'-g')
plot(loc_u_perc, Cp_u_sixteen_deg_35,'-c')
hold off
title('C_{pu} Along Chord at 35 m/s and Various AOA')
xlabel('Percent Along Chord (%)')
ylabel('C_{pu}')
grid('on')
legend('Zero deg', 'Four deg', 'Eight deg', 'Twelve deg', 'Sixteen deg', Location='SouthEast')

figure(8)
hold on
plot(loc_u_perc, Cp_u_zero_deg_50,'-r')
plot(loc_u_perc, Cp_u_four_deg_50,'-m')
plot(loc_u_perc, Cp_u_eight_deg_50,'-b')
plot(loc_u_perc, Cp_u_twelve_deg_50,'-g')
plot(loc_u_perc, Cp_u_sixteen_deg_50,'-c')
hold off
title('C_{pu} Along Chord at 50 m/s and Various AOA')
xlabel('Percent Along Chord (%)')
ylabel('C_{pu}')
grid('on')
legend('Zero deg', 'Four deg', 'Eight deg', 'Twelve deg', 'Sixteen deg', Location='SouthEast')

figure(9)
hold on
plot(loc_l_perc, Cp_l_zero_deg_20,'-r')
plot(loc_l_perc, Cp_l_four_deg_20,'-m')
plot(loc_l_perc, Cp_l_eight_deg_20,'-b')
plot(loc_l_perc, Cp_l_twelve_deg_20,'-g')
plot(loc_l_perc, Cp_l_sixteen_deg_20,'-c')
hold off
title('C_{pl} Along Chord at 20 m/s and Various AOA')
xlabel('Percent Along Chord (%)')
ylabel('C_{pl}')
grid('on')
legend('Zero deg', 'Four deg', 'Eight deg', 'Twelve deg', 'Sixteen deg', Location='SouthEast')

figure(10)
hold on
plot(loc_l_perc, Cp_l_zero_deg_35,'-r')
plot(loc_l_perc, Cp_l_four_deg_35,'-m')
plot(loc_l_perc, Cp_l_eight_deg_35,'-b')
plot(loc_l_perc, Cp_l_twelve_deg_35,'-g')
plot(loc_l_perc, Cp_l_sixteen_deg_35,'-c')
hold off
title('C_{pl} Along Chord at 35 m/s and Various AOA')
xlabel('Percent Along Chord (%)')
ylabel('C_{pl}')
grid('on')
legend('Zero deg', 'Four deg', 'Eight deg', 'Twelve deg', 'Sixteen deg', Location='SouthEast')

figure(11)
hold on
plot(loc_l_perc, Cp_l_zero_deg_50,'-r')
plot(loc_l_perc, Cp_l_four_deg_50,'-m')
plot(loc_l_perc, Cp_l_eight_deg_50,'-b')
plot(loc_l_perc, Cp_l_twelve_deg_50,'-g')
plot(loc_l_perc, Cp_l_sixteen_deg_50,'-c')
hold off
title('C_{pl} Along Chord at 50 m/s and Various AOA')
xlabel('Percent Along Chord (%)')
ylabel('C_{pl}')
grid('on')
legend('Zero deg', 'Four deg', 'Eight deg', 'Twelve deg', 'Sixteen deg', Location='SouthEast')

%% Plotting Cl and Cd
angle_spread = [0 1 2 3 4]*4;
figure(16)
hold on
plot(angle_spread,Cl_20,'-r')
plot(angle_spread,Cl_35,'-b')
plot(angle_spread,Cl_50,'-g')
hold off
grid('on') 
legend('20 m/s', '35 m/s', '50 m/s', location='NorthWest')
title('C_l by AOA at Various Air Speeds')
xlabel('Angle of Attack (deg)')
ylabel('C_l')

figure(17)
hold on
plot(angle_spread,Cd_20,'-r')
plot(angle_spread,Cd_35,'-b')
plot(angle_spread,Cd_50,'-g')
hold off
grid('on')
legend('20 m/s', '35 m/s', '50 m/s', location='NorthWest')
title('C_d by AOA at Various Air Speeds')
xlabel('Angle of Attack (deg)')
ylabel('C_d')

%% Calculating Cm_le
Cm_le_zero_deg_20 = calc_cm_integral(cpu020_pval, cpl020_pval, zero_deg_20(end), 3.5);
Cm_le_zero_deg_35 = calc_cm_integral(cpu035_pval, cpl035_pval, zero_deg_35(end), 3.5);
Cm_le_zero_deg_50 = calc_cm_integral(cpu050_pval, cpl050_pval, zero_deg_50(end), 3.5);

Cm_le_four_deg_20 = calc_cm_integral(cpu420_pval, cpl420_pval, four_deg_20(end), 3.5);
Cm_le_four_deg_35 = calc_cm_integral(cpu435_pval, cpl435_pval, four_deg_35(end), 3.5);
Cm_le_four_deg_50 = calc_cm_integral(cpu450_pval, cpl450_pval, four_deg_50(end), 3.5);

Cm_le_eight_deg_20 = calc_cm_integral(cpu820_pval, cpl820_pval, eight_deg_20(end), 3.5);
Cm_le_eight_deg_35 = calc_cm_integral(cpu835_pval, cpl835_pval, eight_deg_35(end), 3.5);
Cm_le_eight_deg_50 = calc_cm_integral(cpu850_pval, cpl850_pval, eight_deg_50(end), 3.5);

Cm_le_twelve_deg_20 = calc_cm_integral(cpu1220_pval, cpl1220_pval, twelve_deg_20(end), 3.5);
Cm_le_twelve_deg_35 = calc_cm_integral(cpu1235_pval, cpl1235_pval, twelve_deg_35(end), 3.5);
Cm_le_twelve_deg_50 = calc_cm_integral(cpu1250_pval, cpl1250_pval, twelve_deg_50(end), 3.5);

Cm_le_sixteen_deg_20 = calc_cm_integral(cpu1620_pval, cpl1620_pval, sixteen_deg_20(end), 3.5);
Cm_le_sixteen_deg_35 = calc_cm_integral(cpu1635_pval, cpl1635_pval, sixteen_deg_35(end), 3.5);
Cm_le_sixteen_deg_50 = calc_cm_integral(cpu1650_pval, cpl1650_pval, sixteen_deg_50(end), 3.5);

%% Organizing Cm_le
Cm_le_20 = [Cm_le_zero_deg_20 Cm_le_four_deg_20 Cm_le_eight_deg_20 Cm_le_twelve_deg_20 Cm_le_sixteen_deg_20];
Cm_le_35 = [Cm_le_zero_deg_35 Cm_le_four_deg_35 Cm_le_eight_deg_35 Cm_le_twelve_deg_35 Cm_le_sixteen_deg_35];
Cm_le_50 = [Cm_le_zero_deg_50 Cm_le_four_deg_50 Cm_le_eight_deg_50 Cm_le_twelve_deg_50 Cm_le_sixteen_deg_50];

%% Plotting Cm_le
figure(18)
hold on
plot(angle_spread,Cm_le_20,'-r')
plot(angle_spread,Cm_le_35,'-b')
plot(angle_spread,Cm_le_50,'-g')
hold off
title('C_m by AOA at Various Air Speeds')
grid('on')
xlabel('Angle of Attack (deg)')
ylabel('C_m')
legend('20 m/s', '35 m/s', '50 m/s', location='SouthWest')


%% Calculating x_cp
Cn_20 = [Cn_zero_20_fit Cn_four_20_fit Cn_eight_20_fit Cn_twelve_20_fit Cn_sixteen_20_fit];
Cn_35 = [Cn_zero_35_fit Cn_four_35_fit Cn_eight_35_fit Cn_twelve_35_fit Cn_sixteen_35_fit];
Cn_50 = [Cn_zero_50_fit Cn_four_50_fit Cn_eight_50_fit Cn_twelve_50_fit Cn_sixteen_50_fit];

x_cp_20 = -Cm_le_20*3.5./Cn_20;
x_cp_35 = -Cm_le_35*3.5./Cn_35;
x_cp_50 = -Cm_le_50*3.5./Cn_50;

figure(19)
hold on
plot(angle_spread,x_cp_20,'-r')
plot(angle_spread,x_cp_35,'-b')
plot(angle_spread,x_cp_50,'-g')
hold off
grid('on')
title('Center of Pressure x_{cp} over AOA at Various Speeds')
xlabel('Angle of Attack (deg)')
ylabel('Position on Chord from Leading Edge (in)')
legend('20 m/s', '35 m/s', '50 m/s', location='East')

%% Polyfitting Dynamic Pressure Distr. of Rake
q020_pval = polyfit(rake_spread,zero_deg_20(19:36),2);
q035_pval = polyfit(rake_spread,zero_deg_35(19:36),2);
q050_pval = polyfit(rake_spread,zero_deg_50(19:36),2);

q420_pval = polyfit(rake_spread,four_deg_20(19:36),2);
q435_pval = polyfit(rake_spread,four_deg_35(19:36),2);
q450_pval = polyfit(rake_spread,four_deg_50(19:36),2);

q820_pval = polyfit(rake_spread,eight_deg_20(19:36),2);
q835_pval = polyfit(rake_spread,eight_deg_35(19:36),2);
q850_pval = polyfit(rake_spread,eight_deg_50(19:36),2);

q1220_pval = polyfit(rake_spread,twelve_deg_20(19:36),2);
q1235_pval = polyfit(rake_spread,twelve_deg_35(19:36),2);
q1250_pval = polyfit(rake_spread,twelve_deg_50(19:36),2);

q1620_pval = polyfit(rake_spread,sixteen_deg_20(19:36),2);
q1635_pval = polyfit(rake_spread,sixteen_deg_35(19:36),2);
q1650_pval = polyfit(rake_spread,sixteen_deg_50(19:36),2);

%% Calculating Cd
R = 1.7/3.5;

Cd_md_zero_deg_20 = R - 1/3.5/zero_deg_20(end) * calc_cd_md_integral(q020_pval);
Cd_md_zero_deg_35 = R - 1/3.5/zero_deg_35(end) * calc_cd_md_integral(q035_pval);
Cd_md_zero_deg_50 = R - 1/3.5/zero_deg_50(end) * calc_cd_md_integral(q050_pval);

Cd_md_four_deg_20 = R - 1/3.5/four_deg_20(end) * calc_cd_md_integral(q420_pval);
Cd_md_four_deg_35 = R - 1/3.5/four_deg_35(end) * calc_cd_md_integral(q435_pval);
Cd_md_four_deg_50 = R - 1/3.5/four_deg_50(end) * calc_cd_md_integral(q450_pval);

Cd_md_eight_deg_20 = R - 1/3.5/eight_deg_20(end) * calc_cd_md_integral(q820_pval);
Cd_md_eight_deg_35 = R - 1/3.5/eight_deg_35(end) * calc_cd_md_integral(q835_pval);
Cd_md_eight_deg_50 = R - 1/3.5/eight_deg_50(end) * calc_cd_md_integral(q850_pval);

Cd_md_twelve_deg_20 = R - 1/3.5/twelve_deg_20(end) * calc_cd_md_integral(q1220_pval);
Cd_md_twelve_deg_35 = R - 1/3.5/twelve_deg_35(end) * calc_cd_md_integral(q1235_pval);
Cd_md_twelve_deg_50 = R - 1/3.5/twelve_deg_50(end) * calc_cd_md_integral(q1250_pval);

Cd_md_sixteen_deg_20 = R - 1/3.5/sixteen_deg_20(end) * calc_cd_md_integral(q1620_pval);
Cd_md_sixteen_deg_35 = R - 1/3.5/sixteen_deg_35(end) * calc_cd_md_integral(q1635_pval);
Cd_md_sixteen_deg_50 = R - 1/3.5/sixteen_deg_50(end) * calc_cd_md_integral(q1650_pval);

%% Organizing Cd
Cd_md_20 = [Cd_md_zero_deg_20 Cd_md_four_deg_20 Cd_md_eight_deg_20 Cd_md_twelve_deg_20 Cd_md_sixteen_deg_20];
Cd_md_35 = [Cd_md_zero_deg_35 Cd_md_four_deg_35 Cd_md_eight_deg_35 Cd_md_twelve_deg_35 Cd_md_sixteen_deg_35];
Cd_md_50 = [Cd_md_zero_deg_50 Cd_md_four_deg_50 Cd_md_eight_deg_50 Cd_md_twelve_deg_50 Cd_md_sixteen_deg_50];

%% Plotting Cd
figure(20)
hold on
plot(angle_spread,Cd_md_20,'-r')
plot(angle_spread,Cd_md_35,'-b')
plot(angle_spread,Cd_md_50,'-g')
hold off
title('Momentum Deficit C_d by AOA at Various Speeds')
grid('on')
xlabel('Angle of Attack (deg)')
ylabel('C_d')
legend('20 m/s', '35 m/s', '50 m/s', location='NorthWest')

%% Plotting Dynamic Pressure Ratio Over Rake Distr.
dyn_pre_r_zero_deg_20 = zero_deg_20(19:36)/zero_deg_20(end);
dyn_pre_r_zero_deg_35 = zero_deg_35(19:36)/zero_deg_35(end);
dyn_pre_r_zero_deg_50 = zero_deg_50(19:36)/zero_deg_50(end);

dyn_pre_r_four_deg_20 = four_deg_20(19:36)/four_deg_20(end);
dyn_pre_r_four_deg_35 = four_deg_35(19:36)/four_deg_35(end);
dyn_pre_r_four_deg_50 = four_deg_50(19:36)/four_deg_50(end);

dyn_pre_r_eight_deg_20 = eight_deg_20(19:36)/eight_deg_20(end);
dyn_pre_r_eight_deg_35 = eight_deg_35(19:36)/eight_deg_35(end);
dyn_pre_r_eight_deg_50 = eight_deg_50(19:36)/eight_deg_50(end);

dyn_pre_r_twelve_deg_20 = twelve_deg_20(19:36)/twelve_deg_20(end);
dyn_pre_r_twelve_deg_35 = twelve_deg_35(19:36)/twelve_deg_35(end);
dyn_pre_r_twelve_deg_50 = twelve_deg_50(19:36)/twelve_deg_50(end);

dyn_pre_r_sixteen_deg_20 = sixteen_deg_20(19:36)/sixteen_deg_20(end);
dyn_pre_r_sixteen_deg_35 = sixteen_deg_35(19:36)/sixteen_deg_35(end);
dyn_pre_r_sixteen_deg_50 = sixteen_deg_50(19:36)/sixteen_deg_50(end);


figure(21)
hold on
plot(rake_spread,dyn_pre_r_zero_deg_20,'-r')
plot(rake_spread,dyn_pre_r_four_deg_20,'-m')
plot(rake_spread,dyn_pre_r_eight_deg_20,'-b')
plot(rake_spread,dyn_pre_r_twelve_deg_20,'g')
plot(rake_spread,dyn_pre_r_sixteen_deg_20,'-c')
hold off
grid('on')
title('Dynamic Pressure Ratio over Rake Distr. at 20 m/s')
legend('Zero deg', 'Four deg', 'Eight deg', 'Twelve deg', 'Sixteen deg', Location='SouthWest')

figure(22)
hold on
plot(rake_spread,dyn_pre_r_zero_deg_35,'-r')
plot(rake_spread,dyn_pre_r_four_deg_35,'-m')
plot(rake_spread,dyn_pre_r_eight_deg_35,'-b')
plot(rake_spread,dyn_pre_r_twelve_deg_35,'g')
plot(rake_spread,dyn_pre_r_sixteen_deg_35,'-c')
hold off
grid('on')
title('Dynamic Pressure Ratio over Rake Distr. at 35 m/s')
legend('Zero deg', 'Four deg', 'Eight deg', 'Twelve deg', 'Sixteen deg', Location='SouthWest')

figure(23)
hold on
plot(rake_spread,dyn_pre_r_zero_deg_50,'-r')
plot(rake_spread,dyn_pre_r_four_deg_50,'-m')
plot(rake_spread,dyn_pre_r_eight_deg_50,'-b')
plot(rake_spread,dyn_pre_r_twelve_deg_50,'g')
plot(rake_spread,dyn_pre_r_sixteen_deg_50,'-c')
hold off
grid('on')
title('Dynamic Pressure Ratio over Rake Distr. at 50 m/s')
legend('Zero deg', 'Four deg', 'Eight deg', 'Twelve deg', 'Sixteen deg', Location='SouthWest')

%% Organizing Values into Cl and Cd Matrices by their AOA
Cl_bank = [Cl_20; Cl_35; Cl_50];
Cl_zero = Cl_bank(:,1); Cl_four = Cl_bank(:,2); Cl_eight = Cl_bank(:,3); Cl_twelve = Cl_bank(:,4); Cl_sixteen = Cl_bank(:,5);

Cd_bank = [Cd_20; Cd_35; Cd_50];
Cd_zero = Cd_bank(:,1); Cd_four = Cd_bank(:,2); Cd_eight = Cd_bank(:,3); Cd_twelve = Cd_bank(:,4); Cd_sixteen = Cd_bank(:,5);

%% Plotting Cd and Cl across Reynold's Number
ReyNum = [20 35 50]*rho*3.5/mu;

figure(24)
hold on
plot(ReyNum,Cl_zero,'-r')
plot(ReyNum,Cl_four,'-m')
plot(ReyNum,Cl_eight,'-b')
plot(ReyNum,Cl_twelve,'-g')
plot(ReyNum,Cl_sixteen,'-c')
hold off
title('C_l over Reynold''s Number at Various AOA')
grid('on')
xlabel('Reynold''s Number')
ylabel('C_l')
legend('Zero deg', 'Four deg', 'Eight deg', 'Twelve deg', 'Sixteen deg')

figure(25)
hold on
plot(ReyNum,Cd_zero,'-r')
plot(ReyNum,Cd_four,'-m')
plot(ReyNum,Cd_eight,'-b')
plot(ReyNum,Cd_twelve,'-g')
plot(ReyNum,Cd_sixteen,'-c')
hold off
title('C_d over Reynold''s Number at Various AOA')
grid('on')
xlabel('Reynold''s Number')
ylabel('C_d')
legend('Zero deg', 'Four deg', 'Eight deg', 'Twelve deg', 'Sixteen deg')

%% Defined Functions
function output_value = calc_integral(cn,c)
inte = @(x) cn(1)*x.^2 + cn(2)*x + cn(3);
output_value = integral(inte,0,c);
end


function output_values = calculate_cn_values(cn,c,input_values)
calculate_values = @(x) (cn(1)/3*x.^3 + cn(2)/2*x.^2 + cn(3)*x)./c;
output_values = calculate_values(input_values);
end

function output_value = calc_cm_integral(cpu, cpl, dyn, c)
cpul = cpu - cpl;
cm_inte = @(x) (cpul(1)*x.^2 + cpul(2)*x + cpul(3))/dyn.*(x);
output_value = integral(cm_inte,0,c)/c^2;
end

function output = calc_cd_md_integral(pval)
dum_func = @(x) pval(1)*x.^2 + pval(2)*x + pval(3);
output = integral(dum_func,-1.7/2,1.7/2);
end