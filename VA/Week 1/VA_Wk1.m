clear all; close all; clc;

%% Floor 1 & Floor 2

%Read Experimental Files
    %Mode 1 - 8.9 Hz
    F12_M1 = csvread('F12_M1_89.csv',11);
        F12_M1_Time = F12_M1(:,1);
        F12__M1_Ch1 = F12_M1(:,2);
        F12__M1_Ch2 = F12_M1(:,3);
    %Mode 2 - 13.5 Hz
    F12_M2 = csvread('F12_M2_135.csv',11);
        F12_M2_Time = F12_M2(:,1);
        F12__M2_Ch1 = F12_M2(:,2);
        F12__M2_Ch2 = F12_M2(:,3);
    %Mode 3 - 24.5 Hz
    F12_M3 = csvread('F12_M3_245.csv',11);
        F12_M3_Time = F12_M3(:,1);
        F12__M3_Ch1 = F12_M3(:,2);
        F12__M3_Ch2 = F12_M3(:,3);

%Data Smoothing - Mode 1
smooth1 = smoothdata(F12__M1_Ch1,"movmean",50);
smooth2 = smoothdata(F12__M1_Ch2,"movmedian",40);

smooth1 = smoothdata(smooth1,"sgolay");
smooth2 = smoothdata(smooth2,"sgolay");

smooth1 = smoothdata(smooth1,"gaussian");
        
%Plotting - Mode 1
    %Smoothed
        figure()
        hold on
        plot(F12_M1_Time,smooth1)
        plot(F12_M1_Time,smooth2)
        grid on
        legend('Floor 1','Floor 2','Location','southeast')
        title('Mode 1 - Floor 1 & Floor 2 [Smoothed Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])
    %Raw
        figure()
        hold on
        plot(F12_M1_Time,(F12__M1_Ch1))
        plot(F12_M1_Time,(F12__M1_Ch2))
        grid on
        legend('Floor 1','Floor 2','Location','southeast')
        title('Mode 1 - Floor 1 & Floor 2 [Raw Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])

%Data Smoothing - Mode 2
smooth1 = smoothdata(F12__M2_Ch1,"sgolay");
smooth2 = smoothdata(F12__M2_Ch2,"sgolay");

%Data Saving
Floor_2_M2 = smooth2; %Save Mode 1 Data

%Plotting - Mode 2
    %Smoothed
        figure()
        hold on
        plot(F12_M2_Time,smooth1)
        plot(F12_M2_Time,smooth2)
        grid on
        legend('Floor 1','Floor 2','Location','southeast')
        title('Mode 2 - Floor 1 & Floor 2 [Smoothed Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])
    %Raw
        figure()
        hold on
        plot(F12_M2_Time,(F12__M2_Ch1))
        plot(F12_M2_Time,(F12__M2_Ch2))
        grid on
        legend('Floor 1','Floor 2','Location','southeast')
        title('Mode 2 - Floor 1 & Floor 2 [Raw Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])


%Data Smoothing - Mode 3
smooth1 = smoothdata(F12__M3_Ch1,"sgolay");
smooth2 = smoothdata(F12__M3_Ch2,"sgolay");

%Plotting - Mode 3
    %Smoothed
        figure()
        hold on
        plot(F12_M3_Time,smooth1)
        plot(F12_M3_Time,smooth2)
        grid on
        legend('Floor 1','Floor 2','Location','southeast')
        title('Mode 3 - Floor 1 & Floor 2 [Smoothed Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])
    %Raw
        figure()
        hold on
        plot(F12_M3_Time,(F12__M3_Ch1))
        plot(F12_M3_Time,(F12__M3_Ch2))
        grid on
        legend('Floor 1','Floor 2','Location','southeast')
        title('Mode 3 - Floor 1 & Floor 2 [Raw Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])        

%% Floor 1 & Floor 3

%Read Experimental Files
    %Mode 1 - 8.9 Hz
    F13_M1 = csvread('F13_M1_89_V1.csv',11);
        F13_M1_Time = F13_M1(:,1);
        F13__M1_Ch1 = F13_M1(:,2);
        F13__M1_Ch2 = F13_M1(:,3);
    %Mode 2 - 12.7 Hz
    F13_M2 = csvread('F13_M2_127.csv',11);
        F13_M2_Time = F13_M2(:,1);
        F13__M2_Ch1 = F13_M2(:,2);
        F13__M2_Ch2 = F13_M2(:,3);
    %Mode 3 - 24.7 Hz
    F13_M3 = csvread('F13_M3_247_V1.csv',11);
        F13_M3_Time = F13_M3(:,1);
        F13__M3_Ch1 = F13_M3(:,2);
        F13__M3_Ch2 = F13_M3(:,3);

%Data Smoothing - Mode 1
smooth1 = smoothdata(F13__M1_Ch1,"movmean",40);

smooth1 = smoothdata(smooth1,"sgolay");
smooth2 = smoothdata(F13__M1_Ch2,"sgolay");

smooth1 = smoothdata(smooth1,"gaussian");

%Data Saving
Floor_1_M1 = smooth1; %Save Mode 1 Data
        
%Plotting - Mode 1
    %Smoothed
        figure()
        hold on
        plot(F13_M1_Time,smooth1)
        plot(F13_M1_Time,smooth2)
        grid on
        legend('Floor 1','Floor 3','Location','southeast')
        title('Mode 1 - Floor 1 & Floor 3 [Smoothed Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])
    %Raw
        figure()
        hold on
        plot(F13_M1_Time,(F13__M1_Ch1))
        plot(F13_M1_Time,(F13__M1_Ch2))
        grid on
        legend('Floor 1','Floor 3','Location','southeast')
        title('Mode 1 - Floor 1 & Floor 3 [Raw Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])

%Data Smoothing - Mode 2
smooth1 = smoothdata(F13__M2_Ch1,"sgolay");
smooth2 = smoothdata(F13__M2_Ch2,"sgolay");

%Plotting - Mode 2
    %Smoothed
        figure()
        hold on
        plot(F13_M2_Time,smooth1)
        plot(F13_M2_Time,smooth2)
        grid on
        legend('Floor 1','Floor 3','Location','southeast')
        title('Mode 2 - Floor 1 & Floor 3 [Smoothed Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])
    %Raw
        figure()
        hold on
        plot(F13_M2_Time,(F13__M2_Ch1))
        plot(F13_M2_Time,(F13__M2_Ch2))
        grid on
        legend('Floor 1','Floor 3','Location','southeast')
        title('Mode 2 - Floor 1 & Floor 3 [Raw Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])

%Data Smoothing - Mode 3
smooth1 = smoothdata(F13__M3_Ch1,"movmean",40);

smooth1 = smoothdata(smooth1,"sgolay");
smooth2 = smoothdata(F13__M3_Ch2,"sgolay");

%Plotting - Mode 3
    %Smoothed
        figure()
        hold on
        plot(F13_M3_Time,smooth1)
        plot(F13_M3_Time,smooth2)
        grid on
        legend('Floor 1','Floor 3','Location','southeast')
        title('Mode 3 - Floor 1 & Floor 3 [Smoothed Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])
    %Raw
        figure()
        hold on
        plot(F13_M3_Time,(F13__M3_Ch1))
        plot(F13_M3_Time,(F13__M3_Ch2))
        grid on
        legend('Floor 1','Floor 3','Location','southeast')
        title('Mode 3 - Floor 1 & Floor 3 [Raw Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])

%% Floor 2 & Floor 3

%Read Experimental Files
    %Mode 1 - 9.0 Hz
    F23_M1 = csvread('F23_M1_9.csv',11);
        F23_M1_Time = F23_M1(:,1);
        F23__M1_Ch1 = F23_M1(:,2);
        F23__M1_Ch2 = F23_M1(:,3);
    %Mode 2 - 13.2 Hz
    F23_M2 = csvread('F23_M2_132.csv',11);
        F23_M2_Time = F23_M2(:,1);
        F23__M2_Ch1 = F23_M2(:,2);
        F23__M2_Ch2 = F23_M2(:,3);
    %Mode 3 - 25.1 Hz
    F23_M3 = csvread('F23_M3_251.csv',11);
        F23_M3_Time = F23_M3(:,1);
        F23__M3_Ch1 = F23_M3(:,2);
        F23__M3_Ch2 = F23_M3(:,3);

%Data Smoothing
smooth1 = smoothdata(F23__M1_Ch1,"movmean",40);

smooth1 = smoothdata(smooth1,"sgolay");
smooth2 = smoothdata(F23__M1_Ch2,"sgolay");

%Data Saving
Floor_2_M1 = smooth1; %Save Mode 1 Data
Floor_3_M1 = smooth2; %Save Mode 1 Data
        
%Plotting - Mode 1
    %Smoothed
        figure()
        hold on
        plot(F23_M1_Time,smooth1)
        plot(F23_M1_Time,smooth2)
        grid on
        legend('Floor 2','Floor 3','Location','southeast')
        title('Mode 1 - Floor 2 & Floor 3 [Smoothed Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])
    %Raw
        figure()
        hold on
        plot(F23_M1_Time,(F23__M1_Ch1))
        plot(F23_M1_Time,(F23__M1_Ch2))
        grid on
        legend('Floor 2','Floor 3','Location','southeast')
        title('Mode 1 - Floor 2 & Floor 3 [Raw Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])

%Data Smoothing - Mode 2
smooth1 = smoothdata(F23__M2_Ch1,"sgolay");
smooth2 = smoothdata(F23__M2_Ch2,"sgolay");

%Plotting - Mode 2
    %Smoothed
        figure()
        hold on
        plot(F23_M2_Time,smooth1)
        plot(F23_M2_Time,smooth2)
        grid on
        legend('Floor 2','Floor 3','Location','southeast')
        title('Mode 2 - Floor 2 & Floor 3 [Smoothed Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])
    %Raw
        figure()
        hold on
        plot(F23_M2_Time,(F23__M2_Ch1))
        plot(F23_M2_Time,(F23__M2_Ch2))
        grid on
        legend('Floor 2','Floor 3','Location','southeast')
        title('Mode 2 - Floor 2 & Floor 3 [Raw Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])

%Data Smoothing - Mode 3
smooth1 = smoothdata(F23__M3_Ch1,"sgolay");
smooth2 = smoothdata(F23__M3_Ch2,"sgolay");

%Plotting - Mode 3
    %Smoothed
        figure()
        hold on
        plot(F23_M3_Time,smooth1)
        plot(F23_M3_Time,smooth2)
        grid on
        legend('Floor 2','Floor 3','Location','southeast')
        title('Mode 3 - Floor 2 & Floor 3 [Smoothed Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])
    %Raw
        figure()
        hold on
        plot(F23_M3_Time,(F23__M3_Ch1))
        plot(F23_M3_Time,(F23__M3_Ch2))
        grid on
        legend('Floor 2','Floor 3','Location','southeast')
        title('Mode 3 - Floor 2 & Floor 3 [Raw Data]')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])


%% Combined Data
    %Mode 1
        figure()
        hold on
        plot(F13_M1_Time,Floor_1_M1)
        plot(F23_M1_Time,Floor_2_M1)
        plot(F23_M1_Time,Floor_3_M1)
        grid on
        legend('Floor 1','Floor 2','Floor 3','Location','southeast')
        title('Resonance Mode 1')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])