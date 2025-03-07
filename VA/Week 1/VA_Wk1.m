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

%Data Saving
Test_M3_F12 = smooth1;

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
    %Mode 1 - Graph
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

 %% Lilys graphs
%Data Smoothing - Mode 2
smooth1 = smoothdata(F13__M2_Ch1,"sgolay");

%Data Saving
Floor_1_M2 = smooth1; %Save Mode 2 Data

%Data Smoothing - Mode 3
smooth1 = smoothdata(F13__M3_Ch1,"movmean",40);

smooth1 = smoothdata(smooth1,"sgolay");

%Data Saving
Floor_1_M3 = smooth1;

%Data Smoothing - Mode 2
smooth1 = smoothdata(F23__M2_Ch1,"sgolay");
smooth2 = smoothdata(F23__M2_Ch2,"sgolay");

%Data Saving
Floor_2_M2 = smooth1;
Floor_3_M2 = smooth2;

%Data Smoothing - Mode 3
smooth1 = smoothdata(F23__M3_Ch1,"sgolay");
smooth2 = smoothdata(F23__M3_Ch2,"sgolay");

%Data Saving
Floor_2_M3 = smooth1;
Floor_3_M3 = smooth2;

    %Mode 1 floors 1&2
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

    %Mode 2 floors 1&2
        figure()
        hold on
        plot(F13_M2_Time,Floor_1_M2)
        plot(F23_M2_Time,Floor_2_M2)
        plot(F23_M2_Time,Floor_3_M2)
        grid on
        legend('Floor 1','Floor 2','Floor 3','Location','southeast')
        title('Resonance Mode 2')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])
    
    %Mode 3 floors 1&2 
        figure()
        hold on
        plot(F13_M3_Time,Floor_1_M3)
        plot(F23_M3_Time,Floor_2_M3)
        plot(F23_M3_Time,Floor_3_M3)
        grid on
        legend('Floor 1','Floor 2','Floor 3','Location','southeast')
        title('Resonance Mode 3')
        xlabel('Time [s]')
        ylabel('Measured Voltage [V]')
        xlim([-0.1,0.1])



%% Calculated Parameters
    %Mode 1 - Amplitude
        Amp_M1_F1 = (max(Floor_1_M1)-min(Floor_1_M1))/2;
        Amp_M1_F2 = (max(Floor_2_M1)-min(Floor_2_M1))/2;
        Amp_M1_F3 = (max(Floor_3_M1)-min(Floor_3_M1))/2;
    %Mode 1 - Phase Shift
       
    %Mode 2 - Amplitude
        Amp_M2_F1 = (max(Floor_1_M2)-min(Floor_1_M2))/2;
        Amp_M2_F2 = (max(Floor_2_M2)-min(Floor_2_M2))/2;
        Amp_M2_F3 = (max(Floor_3_M2)-min(Floor_3_M2))/2;
    %Mode 2 - Phase Shift
    %Mode 3 - Amplitude
        Amp_M3_F1 = (max(Floor_1_M3)-min(Floor_1_M3))/2;
        Amp_M3_F2 = (max(Floor_2_M3)-min(Floor_2_M3))/2;
        Amp_M3_F3 = (max(Floor_3_M3)-min(Floor_3_M3))/2;
    %Mode 3 - Phase Shift


%% Mode 1 Phase Shift

%Mode 1 - Floor 1 Phase Shift
%Data Smoothing - Mode 1
smooth1 = smoothdata(F12__M1_Ch1,"movmean",175);
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

    PS_M1_F1 = (-0.040075+0.016575)/2;


%Mode 1 - Floor 2 Phase Shift
    LMax = islocalmax(Floor_2_M1); %Locate local maxima
    Peaks = Floor_2_M1(LMax); %Find y-axis value of maxima
    Peak_Time = F23_M1_Time(LMax); %Find x-axis value of maxima

    LMin = islocalmin(Floor_2_M1); %Locate local maxima
    Troughs = Floor_2_M1(LMin); %Find y-axis value of maxima
    Trough_Time = F23_M1_Time(LMin); %Find x-axis value of maxima
    
    %Find Minima and Maxima and Filter Noise
    j = 1;
    for i = 1:numel(Peaks)-1
        if Peaks(i) < 0
            DeleteP(j) = i;
            j = j+1;
        else
            if abs(abs(Peak_Time(i+1))-abs(Peak_Time(i))) < 0.005
                if abs(Peaks(i+1))-abs(Peaks(i)) < 0
                DeleteP(j) = i+1;
                j = j+1;
                else
                DeleteP(j) = i;
                j = j+1;
                end
            end
        end
    end

    j = 1;
    for i = 1:numel(Troughs)-1
        if Troughs(i) > 0
            DeleteT(j) = i;
            j = j+1;
        else
            if abs(abs(Trough_Time(i+1))-abs(Trough_Time(i))) < 0.005
                if abs(Troughs(i+1))-abs(Troughs(i)) > 0
                DeleteT(j) = i+1;
                j = j+1;
                else
                DeleteT(j) = i;
                j = j+1;
                end
            end
        end
    end

    Peaks(DeleteP) = [];
    Peak_Time(DeleteP) = [];

    Troughs(DeleteT) = [];
    Trough_Time(DeleteT) = [];

    clear DeleteT
    clear DeleteP

    %Find Peak and Trough Closest to Zero
    Closest_Peak = Peak_Time(abs(Peak_Time) == min(abs(Peak_Time)));
    Closest_Trough = Trough_Time(abs(Trough_Time) == min(abs(Trough_Time)));
    PS_M1_F2 = (Closest_Peak+Closest_Trough)/2;

%Mode 1 - Floor 3 Phase Shift
    LMax = islocalmax(Floor_3_M1); %Locate local maxima
    Peaks = Floor_3_M1(LMax); %Find y-axis value of maxima
    Peak_Time = F23_M1_Time(LMax); %Find x-axis value of maxima

    LMin = islocalmin(Floor_3_M1); %Locate local maxima
    Troughs = Floor_3_M1(LMin); %Find y-axis value of maxima
    Trough_Time = F23_M1_Time(LMin); %Find x-axis value of maxima
    
    %Find Minima and Maxima and Filter Noise
    for k = 1:5
    
    j = 1;
    for i = 1:numel(Peaks)-1
        if Peaks(i) < 0
            DeleteP(j) = i;
            j = j+1;
        elseif Peaks(i+1) > 0
            DeleteP(j) = i+1;
            j = j+1;
        else
            if abs(abs(Peak_Time(i+1))-abs(Peak_Time(i))) < 0.02
                if abs(Peaks(i+1))-abs(Peaks(i)) < 0
                DeleteP(j) = i+1;
                j = j+1;
                else
                DeleteP(j) = i;
                j = j+1;
                end
            end
        end
    end

    j = 1;
    for i = 1:numel(Troughs)-1
        if Troughs(i) > 0
            DeleteT(j) = i;
            j = j+1;
        elseif Troughs(i+1) > 0
            DeleteT(j) = i+1;
            j = j+1;
        else
            if abs(abs(Trough_Time(i+1))-abs(Trough_Time(i))) < 0.05
                if abs(Troughs(i+1))-abs(Troughs(i)) < 0
                DeleteT(j) = i+1;
                j = j+1;
                else
                DeleteT(j) = i;
                j = j+1;
                end
            end
        end
    end

    if DeleteP ~= 0
        Peaks(DeleteP) = [];
        Peak_Time(DeleteP) = [];
    end

    if DeleteT ~= 0
    Troughs(DeleteT) = [];
    Trough_Time(DeleteT) = [];
    end

    DeleteT = 0;
    DeleteP = 0;
    end

   

    %Find Peak and Trough Closest to Zero
    Closest_Peak = Peak_Time(abs(Peak_Time) == min(abs(Peak_Time)));
    Closest_Trough = Trough_Time(abs(Trough_Time) == min(abs(Trough_Time)));
    PS_M1_F3 = (Closest_Peak+Closest_Trough)/2;

%% Mode 2 Phase Shift

%Mode 2 - Floor 1 Phase Shift
    LMax = islocalmax(Floor_1_M2); %Locate local maxima
    Peaks = Floor_1_M2(LMax); %Find y-axis value of maxima
    Peak_Time = F13_M2_Time(LMax); %Find x-axis value of maxima

    LMin = islocalmin(Floor_1_M2); %Locate local maxima
    Troughs = Floor_1_M2(LMin); %Find y-axis value of maxima
    Trough_Time = F13_M2_Time(LMin); %Find x-axis value of maxima

    %Find Minima and Maxima and Filter Noise
    for k = 1:5

    j = 1;
    for i = 1:numel(Peaks)-1
        if Peaks(i) < 0
            DeleteP(j) = i;
            j = j+1;
        else
            if abs(abs(Peak_Time(i+1))-abs(Peak_Time(i))) < 0.01
                if abs(Peaks(i+1))-abs(Peaks(i)) < 0
                DeleteP(j) = i+1;
                j = j+1;
                else
                DeleteP(j) = i;
                j = j+1;
                end
            end
        end
    end

    j = 1;
    for i = 1:numel(Troughs)-1
        if Troughs(i) > 0
            DeleteT(j) = i;
            j = j+1;
        elseif Troughs(i+1) > 0
            DeleteT(j) = i+1;
            j = j+1;
        else
            if abs(abs(Trough_Time(i+1))-abs(Trough_Time(i))) < 0.01
                if abs(Troughs(i+1))-abs(Troughs(i)) < 0
                DeleteT(j) = i+1;
                j = j+1;
                else
                DeleteT(j) = i;
                j = j+1;
                end
            end
        end
    end

    if DeleteP ~= 0
        Peaks(DeleteP) = [];
        Peak_Time(DeleteP) = [];
    end

    if DeleteT ~= 0
    Troughs(DeleteT) = [];
    Trough_Time(DeleteT) = [];
    end

    DeleteT = 0;
    DeleteP = 0;
    end



    %Find Peak and Trough Closest to Zero
    Closest_Peak = Peak_Time(abs(Peak_Time) == min(abs(Peak_Time)));
    Closest_Trough = Trough_Time(abs(Trough_Time) == min(abs(Trough_Time)));
    PS_M2_F1 = (Closest_Peak+Closest_Trough)/2;

%Mode 2 - Floor 2 Phase Shift
    LMax = islocalmax(Floor_2_M2); %Locate local maxima
    Peaks = Floor_2_M2(LMax); %Find y-axis value of maxima
    Peak_Time = F23_M2_Time(LMax); %Find x-axis value of maxima

    LMin = islocalmin(Floor_2_M2); %Locate local maxima
    Troughs = Floor_2_M2(LMin); %Find y-axis value of maxima
    Trough_Time = F23_M2_Time(LMin); %Find x-axis value of maxima

    %Find Minima and Maxima and Filter Noise
    for k = 1:5

    j = 1;
    for i = 1:numel(Peaks)-1
        if Peaks(i) < 0
            DeleteP(j) = i;
            j = j+1;
        else
            if abs(abs(Peak_Time(i+1))-abs(Peak_Time(i))) < 0.01
                if abs(Peaks(i+1))-abs(Peaks(i)) < 0
                DeleteP(j) = i+1;
                j = j+1;
                else
                DeleteP(j) = i;
                j = j+1;
                end
            end
        end
    end

    j = 1;
    for i = 1:numel(Troughs)-1
        if Troughs(i) > 0
            DeleteT(j) = i;
            j = j+1;
        elseif Troughs(i+1) > 0
            DeleteT(j) = i+1;
            j = j+1;
        else
            if abs(abs(Trough_Time(i+1))-abs(Trough_Time(i))) < 0.01
                if abs(Troughs(i+1))-abs(Troughs(i)) < 0
                DeleteT(j) = i+1;
                j = j+1;
                else
                DeleteT(j) = i;
                j = j+1;
                end
            end
        end
    end

    if DeleteP ~= 0
        Peaks(DeleteP) = [];
        Peak_Time(DeleteP) = [];
    end

    if DeleteT ~= 0
    Troughs(DeleteT) = [];
    Trough_Time(DeleteT) = [];
    end

    DeleteT = 0;
    DeleteP = 0;
    end



    %Find Peak and Trough Closest to Zero
    Closest_Peak = Peak_Time(abs(Peak_Time) == min(abs(Peak_Time)));
    Closest_Trough = Trough_Time(abs(Trough_Time) == min(abs(Trough_Time)));
    PS_M2_F2 = (Closest_Peak+Closest_Trough)/2;

%Mode 2 - Floor 2 Phase Shift
    LMax = islocalmax(Floor_3_M2); %Locate local maxima
    Peaks = Floor_3_M2(LMax); %Find y-axis value of maxima
    Peak_Time = F23_M2_Time(LMax); %Find x-axis value of maxima

    LMin = islocalmin(Floor_3_M2); %Locate local maxima
    Troughs = Floor_3_M2(LMin); %Find y-axis value of maxima
    Trough_Time = F23_M2_Time(LMin); %Find x-axis value of maxima

    %Find Minima and Maxima and Filter Noise
    for k = 1:5

    j = 1;
    for i = 1:numel(Peaks)-1
        if Peaks(i) < 0
            DeleteP(j) = i;
            j = j+1;
        else
            if abs(abs(Peak_Time(i+1))-abs(Peak_Time(i))) < 0.01
                if abs(Peaks(i+1))-abs(Peaks(i)) < 0
                DeleteP(j) = i+1;
                j = j+1;
                else
                DeleteP(j) = i;
                j = j+1;
                end
            end
        end
    end

    j = 1;
    for i = 1:numel(Troughs)-1
        if Troughs(i) > 0
            DeleteT(j) = i;
            j = j+1;
        elseif Troughs(i+1) > 0
            DeleteT(j) = i+1;
            j = j+1;
        else
            if abs(abs(Trough_Time(i+1))-abs(Trough_Time(i))) < 0.01
                if abs(Troughs(i+1))-abs(Troughs(i)) < 0
                DeleteT(j) = i+1;
                j = j+1;
                else
                DeleteT(j) = i;
                j = j+1;
                end
            end
        end
    end

    if DeleteP ~= 0
        Peaks(DeleteP) = [];
        Peak_Time(DeleteP) = [];
    end

    if DeleteT ~= 0
    Troughs(DeleteT) = [];
    Trough_Time(DeleteT) = [];
    end

    DeleteT = 0;
    DeleteP = 0;
    end



    %Find Peak and Trough Closest to Zero
    Closest_Peak = Peak_Time(abs(Peak_Time) == min(abs(Peak_Time)));
    Closest_Trough = Trough_Time(abs(Trough_Time) == min(abs(Trough_Time)));
    PS_M2_F3 = (Closest_Peak+Closest_Trough)/2;

%% Mode 3 Phase Shift

%Mode 3 - Floor 1 Phase Shift
    LMax = islocalmax(Test_M3_F12); %Locate local maxima
    Peaks = Test_M3_F12(LMax); %Find y-axis value of maxima
    Peak_Time = F12_M3_Time(LMax); %Find x-axis value of maxima

    LMin = islocalmin(Test_M3_F12); %Locate local maxima
    Troughs = Test_M3_F12(LMin); %Find y-axis value of maxima
    Trough_Time = F12_M3_Time(LMin); %Find x-axis value of maxima

    %Find Minima and Maxima and Filter Noise
    for k = 1:5

    j = 1;
    for i = 1:numel(Peaks)-1
        if Peaks(i) < 0
            DeleteP(j) = i;
            j = j+1;
        else
            if abs(abs(Peak_Time(i+1))-abs(Peak_Time(i))) < 0.01
                if abs(Peaks(i+1))-abs(Peaks(i)) < 0
                DeleteP(j) = i+1;
                j = j+1;
                else
                DeleteP(j) = i;
                j = j+1;
                end
            end
        end
    end

    j = 1;
    for i = 1:numel(Troughs)-1
        if Troughs(i) > 0
            DeleteT(j) = i;
            j = j+1;
        elseif Troughs(i+1) > 0
            DeleteT(j) = i+1;
            j = j+1;
        else
            if abs(abs(Trough_Time(i+1))-abs(Trough_Time(i))) < 0.01
                if abs(Troughs(i+1))-abs(Troughs(i)) < 0
                DeleteT(j) = i+1;
                j = j+1;
                else
                DeleteT(j) = i;
                j = j+1;
                end
            end
        end
    end

    if DeleteP ~= 0
        Peaks(DeleteP) = [];
        Peak_Time(DeleteP) = [];
    end

    if DeleteT ~= 0
    Troughs(DeleteT) = [];
    Trough_Time(DeleteT) = [];
    end

    DeleteT = 0;
    DeleteP = 0;
    end



    %Find Peak and Trough Closest to Zero
    Closest_Peak = Peak_Time(abs(Peak_Time) == min(abs(Peak_Time)));
    Closest_Trough = Trough_Time(4);
    PS_M3_F1 = (Closest_Peak+Closest_Trough)/2;

%Mode 3 - Floor 2 Phase Shift
    LMax = islocalmax(Floor_2_M3); %Locate local maxima
    Peaks = Floor_2_M3(LMax); %Find y-axis value of maxima
    Peak_Time = F23_M3_Time(LMax); %Find x-axis value of maxima

    LMin = islocalmin(Floor_2_M3); %Locate local maxima
    Troughs = Floor_2_M3(LMin); %Find y-axis value of maxima
    Trough_Time = F23_M3_Time(LMin); %Find x-axis value of maxima

    %Find Minima and Maxima and Filter Noise
    for k = 1:5

    j = 1;
    for i = 1:numel(Peaks)-1
        if Peaks(i) < 0
            DeleteP(j) = i;
            j = j+1;
        else
            if abs(abs(Peak_Time(i+1))-abs(Peak_Time(i))) < 0.01
                if abs(Peaks(i+1))-abs(Peaks(i)) < 0
                DeleteP(j) = i+1;
                j = j+1;
                else
                DeleteP(j) = i;
                j = j+1;
                end
            end
        end
    end

    j = 1;
    for i = 1:numel(Troughs)-1
        if Troughs(i) > 0
            DeleteT(j) = i;
            j = j+1;
        elseif Troughs(i+1) > 0
            DeleteT(j) = i+1;
            j = j+1;
        else
            if abs(abs(Trough_Time(i+1))-abs(Trough_Time(i))) < 0.01
                if abs(Troughs(i+1))-abs(Troughs(i)) < 0
                DeleteT(j) = i+1;
                j = j+1;
                else
                DeleteT(j) = i;
                j = j+1;
                end
            end
        end
    end

    if DeleteP ~= 0
        Peaks(DeleteP) = [];
        Peak_Time(DeleteP) = [];
    end

    if DeleteT ~= 0
    Troughs(DeleteT) = [];
    Trough_Time(DeleteT) = [];
    end

    DeleteT = 0;
    DeleteP = 0;
    end



    %Find Peak and Trough Closest to Zero
    Closest_Peak = Peak_Time(abs(Peak_Time) == min(abs(Peak_Time)));
    Closest_Trough = Trough_Time(abs(Trough_Time) == min(abs(Trough_Time)));
    PS_M3_F2 = (Closest_Peak+Closest_Trough)/2;

%Mode 3 - Floor 3 Phase Shift
    LMax = islocalmax(Floor_3_M3); %Locate local maxima
    Peaks = Floor_3_M3(LMax); %Find y-axis value of maxima
    Peak_Time = F23_M3_Time(LMax); %Find x-axis value of maxima

    LMin = islocalmin(Floor_3_M3); %Locate local maxima
    Troughs = Floor_3_M3(LMin); %Find y-axis value of maxima
    Trough_Time = F23_M3_Time(LMin); %Find x-axis value of maxima

    %Find Minima and Maxima and Filter Noise
    for k = 1:5

    j = 1;
    for i = 1:numel(Peaks)-1
        if Peaks(i) < 0
            DeleteP(j) = i;
            j = j+1;
        else
            if abs(abs(Peak_Time(i+1))-abs(Peak_Time(i))) < 0.01
                if abs(Peaks(i+1))-abs(Peaks(i)) < 0
                DeleteP(j) = i+1;
                j = j+1;
                else
                DeleteP(j) = i;
                j = j+1;
                end
            end
        end
    end

    j = 1;
    for i = 1:numel(Troughs)-1
        if Troughs(i) > 0
            DeleteT(j) = i;
            j = j+1;
        elseif Troughs(i+1) > 0
            DeleteT(j) = i+1;
            j = j+1;
        else
            if abs(abs(Trough_Time(i+1))-abs(Trough_Time(i))) < 0.01
                if abs(Troughs(i+1))-abs(Troughs(i)) < 0
                DeleteT(j) = i+1;
                j = j+1;
                else
                DeleteT(j) = i;
                j = j+1;
                end
            end
        end
    end

    if DeleteP ~= 0
        Peaks(DeleteP) = [];
        Peak_Time(DeleteP) = [];
    end

    if DeleteT ~= 0
    Troughs(DeleteT) = [];
    Trough_Time(DeleteT) = [];
    end

    DeleteT = 0;
    DeleteP = 0;
    end



    %Find Peak and Trough Closest to Zero
    Closest_Peak = Peak_Time(abs(Peak_Time) == min(abs(Peak_Time)));
    Closest_Trough = Trough_Time(abs(Trough_Time) == min(abs(Trough_Time)));
    PS_M3_F3 = (Closest_Peak+Closest_Trough)/2;


%% Table Outputs

Mode = {'Mode 1';'Mode 1';'Mode 1';'Mode 2';'Mode 2';'Mode 2';'Mode 3';'Mode 3';'Mode 3'};
Floor = {'Floor 1';'Floor 2';'Floor 3';'Floor 1';'Floor 2';'Floor 3';'Floor 1';'Floor 2';'Floor 3'};
Amplitude = [Amp_M1_F1; Amp_M1_F2; Amp_M1_F3; Amp_M2_F1; Amp_M2_F2; Amp_M2_F3; Amp_M3_F1; Amp_M3_F2; Amp_M3_F3];

Column_Name_Amp = {'Mode';'Floor';'Amplitude'};

Amp_Table = table(Mode, Floor, Amplitude,'VariableNames',Column_Name_Amp)

Phase_Shift = [PS_M1_F1; PS_M1_F2; PS_M1_F3; PS_M2_F1; PS_M2_F2; PS_M2_F3; PS_M3_F1; PS_M3_F2; PS_M3_F3];

Column_Name_PS = {'Mode';'Floor';'Phase Shift'};

PS_Table = table(Mode, Floor, Phase_Shift,'VariableNames',Column_Name_PS)

%% Amplification
Amp_M1_F1_F2 = Amp_M1_F2/Amp_M1_F1
Amp_M1_F1_F3 = Amp_M1_F3/Amp_M1_F1

Amp_M2_F1_F2 = Amp_M2_F2/Amp_M2_F1
Amp_M2_F1_F3 = Amp_M2_F3/Amp_M2_F1

Amp_M3_F1_F2 = Amp_M3_F2/Amp_M3_F1
Amp_M3_F1_F3 = Amp_M3_F3/Amp_M3_F1


%% Error Analysis

%Human Error Analysis - Floors 1 & 2
M1 = [8.7, 9.0, 8.9, 8.9];
M2 = [13.5, 13.5, 13.6, 13.8];
M3 = [25, 25, 24.5, 24];

%Mean
Mean_M1 = mean(M1);
Mean_M2 = mean(M2);
Mean_M3 = mean(M3);

%Standard Deviation
SD_M1 = std(M1);
SD_M2 = std(M2);
SD_M3 = std(M3);

%Confidence Interval (95%)
CI_M1 = SD_M1*1.96;
CI_M2 = SD_M2*1.96;
CI_M3 = SD_M3*1.96;

%Lower Bounds
LB_M1 = Mean_M1 - CI_M1;
LB_M2 = Mean_M2 - CI_M2;
LB_M3 = Mean_M3 - CI_M3;

%Upper Bounds
UB_M1 = Mean_M1 + CI_M1;
UB_M2 = Mean_M2 + CI_M2;
UB_M3 = Mean_M3 + CI_M3;