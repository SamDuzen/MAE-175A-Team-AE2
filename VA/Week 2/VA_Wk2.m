clear all; close all; clc;

%% Floor 1 and Floor 2

%Read Experimental Files
    %Floor 1 and Floor 2
    F12 = csvread('F12_Wk2.csv',21);
        F12_Freq = F12(:,1);
        F12_Ch1_Mag = F12(:,2);
        F12_Ch2_Mag = F12(:,3);
        F12_Ch2_Phase = F12(:,4);

%Modeling Parameters
w1 = 13.3548 * 2*pi; %Mode 2
w2 = 9.0985 * 2*pi; %Mode 1
w3 = 24.6776 * 2*pi; %Mode 3

b1 = 0.01; %Mode 2
b2 = 0.01; %Mode 1
b3 = 0.01; %Mode 3

%Transfer Function
Ci = (w2^2*w3^2)/(w1^2); %Gain
num = Ci*[1, 2*b1*w1, w1^2];
den = conv([1,2*b2*w2,w2^2],[1,2*b3*w3,w3^2]);
H12 = tf(num,den); %Transfer function from floor 2 to floor 1

%Plotting
myf = logspace(0,2,500);
[m,p] = bode(H12,2*pi*myf);
figure()
subplot(2,1,1)
    hold on
    semilogx(myf,20*log10(abs(squeeze(m))))
    plot(F12_Freq,F12_Ch2_Mag)
    xlim([5,50])
    title('Floor 1 & Floor 2 - Magnitude')
    grid on
    legend('Simulated','Experimental','Location','best')
subplot(2,1,2)
    hold on
    semilogx(myf,squeeze(p))
    plot(F12_Freq,F12_Ch2_Phase)
    xlim([5,50])
    title('Floor 1 & Floor 2 - Phase')
    grid on
    legend('Simulated','Experimental','Location','northwest')

%% Floor 1 and Floor 3

%Read Experimental Files
    %Floor 1 and Floor 3
    F13 = csvread('F13_Wk2.csv',21);
        F13_Freq = F13(:,1);
        F13_Ch1_Mag = F13(:,2);
        F13_Ch2_Mag = F13(:,3);
        F13_Ch2_Phase = F13(:,4);

%Modeling Parameters
w1 = 13.3548 * 2*pi; %Mode 2
w2 = 9.0985 * 2*pi; %Mode 1
w3 = 24.6776 * 2*pi; %Mode 3

b1 = 0.01; %Mode 2
b2 = 0.01; %Mode 1
b3 = 0.01; %Mode 3

%Transfer Function
Ci = (w2^2*w3^2)/(w1^2); %Gain
num = Ci*[1, 2*b1*w1, w1^2];
den = conv([1,2*b2*w2,w2^2],[1,2*b3*w3,w3^2]);
H12 = tf(num,den); %Transfer function from floor 2 to floor 1

%Plotting
myf = logspace(0,2,500);
[m,p] = bode(H12,2*pi*myf);
figure()
subplot(2,1,1)
    hold on
    semilogx(myf,20*log10(abs(squeeze(m))))
    plot(F13_Freq,F13_Ch2_Mag)
    xlim([5,50])
    title('Floor 1 & Floor 3 - Magnitude')
    grid on
    legend('Simulated','Experimental','Location','best')
subplot(2,1,2)
    hold on
    semilogx(myf,squeeze(p))
    plot(F13_Freq,F13_Ch2_Phase)
    xlim([5,50])
    title('Floor 1 & Floor 3 - Phase')
    grid on
    legend('Simulated','Experimental','Location','northwest')

%% Floor 2 and Floor 3

%Read Experimental Files
    %Floor 2 and Floor 3
    F23 = csvread('F23_Wk2.csv',21);
        F23_Freq = F23(:,1);
        F23_Ch1_Mag = F23(:,2);
        F23_Ch2_Mag = F23(:,3);
        F23_Ch2_Phase = F23(:,4);

%Modeling Parameters
w1 = 13.3548 * 2*pi; %Mode 2
w2 = 9.0985 * 2*pi; %Mode 1
w3 = 24.6776 * 2*pi; %Mode 3

b1 = 0.01; %Mode 2
b2 = 0.01; %Mode 1
b3 = 0.01; %Mode 3

%Transfer Function
Ci = (w2^2*w3^2)/(w1^2); %Gain
num = Ci*[1, 2*b1*w1, w1^2];
den = conv([1,2*b2*w2,w2^2],[1,2*b3*w3,w3^2]);
H12 = tf(num,den); %Transfer function from floor 2 to floor 1

%Plotting
myf = logspace(0,2,500);
[m,p] = bode(H12,2*pi*myf);
figure()
subplot(2,1,1)
    hold on
    semilogx(myf,20*log10(abs(squeeze(m))))
    plot(F23_Freq,F23_Ch2_Mag)
    xlim([5,50])
    title('Floor 2 & Floor 3 - Magnitude')
    grid on
    legend('Simulated','Experimental','Location','best')
subplot(2,1,2)
    hold on
    semilogx(myf,squeeze(p))
    plot(F23_Freq,F23_Ch2_Phase)
    xlim([5,50])
    title('Floor 2 & Floor 3 - Phase')
    grid on
    legend('Simulated','Experimental','Location','northwest')

