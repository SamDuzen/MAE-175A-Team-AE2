%establish parameters
omega_2 = 9 * 2 *pi;
omega_1 = 13.8 *2 *pi;
omega_3 = 24.5 *2 *pi;
beta_1 = 0.005;
beta_2 = 0.003;
beta_3 = 0.01;
K = 1;

%retrieve data
freq = F13Wk2{:,1};
ch1_mag = F13Wk2{:,2};
ch2_mag = F13Wk2{:,3};
ch2_phase = F12Wk2{:,4};

%design sinusoidal model
num = [1 2*beta_1*omega_1 omega_1^2];
den = conv([1 2*beta_2*omega_2 omega_2^2], [1 2*beta_3*omega_3 omega_3^2]);
Hij = K*omega_2^2*omega_3^2/omega_1^2*tf(num,den);

%plot together
myf = logspace(0,2,500);
[m,p]=bode(Hij,2*pi*myf);
subplot(2,1,1),semilogx(myf,20*log10(abs(squeeze(m))));
hold on
plot(freq,ch2_mag)
subplot(2,1,2),semilogx(myf,squeeze(p))
hold on 
plot(freq,ch2_phase)
