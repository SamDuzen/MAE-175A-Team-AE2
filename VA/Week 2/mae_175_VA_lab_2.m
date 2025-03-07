w2 = 9.5 * 2 *pi;
w1 = 13.8 *2 *pi;
w3 = 24.5 *2 *pi;
b1 = 0.005;
b2 = 0.003;
b3 = 0.01;
K = 1;

num = [1 2*b1*w1 w1^2];
den = conv([1 2*b2*w2 w2^2], [1 2*b3*w3 w3^2]);
Hij = K*w2^2*w3^2/w1^2*tf(num,den)


myf = logspace(0,2,500);
[m,p]=bode(Hij,2*pi*myf);
subplot(2,1,1),semilogx(myf,20*log10(abs(squeeze(m))));
hold on
plot(freq,ch2)
subplot(2,1,2),semilogx(myf,squeeze(p))
hold on 
plot(freq,phase)

