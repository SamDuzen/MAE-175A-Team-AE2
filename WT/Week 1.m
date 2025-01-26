close all; clear all; clc;

speed = [10,20,30,40,50,60,70,80];
pressure = [0.082,0.381,0.929,1.685,2.659,3.833,5.289,6.81];

p_si = 248.84.*pressure;

U = sqrt((2*p_si)/(1.196384077));

figure()
plot(speed,U)

percent = @(x) 1.251 + 1.471*(x);