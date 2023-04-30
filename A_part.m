clear all
close all


T=0.01;
over=10;
Ts=T/over;
A=4;
a=0.5;


 

%A.1
%create SRRC
[phi, t] = srrc_pulse(T, over, A, a);



%FFT SRRC
figure()
Nf=2048;
Fs = 1/Ts;               % sampling frequency
freq = (-Fs/2:Fs/Nf:Fs/2-1/Nf); % zero-centered frequency range
%fft SRRC
fftshift_SRRC = fftshift(fft(phi,Nf)*Ts);
power_fftshift_SRRC = abs(fftshift_SRRC).^2;     % zero-centered power

semilogy(freq,power_fftshift_SRRC)
grid on;










