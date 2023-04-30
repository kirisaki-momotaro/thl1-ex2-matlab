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

%A2
N=100;
b = (sign(randn(N, 1)) + 1)/2;
X = bits_to_2PAM(b);

X_delta = 1/Ts * upsample(X, over);
figure(1)
%X_delta_time = 0:N*over-1;
X_delta_time = 0:Ts:N*Ts*over-Ts;
plot(X_delta_time,X_delta);
grid on;
title('original signal')




%create signal to be sent by sender
signal = conv(X_delta,phi)*Ts;
signal_t = [X_delta_time(1)+t(1):Ts:X_delta_time(end)+t(end)];
figure(2)
plot(signal_t,signal);
grid on;
title('modulated signal')









