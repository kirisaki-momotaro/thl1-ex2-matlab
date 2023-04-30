clear all
close all


T=0.001;
over=10;
Ts=T/over;
A=4;



 

%A.1
%plot SRRC with a=0,0.5,1
[phi_0, t] = srrc_pulse(T, over, A, 0);

plot(t,phi_0)


[phi_0_5, t] = srrc_pulse(T, over, A, 0.5);
hold on;
plot(t,phi_0_5)
grid on;


[phi_1, t] = srrc_pulse(T, over, A, 1);
hold on;
plot(t,phi_1)
grid on;
title('SRRC a = 0,0.5,1' )

legend('a=0','a=0.5','a=1')

%A.2
figure()
Nf=1024;
Fs = 1/Ts;               % sampling frequency
freq = (-Fs/2:Fs/Nf:Fs/2-1/Nf); % zero-centered frequency range
%fft SRRC 0
fftshift_SRRC_0 = fftshift(fft(phi_0,Nf)*Ts);
power_fftshift_SRRC_0 = abs(fftshift_SRRC_0).^2;     % zero-centered power

plot(freq,power_fftshift_SRRC_0)
grid on;


%fft SRRC 0.5
fftshift_SRRC_0_5 = fftshift(fft(phi_0_5,Nf)*Ts);
power_fftshift_SRRC_0_5 = abs(fftshift_SRRC_0_5).^2;     % zero-centered power
hold on;
plot(freq,power_fftshift_SRRC_0_5)


%fft SRRC 1
fftshift_SRRC_1 = fftshift(fft(phi_1,Nf)*Ts);
power_fftshift_SRRC_1 = abs(fftshift_SRRC_1).^2;     % zero-centered power
hold on;
plot(freq,power_fftshift_SRRC_1)
title('Spectral density SRRC a = 0,0.5,1' )

legend('a=0','a=0.5','a=1')

%semilogy 

%semilogy 
figure()
semilogy(freq,power_fftshift_SRRC_0)
hold on;
grid on;
semilogy(freq,power_fftshift_SRRC_0_5)
hold on;
grid on;
semilogy(freq,power_fftshift_SRRC_1)
hold on;
grid on;
legend('a=0','a=0.5','a=1')
title('semilogy')


figure()

semilogy(freq,power_fftshift_SRRC_0)
hold on;
grid on;
semilogy(freq,power_fftshift_SRRC_0_5)
hold on;
grid on;
semilogy(freq,power_fftshift_SRRC_1)
hold on;
grid on;


cb=T/(100000) +0*freq;
ca=T/(1000) +0*freq;
semilogy(freq,ca)
semilogy(freq,cb)
legend('a=0','a=0.5','a=1','ca','cb')
title('semilogy with cut frequencies')

%A3 b 



power_fftshift_SRRC_0_5( power_fftshift_SRRC_0_5 <= cb ) = 0;
power_fftshift_SRRC_0( power_fftshift_SRRC_0 <= cb ) = 0;
power_fftshift_SRRC_1( power_fftshift_SRRC_1 <= cb ) = 0;
figure()
semilogy(freq,power_fftshift_SRRC_0)
hold on;
grid on;
semilogy(freq,power_fftshift_SRRC_0_5)
hold on;
grid on;
semilogy(freq,power_fftshift_SRRC_1)
hold on;
grid on;



semilogy(freq,cb)
title('c = 10^{-8}')
hold on;
grid on;
legend('a=0','a=0.5','a=1','c')

%A3 b 


power_fftshift_SRRC_0_5( power_fftshift_SRRC_0_5 <= ca )=0;
power_fftshift_SRRC_0( power_fftshift_SRRC_0 <= ca ) = 0;
power_fftshift_SRRC_1( power_fftshift_SRRC_1 <= ca ) = 0;
figure()
semilogy(freq,power_fftshift_SRRC_0)
hold on;
grid on;
semilogy(freq,power_fftshift_SRRC_0_5)
hold on;
grid on;
semilogy(freq,power_fftshift_SRRC_1)
hold on;
grid on;


semilogy(freq,ca)
title('c = 10^{-6}')
hold on;
grid on;

legend('a=0','a=0.5','a=1','c')




