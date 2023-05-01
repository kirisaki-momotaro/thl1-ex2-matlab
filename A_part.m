clear all
close all

%initial parameters
T=0.01;
over=10;
Ts=T/over;
A=4;
a=0.5;


 

%A.1
%create SRRC pulse
[phi, t] = srrc_pulse(T, over, A, a);



%FFT SRRC
figure(20)
Nf=2048; %number of samples
Fs = 1/Ts;               % sampling frequency
freq = (-Fs/2:Fs/Nf:Fs/2-1/Nf); % zero-centered frequency range
%fft SRRC
fftshift_SRRC = fftshift(fft(phi,Nf)*Ts);
power_fftshift_SRRC = abs(fftshift_SRRC).^2;     % zero-centered power spectral density

semilogy(freq,power_fftshift_SRRC)
title('energy spectral density of SRRC pulse')
grid on;

%A2 2PAM
N=100; %number of symbols
b = (sign(randn(N, 1)) + 1)/2; 
X = bits_to_2PAM(b); %random symbols

X_delta = 1/Ts * upsample(X, over); %upsample
figure(1)
%X_delta_time = 0:N*over-1;
X_delta_time = 0:Ts:N*Ts*over-Ts;
plot(X_delta_time,X_delta);
grid on;
title('2PAM symbols waveform')




%create modulated signal 
signal = conv(X_delta,phi)*Ts; %convolute symbols waveform with SRRC pulse
signal_t = [X_delta_time(1)+t(1):Ts:X_delta_time(end)+t(end)];
figure(2)
plot(signal_t,signal);
grid on;
title('modulated instance of 2PAM waveform')


%A3

%fft of modulated waveform
figure(3)
fftshift_signal = fftshift(fft(signal,Nf)*Ts);
power_fftshift_signal = abs(fftshift_signal).^2;     % zero-centered power

semilogy(freq,power_fftshift_signal)
title('Spectral density of modulated 2PAM waveform instance')
grid on;



%approximate spectral density by adding up numerous different experiments
%to get an approximation 
K=500; %number of experiments
power_fftshift_signal_sum=zeros(2048,1);
for i = 1:K
    b = (sign(randn(N, 1)) + 1)/2;
    X = bits_to_2PAM(b);

    X_delta = 1/Ts * upsample(X, over);    
    X_delta_time = 0:Ts:N*Ts*over-Ts;

    %create signal to be sent by sender
    signal = conv(X_delta,phi)*Ts;
    signal_t = [X_delta_time(1)+t(1):Ts:X_delta_time(end)+t(end)]; 


    %fft signal
   
    fftshift_signal = fftshift(fft(signal,Nf)*Ts);
    power_fftshift_signal = abs(fftshift_signal).^2;     % zero-centered power
    power_fftshift_signal_sum=power_fftshift_signal_sum+power_fftshift_signal; %add up all waveforms
    
end

%plot spectral density approximation
figure(4)
power_fftshift_signal_sum_normal=power_fftshift_signal_sum/K; %divide by K to normalize
semilogy(freq,power_fftshift_signal_sum_normal)
grid on;
hold on;


%theoretical spectral density according to provided equation

b = (sign(randn(N, 1)) + 1)/2;
X = bits_to_2PAM(b);
theoretical_spectral_density=((var(X)^2)/T)*power_fftshift_SRRC
semilogy(freq,theoretical_spectral_density)
title('Avg spectral density of modulated 2PAM waveform')
legend('approx','theoretical')






%A4 4PAM
N=100;
b = randi(4,(N/2)-1,1); %generate uniformely numbers 1-4
X = bits_to_4PAM(b);

X_delta = 1/Ts * upsample(X, over);
figure(6)

X_delta_time = 0:Ts:((N/2)-1)*Ts*over-Ts;
plot(X_delta_time,X_delta);
grid on;
title('original waveform 4PAM symbols')

%create signal to be sent by sender
signal = conv(X_delta,phi)*Ts;
signal_t = [X_delta_time(1)+t(1):Ts:X_delta_time(end)+t(end)];
figure(7)
plot(signal_t,signal);
grid on;
title('modulated instance of 4PAM waveform')

%fft signal
figure(8)
fftshift_signal = fftshift(fft(signal,Nf)*Ts);
power_fftshift_signal = (abs(fftshift_signal).^2)/0.5;     % zero-centered power

semilogy(freq,power_fftshift_signal)
title('Spectral density of modulated 4PAM waveform instance')
grid on;

%approximate spectral density
K=500;
power_fftshift_signal_sum=zeros(2048,1);
for i = 1:K
    b = randi(4,(N/2)-1,1); %generate uniformely numbers 1-4
    X = bits_to_4PAM(b);

    X_delta = 1/Ts * upsample(X, over);
    X_delta_time = 0:Ts:((N/2)-1)*Ts*over-Ts;

    %create signal to be sent by sender
    signal = conv(X_delta,phi)*Ts;
    signal_t = [X_delta_time(1)+t(1):Ts:X_delta_time(end)+t(end)]; 


    %fft signal
   
    fftshift_signal = fftshift(fft(signal,Nf)*Ts);
    power_fftshift_signal = (abs(fftshift_signal).^2)/0.5;     % zero-centered power
    power_fftshift_signal_sum=power_fftshift_signal_sum+power_fftshift_signal;
    
end

%plot spectral density approximation
figure(9)
power_fftshift_signal_sum_normal=power_fftshift_signal_sum/K;
semilogy(freq,power_fftshift_signal_sum_normal)
grid on;
hold on;


%theoretical spectral density according to provided equation
b = randi(4,(N/2)-1,1); %generate uniformely numbers 1-4
X = bits_to_4PAM(b);
theoretical_spectral_density=((var(X)^2)/T)*power_fftshift_SRRC
semilogy(freq,theoretical_spectral_density)
title('Avg spectral density of modulated 4PAM waveform')
legend('approx','theoretical')




%A5 repeat A3 with over,T *=2
T=0.02;
over=20;
Ts=T/over;
A=4;
a=0.5;

%create NEW SRRC
[phi, t] = srrc_pulse(T, over, A, a);

%NEW FFT SRRC
figure()
Nf=2048;
Fs = 1/Ts;               % sampling frequency
freq = (-Fs/2:Fs/Nf:Fs/2-1/Nf); % zero-centered frequency range
%fft SRRC
fftshift_SRRC = fftshift(fft(phi,Nf)*Ts);
power_fftshift_SRRC = abs(fftshift_SRRC).^2;     % zero-centered power

semilogy(freq,power_fftshift_SRRC)
title('energy spectral density of new SRRC pulse Tnew=2T')
grid on;
%create signal to be sent by sender
signal = conv(X_delta,phi)*Ts;
signal_t = [X_delta_time(1)+t(1):Ts:X_delta_time(end)+t(end)];
figure()
plot(signal_t,signal);
grid on;
title('modulated 2PAM signal Tnew=2T')
%fft signal
figure(10)
fftshift_signal = fftshift(fft(signal,Nf)*Ts);
power_fftshift_signal = abs(fftshift_signal).^2;     % zero-centered power

semilogy(freq,power_fftshift_signal)
title('Spectral density of modulated 2PAM waveform instance Tnew=2T')
grid on;



%approximate spectral density
K=500;
power_fftshift_signal_sum=zeros(2048,1);
for i = 1:K
    b = (sign(randn(N, 1)) + 1)/2;
    X = bits_to_2PAM(b);

    X_delta = 1/Ts * upsample(X, over);    
    X_delta_time = 0:Ts:N*Ts*over-Ts;

    %create signal to be sent by sender
    signal = conv(X_delta,phi)*Ts;
    signal_t = [X_delta_time(1)+t(1):Ts:X_delta_time(end)+t(end)]; 


    %fft signal
   
    fftshift_signal = fftshift(fft(signal,Nf)*Ts);
    power_fftshift_signal = abs(fftshift_signal).^2;     % zero-centered power
    power_fftshift_signal_sum=power_fftshift_signal_sum+power_fftshift_signal;
    
end

%plot spectral density approximation
figure(11)
power_fftshift_signal_sum_normal=power_fftshift_signal_sum/K;
semilogy(freq,power_fftshift_signal_sum_normal)
grid on;
hold on;




%theoretical spectral density according to provided equation
b = (sign(randn(N, 1)) + 1)/2;
X = bits_to_2PAM(b);
theoretical_spectral_density=((var(X)^2)/T)*power_fftshift_SRRC
semilogy(freq,theoretical_spectral_density)
title('Avg Spectral density of modulated 2PAM waveform Tnew=2T')
legend('approx','theoretical')



%save images
% FolderName = ('C:\Users\chris\Desktop\THL I\ex2\thl1-ex2-matlab\images');   % using my directory
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = num2str(get(FigHandle, 'Number'));
%   set(0, 'CurrentFigure', FigHandle);
% %   saveas(FigHandle, strcat(FigName, '.png'));
%   saveas(FigHandle, fullfile(FolderName,strcat(FigName, '.png'))); % specify the full path
% end




