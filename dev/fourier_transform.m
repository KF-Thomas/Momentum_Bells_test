Fs = 1e9;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 17001;             % Length of signal
t = (0:L-1)*T;        % Time vector

Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
stfig('fourier transform')
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')