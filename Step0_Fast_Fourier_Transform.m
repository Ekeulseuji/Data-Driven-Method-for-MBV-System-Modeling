%% already got the measured data（MBV_X, MBV_t）
% Need FFT for frequency ω, amplitude A
clear all;
close all;
clc

set(0, 'DefaultTextFontName', 'Times New Roman')

%% real X and t
mbv_X = readtable('Real_X.csv'); %  X = [ x y ] 
mbv_t = readtable('Real_t.csv'); % t = time points


%% frequency ω

% sampled points
Fs = 199;

% real data
S = mbv_X(1:Fs,1);

% corresponding real t
t = mbv_t(1:Fs,1);

% number of sampling points 
% when sampling in the time domain N points, the FFT results in N frequency bins
Nf = size(t,1);  


%% FFT
figure;

% table S to array S
S_array = table2array(S);
S_double = double(S_array);

% fft( ), do Fast Fourier Transform on signal S
FFT = fft(S_double, Nf);

% FFT = fft(S,Nf); 

% The magnitude of the FFT transformation results is stored in the Spectrum. 
% The Spectrum represents the frequency values obtained by computing the magnitude of the FFT results, 
% which are used to display the signal spectrum after FFT transformation.
Spectrum = (abs(FFT));

% Display the original FFT magnitude values (the signal spectrum with many frequency bins)
plot(Spectrum); 
title('The signal spectrum obtained after the data undergoes FFT magnitude');

% Calculate the actual frequency values corresponding to each frequency point in the FFT results.
% Fs is the sampling frequency, for example, 500, and Nf is the number of points used in the FFT transformation, which is the number of t.
% index-1 is to count from 0

% After performing the FFT calculation, the first frequency point in the spectrum result (the bin with index 0) 
% corresponds to the average value or offset in the signal, making it convenient to calculate the direct current component at frequency 0 (i.e., calculating the offset).

% Frequency represents the frequency values on the frequency axis, 
% used to indicate the actual frequencies corresponding to each frequency point in the spectrum Spectrum.

% index= 1:Nf;
% Frequency = (index - 1) * Fs /Nf;             

% Perform correction on the FFT transformation results 
% (convert to actual frequency values where N*Fs/2 corresponds to w = pi).

% Due to the fact that the transformation matrix in MATLAB's FFT is not a unitary matrix, 
% the signal is amplified after the FFT transformation. 
% 
% Dividing the matrix by 1/sqrt(N) turns it into a unitary matrix to correct this amplification effect.
% Therefore, after the FFT transformation, the result needs to be divided by N/2 to correct this "amplification" effect.
% The final result is the actual spectrum values, after correction.
% Actual_Specrum = 2/Nf * fft(S_double,Nf);

% the DC component should be divided by N for correction

figure;
plot(Frequency, abs(Actual_Specrum));
title('real intensity - frequency (correction before abs)');

% correction after abs()
figure;
Spectrum = Spectrum / (Nf /2); % Spectrum dividng by Nf/2（multiplying by 2/Nf）
Spectrum(1) = Spectrum(1) /2; % DC/N

% After FFT：Fn=(n-1)*Fs/N；(N starting from 1)
% real intensity - frequency (correction after abs)
plot(Frequency(1:Nf/2), Spectrum(1:Nf/2),'Color', [0.18434 0.3098 0.3098 ], 'linewidth', 1.5); 
set(gca, 'FontName', 'Times New Roman')

figure(1)
plot(Frequency(1:50),Spectrum(1:50),'Color', [0.18434 0.3098 0.3098 ], 'linewidth', 1.5);
set(gca, 'FontName', 'Times New Roman')
ylabel('Intensity', 'FontName', 'Times New Roman')
xlabel('Frequency', 'FontName', 'Times New Roman')
set(gca,'FontName','Times New Roman','FontSize',17)
figure(2)
% 0~N*Fs/2
plot(Frequency(1:Nf/2),Spectrum(1:Nf/2),'Color', [0.18434 0.3098 0.3098 ], 'linewidth', 1.5);    
ylabel('Intensity', 'FontName', 'Times New Roman')
xlabel('Frequency', 'FontName', 'Times New Roman')
set(gca,'FontName','Times New Roman','FontSize',17)

% save the data
save('sMBV.mat')
% 
% csvwrite('MBV_sampled_X.csv', S)
% csvwrite('MBV_sampled_t.csv', t)
% csvwrite('MBV_sampled_amplitude.csv', Spectrum(1:Nf/2))
% csvwrite('MBV_sampled_frequency.csv', Frequency(1:Nf/2))