%% Glossary     

% Feed per tooth (fz) = [0.05 0.075 0.1 0.125 0.15 0.175 0.2] ae/mm 
% Radial depth (ae) = [8 16 24 32 40] mm
% Diameter (D) = 80 mm
% Ratio of radial depth (ae/D)  = [0.1 0.2 0.3 0.4 0.5] 
% Cutting speed (vc) = 120 m/min

%% Prepare Data
clc; clear; close all; rng(67); 
load("preproc_DBS900_F01_ae05_fz0175.mat") % change data 

% Start 1 s earlier 
acc1f_cut=acc1f(round((T_begin-1)*fsd):round((T_end)*fsd));
acc2f_cut=acc2f(round((T_begin-1)*fsd):round((T_end)*fsd));
Td_cut=Td(round((T_begin-1)*fsd):round((T_end)*fsd))';

fs = fsd; % Hz 
win_dur = 1.64; % window duration (s) 
L = win_dur*fs; % window length
window = hann(L);
overlap = round(L*0.9); % welch method 90% overlap

% Bandpass 
d = designfilt('bandpassiir', 'FilterOrder', 6,'HalfPowerFrequency1', 20, ...
    'HalfPowerFrequency2', 30, 'SampleRate', fs);
NFFT = 4*L; % could try other values 
tol = 1e-6; % numerical stability for PS normalization 

%% Spectrograms (full data) 

% Bandpass filters 
band_acc1f = filtfilt(d, acc1f); % accelerator 1
band_acc2f= filtfilt(d, acc2f); % accelaerator 2

% Spectrogram accelerator 1
[s1,f1,t1] = spectrogram(band_acc1f, window, overlap, NFFT, fs);
p1 = (abs(s1).^2)./ max(abs(s1).^2); % normalize
SD_acc1f = 10*log10(p1+tol); % convert to dB
    
% Spectrogram accelerator 2
[s2,f2,t2] = spectrogram(band_acc2f, window, overlap, NFFT, fs);
p2 = (abs(s2).^2)./ max(abs(s2).^2); % normalize
SD_acc2f = 10*log10(p2+tol); % convert to dB


%% Spectrograms (reduced data) 
% Bandpass filter 20-30 Hz 
d = designfilt('bandpassiir', 'FilterOrder', 6,'HalfPowerFrequency1', 20, ...
    'HalfPowerFrequency2', 30, 'SampleRate', fs);
band_acc1f_cut = filtfilt(d, acc1f_cut); % accelerator 1
band_acc2f_cut = filtfilt(d, acc2f_cut); % accelerator 2

% Spectrogram accelerator 1
[s1c,f1c,t1c] = spectrogram(band_acc1f_cut, window, overlap, NFFT, fs);
p1c = (abs(s1c).^2)./ max(abs(s1c).^2); % normalize
SD_acc1f_cut = 10*log10(p1c+tol); % convert to dB
    
% Spectrogram accelerator 2
[s2c,f2c,t2c] = spectrogram(band_acc2f_cut, window, overlap, NFFT, fs);
p2c = (abs(s2c).^2)./ max(abs(s2c).^2); % normalize
SD_acc2f_cut = 10*log10(p2c+tol); % convert to dB

%% Plots 

% Acelerometer data
figure; 
subplot(211)
plot(Td,acc1f)
title('Total data, Data1')
subplot(212)
plot(Td,acc2f)
title('Total data, Data2')
xlabel('Time (s)')

% Renyi Entropy
figure(2)
subplot(111)
plot(TI,R1,'-b',TI,R2,'-r')
hold 
plot([T_begin-1 T_begin-1],[min(R1) max(R1)],'--g')
plot([T_end T_end],[min(R1) max(R1)],'--g')
hold
title('Renyi Data1-blue, Renyi Data2-red')
xlabel('Time (s)')

% Power Spectra (full)
figure(3)
subplot(211)
imagesc(t1, f1, SD_acc1f);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Power Spectrum acc1 (dB)');
shading interp;  
colormap("turbo");   
axis tight;
ylim([0 152])
subplot(212)
imagesc(t2, f2, SD_acc2f);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Power Spectrum acc2 (dB)');
shading interp; 
colormap("turbo");   
axis tight; 
ylim([0 152])

% Power Spectra (reduced)
figure(4)
subplot(211)
imagesc(t1c, f1c, SD_acc1f_cut);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Power Spectrum acc1-cut (dB)');
shading interp;  
colormap("turbo");   
axis tight;
ylim([0 152])
subplot(212)
imagesc(t2c, f2c, SD_acc2f_cut);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Power Spectrum acc2-cut (dB)');
shading interp; 
colormap("turbo");   
axis tight; 
ylim([0 152])
