%% Glossary     

% Feed per tooth (fz) = [0.05 0.075 0.1 0.125 0.15 0.175 0.2] ae/mm 
% Radial depth (ae) = [8 16 24 32 40] mm
% Diameter (D) = 80 mm
% Ratio of radial depth (ae/D)  = [0.1 0.2 0.3 0.4 0.5] 
% Cutting speed (vc) = 120 m/min

%% Prepare Data
clc; clear; close all; rng(67); 
load("preproc_DBS900_F01_ae05_fz0150.mat")

% Start 1 s earlier 
acc1f_cut=acc1f(round((T_begin-1)*fsd):round((T_end)*fsd));
acc2f_cut=acc2f(round((T_begin-1)*fsd):round((T_end)*fsd));
Td_cut=Td(round((T_begin-1)*fsd):round((T_end)*fsd))';

% Parameters 
factor = 10; 
fs = fsd/factor; % Hz 
win_dur = 1.64; % window duration (s) 
L = round(win_dur*fs); % window length
window = hann(L);
overlap = 0.9; % welch method 90% overlap
noverlap = round(L*overlap); 
NFFT = 4*L; % could try other values 
tol = 1e-6; % numerical stability for Power Spectrum normalization 
P = NFFT; 
Np = 5; % Number of peaks 

% Bandpass Filter 20-30Hz
d = designfilt('bandpassiir', 'FilterOrder', 6,'HalfPowerFrequency1', 20, ...
    'HalfPowerFrequency2', 30, 'SampleRate', fs);
band_acc1f = filtfilt(d, acc1f); % accelerator 1 full
band_acc2f= filtfilt(d, acc2f); % accelaerator 2 full 
band_acc1f_cut = filtfilt(d, acc1f_cut); % accelerator 1 reduced
band_acc2f_cut = filtfilt(d, acc2f_cut); % accelerator 2 reduced 

%% Spectral Estimates (full data)

% Accelerator 1 ------------
N1 = length(band_acc1f_cut); 
% Periodogram
[per_acc1f, ffgrid, fper_acc1f] = per_est(band_acc1f, fs, NFFT, Np);
% RELAX
fREL_acc1f  = sort( relax(band_acc1f(:),Np)*fs/(2*pi)); 

% Accelerator 2 ------------
N2 = length(band_acc2f_cut);
% Periodogram
[per_acc2f, ~, fper_acc2f] = per_est(band_acc2f, fs, NFFT, Np);
% RELAX
fREL_acc2f  = sort( relax(band_acc2f(:),Np)*fs/(2*pi));

%% Spectral Estimates (reduced data)

% Accelerator 1 ------------
% Periodogram
[per_acc1f_cut, ~, fper_acc1f_cut] = per_est(band_acc1f_cut, fs, NFFT, Np);
% RELAX
fREL_acc1f_cut  = sort( relax(band_acc1f_cut(:),Np)*fs/(2*pi));

% Accelerator 2 ------------
% Periodogram
[per_acc2f_cut, ~, fper_acc2f_cut] = per_est(band_acc2f_cut, fs, NFFT, Np);
% RELAX
fREL_acc2f_cut  = sort( relax(band_acc2f_cut(:),Np)*fs/(2*pi));

%% Results

% Full Data
fprintf('----------------------------------------------\n')
fprintf('\nFrequency Estimates Full Data\n')
fprintf('----------------------------------------------\n')
fprintf('\nPeriodogram Acc1:\n')
fprintf('%g ', fper_acc1f)
fprintf('\n')
fprintf('\nPeriodogram Acc2:\n')
fprintf('%g ', fper_acc2f)
fprintf('\n')
fprintf('\nRELAX Acc1:\n')
fprintf('%g ', fREL_acc1f)
fprintf('\n')
fprintf('\nRELAX Acc2:\n')
fprintf('%g ', fREL_acc2f)
fprintf('\n')
fprintf('----------------------------------------------\n')

% Reduced Data
fprintf('----------------------------------------------\n')
fprintf('\nFrequency Estimates Reduced Data\n')
fprintf('----------------------------------------------\n')
fprintf('\nPeriodogram Acc1:\n')
fprintf('%g ', fper_acc1f_cut)
fprintf('\n')
fprintf('\nPeriodogram Acc2:\n')
fprintf('%g ', fper_acc2f_cut)
fprintf('\n')
fprintf('\nRELAX Acc1:\n')
fprintf('%g ', fREL_acc1f_cut)
fprintf('\n')
fprintf('\nRELAX Acc2:\n')
fprintf('%g ', fREL_acc2f_cut)
fprintf('\n')
fprintf('----------------------------------------------\n')

%% Plots
% Periodograms
figure;
subplot(211)
semilogy(ffgrid, per_acc1f)
hold on; 
semilogy(ffgrid, per_acc2f)
hold off; 
title('Full Data Periodogram')
legend('Acc1', 'Acc2')
grid on; 
subplot(212)
semilogy(ffgrid, per_acc1f_cut)
hold on; 
semilogy(ffgrid, per_acc2f_cut)
hold off; 
title('Reduced Data Periodogram')
legend('Acc1', 'Acc2')
grid on; 
