function [est, ff, fpeaks] = per_est(x,fs, P, Np)
    %----------------------------------------------
    % INPUT
    %----------------------------------------------
    % x         : Signal 
    % fs        : Sampling frequency
    % P         : Zero-padding 
    % Np        : Number of peaks to find
    %----------------------------------------------
    % OUTPUT
    %----------------------------------------------
    % est       : Spectral estimate
    % ff        : Frequency grid
    % peaksval  : Power of peaks
    % fpeaks    : Location of peaks on frequency grid 

    N = length(x); 
    ff = (0:P-1).*(fs/P) - fs/2; 
    est = fftshift(abs(fft(x, P)).^2)/(N*fs);
    [peaksidx,~] = findpeaks(est, Np);
    fpeaks = ff(peaksidx) ;  fpeaks = sort(fpeaks); 
end
