function [WigVil, fvec, tvec] = chunk_wvd(x, fs, Nc)
%----------------------------------------------
% INPUT
%----------------------------------------------
% x         : Signal 
% fs        : Sampling Frequency
% Nc > 1    : Number of Chunks
%
%----------------------------------------------
% OUTPUT
%----------------------------------------------
% WigVil    : Chunkwise Pseudo Wigner-Ville distribution
% fvec      : Frequency bins
% tvec      : Time bins converted to actual time 

Xchunk = floor(length(x)/Nc); % data points per chunk (integer)
wvd_list = zeros(Xchunk, Xchunk*2, Nc); 
fvals = zeros(Xchunk, Nc);
tvals = zeros(Xchunk*2, Nc);
tol = 1e-6; 
x = hilbert(x);  % force analytic signal (diagnostic)
for k = 1:Nc
   block = x(1+(k-1)*Xchunk:k*Xchunk); 
   [power, fvals(:,k), tvals(:,k)] = wvd(block,fs); 
   wvd_list(:,:,k) = power; 
end
fvec = fvals(:,1); 
tvec = tvals(end,:); 
WigVil =  wvd_list(:,:,1);
for i = 1:Nc-1
    WigVil = horzcat(WigVil,wvd_list(:,:,i+1)); 
    tvec(i+1) = tvec(i) + tvals(end,i+1); 
end 
WigVil(:) = WigVil(:)./max(WigVil(:)); % normalize
WigVil(:) = 10*log10(WigVil(:)+tol); % convert to dB
