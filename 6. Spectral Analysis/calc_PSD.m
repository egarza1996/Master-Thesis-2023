function [DSA, freq] = calc_PSD(signal, NFFT, sr, window, overlap)
% created by Julian Ostertag / . Matthias Kreuzer /. Anesthesiology r.d.I
% 05.05.2023
% %%%%%%%%%%%%%%%all rights reserved
% 
% CALC_PSD_DSA_NAN Calculates the power spectral density (PSD) of a multi-channel time series signal
%
% [DSA, freq] = calc_PSD(signal, NFFT, sr, window, overlap)
%
% Inputs:
% - signal: A 2D matrix of size number of time points x number of channels, representing the multi-channel time series signal.
% - NFFT: The number of frequency bins to use when calculating the PSD. This determines the frequency resolution of the resulting PSD.
% - sr: The sample rate of the signal, in Hz.
% - window: The length of the window to use for calculating the PSD, in seconds.
% - overlap: The length of the overlap between windows, in seconds.
%
% Outputs:
% - DSA: A 3D matrix of size number of frequency bins x number of channels x number of windows, where each frequency bin represents the power at a specific frequency.
% - freq: A frequency vector of size number of frequency bins x 1, giving the corresponding frequencies for each element in the PSD.
%
% Calculate the shift in samples between windows
shift = window - overlap;
shift = shift * sr; % Convert shift from seconds to samples
% Calculate the starting time points for each window
start_timepoints = [1, shift:shift:(size(signal, 1) - window *sr)];
num_windows = length(start_timepoints);
% Preallocate DSA matrix
DSA = nan(floor(NFFT/2)+1, size(signal, 2), num_windows);
freq = 0:(sr/NFFT):(sr/2);
% Calculate PSD for each window and channel
for ch = 1:size(signal, 2)
for i = 1:num_windows
start = start_timepoints(i);
stop = start + window * sr - 1;
cur_EEG = signal(start:stop, ch);
% If current window has no NaN values, calculate its PSD
if ~any(isnan(cur_EEG))
[cur_PSD, freq] = pwelch(cur_EEG, [], [], NFFT, sr, 'onesided');
DSA(:, ch, i) = cur_PSD;
end
end
end
