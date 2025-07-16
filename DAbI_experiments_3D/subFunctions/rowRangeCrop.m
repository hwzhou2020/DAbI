function [row1_new, row2_new] = rowRangeCrop(overlap_crop,Nk_est,row1,row2)
%  rowRangeCrop: Refines vertical cropping range based on fitted fringe contrast.
%
%
%   This function uses high-pass filtering on the averaged vertical fringe profile
%   to detect fringe edges according to local contrast and refine the cropping window.
%
%   Inputs:
%       overlap_crop - Cropped image region
%       Nk_est       - Estimated fringes intensity (used to tune filters)
%       row1, row2   - Original vertical crop range in the full image
%
%   Outputs:
%       row1_new, row2_new - Refined crop range based on edge detection
    % Parameters
    alpha = 0.7; % Ratio of corner frequency to fringe frequency
    extFac = 1/30; % Margin factor to extend cropping range
    fs = 1;  % Sampling frequency (normalized to 1)
    
    % Create a blur to suppress noise and small details
    filter_size = max(floor(Nk_est/4)*2+1, 11);
    imblur = imgaussfilt(overlap_crop, 5,'FilterSize', filter_size);
    % Define high-pass filter based on estimated fringe frequency
    x = imblur(:,1);
    f0 = 1/Nk_est; % Fringe frequency                   
    fc = alpha * f0; % High-pass cutoff frequency
    Wn = fc / (fs/2);  % Normalized cutoff frequency              
    [b,a] = butter(4, Wn, 'high'); % 4th order Butterworth high-pass filter

    % Apply zero-phase high-pass filtering to detect edges
    x_hp  = filtfilt(b, a, double(-x(:))); 
    x_hp = x_hp(:);  
    N    = numel(x_hp);
    mid  = round(N / 2)+1;
    % Find peaks indicating content edges
    [pks, locs] = findpeaks(x_hp, 'MinPeakProminence', 0.3 * max(x_hp));

    % Fallback: if no peaks are found, return full original range
    if isempty(pks)
        row1_new = row1;
        row2_new = row2;
        return
    end
   
    % Find outermost peaks and define symmetric cropping range around center
    left_idx = min(locs);
    right_idx = max(locs);

    half_width = max([mid - left_idx, right_idx - mid]);
    half_width = half_width + round(N * extFac);

    row1_new = row1 + max(1, mid - half_width)-1;
    row2_new = row2 - (N - min(N, mid + half_width));

end
