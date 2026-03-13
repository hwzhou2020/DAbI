function z = findDefocus_DAbI_FFT(overlap_crop_all, angles_sep, xsize, dpix_c, z_est, z_threshold_low, mag, na_illu, sub_pixel_resolve, z_shift, defocus_level)
%   findDefocus_DAbI_FFT estimate defocus distance from the summed image spectrum fringes using FFT.
%
%   overlap_crop_all
%
%   Inputs:
%       overlap_crop_all    - Summed image spectra for different image pairs (With/without virtual modification and high-order compensation)
%       angles_sep          - Illumination angle separation
%       xsize               - Field of view size (pixels)
%       dpix_c              - Camera pixel size (microns)
%       z_threshold_low     - Lower bound of defocus search range
%       mag                 - Imaging system magnification
%       na_illu             - Illumination NA
%       sub_pixel_resolve   - Boolean: whether to refine fringe spacing to subpixel accuracy
%       z_shift             - Virtual defocus added
%       defocus_level       - Defocus level, 0,1,2,3 from small to large.
%
%   Output:
%       z                   - Estimated defocus distance (absolute value, microns)

z_values = [];  % Store z estimates from each multi-frame
nBrightField_multi = size(overlap_crop_all, 3);

for n_multi = 1:nBrightField_multi
    overlap_crop = overlap_crop_all(:,:,n_multi); % Extract current overlap region
    % Pad vertically to 512 rows (for stable performance)
    pad_size = 512 - size(overlap_crop, 1);
    if pad_size > 10000
        pad_before = floor(pad_size / 2);
        pad_after = pad_size - pad_before;
        top_pad = repmat(overlap_crop(1, :), pad_before, 1);
        bot_pad = repmat(overlap_crop(end, :), pad_after, 1);
        overlap_crop_pad = [top_pad; overlap_crop; bot_pad];
    else
        overlap_crop_pad = overlap_crop;
    end

    % Compute 2D FFT and extract central vertical line
    spectrum_pad = abs(fftshift(fft2(overlap_crop_pad)));
    central_line = mean(spectrum_pad(:, floor(end/2)+1-2 : floor(end/2)+1+2), 2);

    % Detect peaks with moderate prominence
    [~, locs] = findpeaks(central_line, 'MinPeakProminence', max(central_line) * 0.03);% Factors may need to be changed
    center = floor(length(central_line)/2 + 1);
    locs   = locs(locs == center | ismembertol(2*center - locs, locs, 3, 'DataScale', 1)); % Filter out all insymmetric peaks
    
    % Compute thresholds for delta_k (Separation of two peaks in FFT originated from the fringes) for z_threshold_low
    % delta_k_theshold_low is used to filter out all low-frequency peaks.
    Nk = xsize * dpix_c / z_threshold_low / mag / na_illu / sin(angles_sep(1))/2;
    delta_k_theshold_low = 1 / Nk * size(spectrum_pad, 1);

    % Compute thresholds for delta_k (Separation of two peaks in FFT originated from the fringes) for z_shift
    % delta_k_theshold_shift is used to filter out the contribution from the added virtual defocus.
    Nk = xsize * dpix_c / z_shift / mag / na_illu / sin(angles_sep(1))/2;
    delta_k_theshold_shift = 1 / Nk * size(spectrum_pad, 1);


    if z_shift == 0  %large defocus distance
        Nk = xsize * dpix_c / (z_est*0.8) / mag / na_illu / sin(angles_sep(1))/2;
        delta_k_theshold_down = 1 / Nk * size(spectrum_pad, 1);
        Nk = xsize * dpix_c / (z_est*1.2) / mag / na_illu / sin(angles_sep(1))/2;
        delta_k_theshold_up = 1 / Nk * size(spectrum_pad, 1);
        locs = locs(abs(locs - center) > delta_k_theshold_down & ...
            abs(locs - center) < delta_k_theshold_up);
    else    % small defocus distance
        if length(locs)>=5 & defocus_level==1
            lower_bound = center - delta_k_theshold_shift;
            upper_bound = center + delta_k_theshold_shift;
            
            locs = locs(abs(locs - center) > delta_k_theshold_low & ...
                        abs(locs - lower_bound) > 2 & abs(locs - upper_bound) > 2);
        else
            locs = locs(abs(locs - center) > delta_k_theshold_low);
        end
    end

    locs = [locs; center];
    

    if length(locs) < 3 % If insufficient peaks, relax peak detection criteria
        [~, locs] = findpeaks(central_line, 'MinPeakProminence', max(central_line) * 0.01);
        locs   = locs(locs == center | ismembertol(2*center - locs, locs, 3, 'DataScale', 1));
        if z_shift == 0     %large defocus distance
            Nk = xsize * dpix_c / (z_est*0.8) / mag / na_illu / sin(angles_sep(1))/2;
            delta_k_theshold_down = 1 / Nk * size(spectrum_pad, 1);
            Nk = xsize * dpix_c / (z_est*1.2) / mag / na_illu / sin(angles_sep(1))/2;
            delta_k_theshold_up = 1 / Nk * size(spectrum_pad, 1);
            locs = locs(abs(locs - center) > delta_k_theshold_down & ...
                abs(locs - center) < delta_k_theshold_up);
        else    % small defocus distance
            if length(locs)>=5
                lower_bound = center - delta_k_theshold_shift;
                upper_bound = center + delta_k_theshold_shift;
                
                locs = locs(abs(locs - center) > delta_k_theshold_low & ...
                            abs(locs - lower_bound) > 2 & abs(locs - upper_bound) > 2);
            else
                locs = locs(abs(locs - center) > delta_k_theshold_low);
            end
        end
        locs = [locs; center];
        
        if length(locs) >= 3 
            % Estimate delta_k from the most prominent symmetric pair
            cand = locs(abs(locs-center)>1 & ismembertol(2*center-locs, locs, 3, 'DataScale', 1));
            [~,k] = max( central_line(cand) + central_line(2*center - cand) );
            locs  = [center;cand(k); 2*center - cand(k);];       
            if sub_pixel_resolve
                loc_p_1 = locs(2);
                loc_p_2 = locs(3);
                loc_p_1_sub = subPixelFit(central_line, loc_p_1, 7, 0.01);
                loc_p_2_sub = subPixelFit(central_line, loc_p_2, 7, 0.01);
                delta_k = abs((loc_p_2_sub - loc_p_1_sub) / 2);
            else
                delta_k = abs((locs(3) - locs(2)) / 2);
            end
        else
            delta_k = 0;
        end
    else
        % Enough peaks found; estimate delta_k from the most prominent symmetric pair
        cand = locs(abs(locs-center)>1 & ismembertol(2*center-locs, locs, 3, 'DataScale', 1));
        [~,k] = max( central_line(cand) + central_line(2*center - cand) );
        locs  = [center;cand(k); 2*center - cand(k);];       
        if sub_pixel_resolve
            loc_p_1 = locs(2);
            loc_p_2 = locs(3);
            loc_p_1_sub = subPixelFit(central_line, loc_p_1, 7, 0.01);
            loc_p_2_sub = subPixelFit(central_line, loc_p_2, 7, 0.01);
            delta_k = abs((loc_p_2_sub - loc_p_1_sub) / 2);
        else
            delta_k = abs((locs(3) - locs(2)) / 2);
        end
    end

    % Compute z from delta_k using physical models
    if delta_k ~= 0
        Nk_est_multi = (1 / delta_k * size(spectrum_pad, 1));
        z_value =  xsize * dpix_c / Nk_est_multi / mag / na_illu / sin(angles_sep(1))/2;
        z_values = [z_values; z_value];
    end
end

z = mean(z_values);

end
