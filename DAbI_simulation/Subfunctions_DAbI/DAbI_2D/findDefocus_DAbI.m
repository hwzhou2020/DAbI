function [z] = findDefocus_DAbI(imlow_defocus, k_illu, CTF, na_illu, na_obj, wavelength, mag, dpix_c, sub_pixel_resolve, use_GPU, precision, max_iter)
    % findDefocus_DAbI is the main function to estimate the defocus distance 
    % using digital Defocus Aberration Inference (DAbI) methods.
    
    % Inputs:
    % - imlow_defocus:   Captured two images from LED illuminations.
    % - k_illu:          Illumination wavevectors k_ (already converted to pixel shift in k space).
    % - CTF:             System Coherent Transfer Function (Including other aberrations that we want to compensate).
    % - na_illu:         Illumination NA (Symbol: sinβ).
    % - na_obj:          NA of the objective.
    % - wavelength:      Illumination wavelength (um).
    % - mag:             System magnification.
    % - dpix_c:          Camera pixel pitch (um).
    % - zernike_abe:     Zernike coefficients for system aberrations.
    % - sub_pixel_resolve: Whether to use sub-pixel precision to locate the peak.
    % - use_GPU:         Whether to enable GPU acceleration.
    % - precision        Iteration precision (stopping criterion)
    % - max_iter         Maximum iteration number
    % Output:
    % - z: Estimated defocus distance (in um).
    
    %% Part 1: Initialization and system parameters calculation
    k_illu = round(k_illu);
    [xsize, ysize, nBrightField] = size(imlow_defocus);
    
    k0 = 2 * pi / wavelength;
    dk = 2 * pi * mag / (dpix_c * xsize);
    kmax = pi * mag / dpix_c;
    mask2use = abs(CTF) > 1e-2; % CTF support (k0*na_obj)
    HDoF = wavelength / na_obj^2 /2; % System half Depth of Field (DoF).
    
    % Define thresholds for different defocus levels (in um)
    z_threshold_low = HDoF * 6; % Small
    z_threshold_medium = HDoF * 15;  % Medium
    z_threshold_high = HDoF * 100;  % Large
    
    % Virtual defocus distance added for low-defocus cases (Symbol: z_v in Supplementary Note 1.4)
    z_virtual = HDoF * 50;
    
    % This is for inputting more two images
    % By default, adjacent image pairs in imlow_defocus are used to estimate ∆z via the DAbI method. 
    % If desired, specific image pairs can be selected by modifying the 'multiplex_order' variable below.
    multiplex_order = 1:size(imlow_defocus, 3);
    imlow_defocus_mo = single(imlow_defocus(:, :, multiplex_order));
    mycoord_mo = k_illu(multiplex_order, :);
    
    
    nBrightField_multi = floor(nBrightField/2); % # of image pairs to use.
    BF_fft = fftshift(fftshift(fft2(imlow_defocus_mo), 1), 2); % Fourier transform of the captured intensity images.
    
    if use_GPU
        imStack_FFT = gpuArray(zeros(xsize, ysize, nBrightField_multi,'single')); % Summed FFT image stacks
        measurement_mask = gpuArray(false(xsize, ysize, nBrightField_multi)); % Overlapping regions S
    else
        imStack_FFT = zeros(xsize, ysize, nBrightField_multi,'single'); % Summed FFT image stacks
        measurement_mask = false(xsize, ysize, nBrightField_multi); % Overlapping regions S
    end
    
    % Create overlap masks and sum spectra for each image pair 
    for i = 1:nBrightField_multi
        idx1 = 2*i - 1;
        idx2 = 2*i;
    
        kxc1 = mycoord_mo(idx1, 1);
        kyc1 = mycoord_mo(idx1, 2);
        mask1 = circshift(mask2use, [-kxc1, -kyc1]);
    
        kxc2 = mycoord_mo(idx2, 1);
        kyc2 = mycoord_mo(idx2, 2);
        mask2 = circshift(mask2use, [-kxc2, -kyc2]);
        
        measurement_mask(:, :, i) = mask1 & mask2;
        imStack_FFT(:,:,i) = BF_fft(:,:,idx1) + BF_fft(:,:,idx2);
    end
    
    % Compute angles and centers of two illumination angles
    angles = atan2(-mycoord_mo(:, 1), -mycoord_mo(:, 2));
    angles_sep = abs(angles(2:2:end) - angles(1:2:end))/2;
    k_center = (mycoord_mo(1:2:end, :) + mycoord_mo(2:2:end, :)) / 2;
    angles_rad = atan2(-k_center(:, 1), -k_center(:, 2)) / pi * 180;
    
    
    %% Part 2: Rough estimation of defocus level using FFT 
    
    % Crop a sub rectangular region in set S to estimate the defocus level
    crop_ratio_in = 1 / 6;
    crop_ratio_out = 1 / 10;
    
    % Crop S and rotate to rotate it to align the fringes horizontally
    overlap_spectrum = imrotate(abs(imStack_FFT(:, :, 1)) .* measurement_mask(:, :, 1), angles_rad(1), 'bilinear', 'crop');
    
    % Find inscribed rectangle withn S
    [~, colIdx] = find(overlap_spectrum);
    col_min = min(colIdx); col_max = max(colIdx);
    col_idx = round(col_min + (col_max - col_min) * crop_ratio_in);
    nonzero_rows = find(overlap_spectrum(:, col_idx));
    row_min = min(nonzero_rows);
    row_max = max(nonzero_rows);
    row_crop = min(size(overlap_spectrum,1)/2-row_min, row_max-size(overlap_spectrum,1)/2);
    row1 = size(overlap_spectrum,1)/2 - row_crop +1;
    row2 = size(overlap_spectrum,1)/2 + row_crop;
    col1 = floor(col_min + (col_max - col_min) * crop_ratio_in);
    col2 = floor(col_max - (col_max - col_min) * crop_ratio_out);
    
    % Crop subregion based on the SNR
    SNR_tolerance = 3.5; % May need to be modified according to image SNR
    spectrum_sum = abs(fftshift(fft2(imlow_defocus_mo(:, :, 1)))) + abs(fftshift(fft2(imlow_defocus_mo(:, :, 2))));
    spectrum_sum_rot = imrotate(spectrum_sum .* measurement_mask(:, :, 1), angles_rad(1), 'bilinear', 'crop');
    mask_rot = imrotate(measurement_mask(:, :, 1), angles_rad(1), 'bilinear', 'crop');
    spectrum_avg = sum(spectrum_sum_rot,1)./(sum(mask_rot,1)+eps).*(sum(mask_rot,1)>0);
    spectrum_avg_cut = log(spectrum_avg(col1:col2)+1);
    sig_max = max(spectrum_avg_cut);
    spectrum_avg_cut(spectrum_avg_cut>sig_max-SNR_tolerance) = 0;
    col2_relative_idx = find(spectrum_avg_cut, 1, 'first');
    if ~isempty(col2_relative_idx)
        col2 = col1 + col2_relative_idx;
    end
    overlap_crop = overlap_spectrum(row1:row2, col1:col2);
    
    if na_obj < 0.5
        insymmetric_tolerance = 3;
    else
        insymmetric_tolerance = 15;
    end
    
    % Apply FFT and extract central FFT line for peak analysis
    spectrum_crop = abs(fftshift(fft2(overlap_crop)));
    central_line = mean(spectrum_crop(:, floor(end/2)-2 : floor(end/2)+2), 2);
    [~, locs] = findpeaks(central_line, 'MinPeakProminence', max(central_line) * 0.05); % Factors may need to be changed
    center = floor(length(central_line)/2+1);
    locs   = locs(locs == center | ismembertol(2*center - locs, locs, insymmetric_tolerance, 'DataScale', 1)); % Filter out all insymmetric peaks
    
    % Compute thresholds for delta_k (Separation of two peaks in FFT originated from the fringes) for different defocus levels
    Nk = xsize * dpix_c ./ [z_threshold_low, z_threshold_medium, z_threshold_high] ./ mag ./ na_illu ./ sin(angles_sep(1))/2;
    delta_k_thesholds = 1 ./ Nk .* size(spectrum_crop, 1);
    
    locs = locs(abs(locs - center) > delta_k_thesholds(1));
    locs = [locs; center];
    
    if length(locs) < 3 % If insufficient peaks, relax peak detection criteria
        [~, locs] = findpeaks(central_line, 'MinPeakProminence', max(central_line) * 0.03); % Factors may need to be changed
        locs   = locs(locs == center | ismembertol(2*center - locs, locs, insymmetric_tolerance, 'DataScale', 1)); % Filter out all insymmetric peaks
        locs = locs(abs(locs - center) > delta_k_thesholds(1));
        locs = [locs; center];
    
        if length(locs) < 3 % If still no preaks are found, the defocus level should be really low
            defocus_level = 0;
        else
            cand = locs(abs(locs-center)>1 & ismembertol(2*center-locs, locs, insymmetric_tolerance, 'DataScale', 1));
            [~,k] = max( central_line(cand) + central_line(2*center - cand) );
            locs  = [center;cand(k); 2*center - cand(k);];       
            delta_k = abs((locs(3) - locs(2)) / 2);
            % Decide the defocus level based on delta_k
            if delta_k > delta_k_thesholds(3)
                defocus_level = 3;
            elseif delta_k > delta_k_thesholds(2) && delta_k < delta_k_thesholds(3)
                defocus_level = 2;
            else
                defocus_level = 1;
            end
        end
    else
        cand = locs(abs(locs-center)>1 & ismembertol(2*center-locs, locs, insymmetric_tolerance, 'DataScale', 1));
        [~,k] = max( central_line(cand) + central_line(2*center - cand) );
        locs  = [center;cand(k); 2*center - cand(k);];       
        delta_k = abs((locs(3) - locs(2)) / 2);
        % Decide the defocus level based on delta_k
        if delta_k > delta_k_thesholds(3)
            defocus_level = 3;
        elseif delta_k > delta_k_thesholds(2) && delta_k < delta_k_thesholds(3)
            defocus_level = 2;
        else
            defocus_level = 1;
        end
    end
    
    
    %% Part 3-1: Accurate defocus estimation; Condition 1: Large defocus level
    if defocus_level == 3 || defocus_level == 2 
        z_virtual = 0;
        % Estimate defocus magnitude
        Nk_est_multi = (1 / delta_k * size(spectrum_crop, 1));
        z_est =  xsize * dpix_c / Nk_est_multi / mag / na_illu / sin(angles_sep(1))/2;
    
        cropH = row2 - row1 + 1;
        cropW = col2 - col1 + 1;
        if use_GPU
            overlap_crop_all = gpuArray(zeros(cropH, cropW, nBrightField_multi,'single')); 
            imStack_FFT_shift = gpuArray(zeros(xsize, ysize, nBrightField/2,'single'));
        else
            overlap_crop_all = zeros(cropH, cropW, nBrightField_multi,'single'); 
            imStack_FFT_shift = zeros(xsize, ysize, nBrightField/2,'single');
        end
    
        % Extract rotated spectrum from all image pairs
        for n_multi = 1:nBrightField_multi
            rotated = imrotate(abs(imStack_FFT(:, :, n_multi)) .* measurement_mask(:, :, n_multi), ...
                               angles_rad(n_multi), 'bilinear', 'crop');
            overlap_crop_all(:, :, n_multi) = rotated(row1:row2, col1:col2);
        end
    
        % Find an estimated value of defocus from the defocus aberration (Without high-order compensation)
        z_d = findDefocus_DAbI_FFT(overlap_crop_all, angles_sep, xsize, dpix_c,z_est,z_threshold_low, mag, na_illu,na_obj, sub_pixel_resolve, z_virtual,defocus_level);
        if isnan(z_d)
            z_d = z_est;
        end
    
        % Use virtual defocus aberration method to decide defocus direction (Supplementary Note 1.4)
        [kx_coord, ky_coord] = meshgrid( linspace(-kmax, kmax - dk, xsize), ...
                                         linspace(-kmax, kmax - dk, xsize) );
        sq_term = k0^2 - kx_coord.^2 - ky_coord.^2;
        sq_term(sq_term < 0) = 0;                   
        sqrt_term = sqrt(sq_term);              
    
        z_change = z_d * 0.4; % Small defocus change.
        % Estimate the fringe densities if added z_change on both side.
        Nk_est_up  = xsize * dpix_c / (z_d-z_change) / mag / na_illu / angles_sep(1)/2;
        Nk_est_down  = xsize * dpix_c / (z_d+z_change) / mag / na_illu / angles_sep(1)/2;
        CTF_vir = CTF .* exp(1i * (sqrt_term - k0) * (z_change)); 
        z_f = findDefocus_DAbI_Direction( ...
                     Nk_est_up, ...
                     Nk_est_down, ...
                     CTF_vir, ...
                     mycoord_mo(1:2,:), ...
                     imlow_defocus_mo(:,:,1:2), ...
                     measurement_mask(:,:,1), ...
                     angles_rad(1),z_d, row1, row2, col1, col2);
    
    
    
        % From the expansion of defocus kernel 
        % sqrt(k0^2-kx^2-ky^2) = k0 - (kx^2+ky^2)/(2*k0) + high-order terms,
        % the defocus aberration have basic terms and high-order terms.    
        
        % We need to compensate for the higher-order terms in defocus to 
        % enhance DAbI accuracy. 
        % Here, we use iterations until (1) the higher-order aberration terms 
        % from defocus kernel are fully compensated, and 
        % (2) the defocus distance converges within the predefined precision.
        
        % Start compensating for these higher-order aberration terms
        z_0 = z_f;
        z_1 = z_f;
        count = 1;
        while (abs(z_0-z_1) > precision && count < max_iter) || count==1
            CTF_def = CTF .* exp(-1i*sqrt_term*z_1 - 1i*(0.5*(kx_coord.^2+ky_coord.^2)./k0.*z_1));   
            for i = 1:nBrightField_multi
                idx1 = 2*i - 1;
                idx2 = 2*i;
            
                kxc1 = mycoord_mo(idx1, 1);
                kyc1 = mycoord_mo(idx1, 2);
                spectrum1 = circshift(CTF_def, [-kxc1, -kyc1]) .* BF_fft(:,:,idx1);
            
                kxc2 = mycoord_mo(idx2, 1);
                kyc2 = mycoord_mo(idx2, 2);
                spectrum2 = circshift(CTF_def, [-kxc2, -kyc2]) .* BF_fft(:,:,idx2);
                imStack_FFT_shift(:,:,i) = spectrum1 + spectrum2;
            end
            for n_multi = 1:nBrightField_multi
                rotated = imrotate(abs(imStack_FFT_shift(:, :, n_multi)) .* measurement_mask(:, :, n_multi), ...
                                   angles_rad(n_multi), 'bilinear', 'crop');
                overlap_crop_all(:, :, n_multi) = rotated(row1:row2, col1:col2);
            end
    
            % Recalculate the defocus distance after compensating the higher-order aberration terms
            z_2 = findDefocus_DAbI_FFT(overlap_crop_all, angles_sep, xsize, dpix_c,z_est,z_threshold_low, mag, na_illu,na_obj, sub_pixel_resolve, z_virtual,defocus_level);
            z_2 = abs(z_2)*sign(z_f);
            if isnan(z_2)
                if isnan(z_f)
                    z = z_est;
                else
                    z = z_f;
                end
                break
            else
                z_0 = z_1;
                z_1 = z_2;
                z = z_1;
            end
            count = count + 1;
        end
    
    
    %% Part 3-2: Accurate defocus estimation; Condition 2: Small defocus level
    else
        % If the defocus level is small, the fringes are not dense enough for us 
        % to accurately estimate defocus distance. 
        % Therefore, we employed virtual defocus aberration strategy: 
        % virtually adding a known, large defocus aberration to the each spectrum.
        % Then we can calculate defocus distanace from virtual defocus aberration
        % (Supplementary Note 1.4)
        z_est = 0;
        [kx_coord, ky_coord] = meshgrid( linspace(-kmax, kmax - dk, xsize), ...
                                     linspace(-kmax, kmax - dk, xsize) );
        sq_term = k0^2 - kx_coord.^2 - ky_coord.^2;
        sq_term(sq_term < 0) = 0;                   
        sqrt_term = sqrt(sq_term);   
        CTF_support = abs(CTF);
        CTF_vir = CTF_support .* exp(1i * (sqrt_term - k0) * z_virtual); 
        
        cropH = row2 - row1 + 1;
        cropW = col2 - col1 + 1;
    
        if use_GPU
            overlap_crop_all = gpuArray(zeros(cropH, cropW, nBrightField_multi,'single')); 
            imStack_FFT_shift = gpuArray(zeros(xsize, ysize, nBrightField/2,'single'));
        else
            overlap_crop_all = zeros(cropH, cropW, nBrightField_multi,'single'); 
            imStack_FFT_shift = zeros(xsize, ysize, nBrightField/2,'single');
        end
    
        % Add the known, large defocus aberration to the each spectrum
        for i = 1:nBrightField_multi
            idx1 = 2*i - 1;
            idx2 = 2*i;
        
            kxc1 = mycoord_mo(idx1, 1);
            kyc1 = mycoord_mo(idx1, 2);
            spectrum1 = circshift(CTF_vir, [-kxc1, -kyc1]) .* BF_fft(:,:,idx1);
        
            kxc2 = mycoord_mo(idx2, 1);
            kyc2 = mycoord_mo(idx2, 2);
            spectrum2 = circshift(CTF_vir, [-kxc2, -kyc2]) .* BF_fft(:,:,idx2);
            imStack_FFT_shift(:,:,i) = spectrum1 + spectrum2;
        end
    
        % Extract rotated spectrum from all image pairs
        for n_multi = 1:nBrightField_multi
            rotated = imrotate(abs(imStack_FFT_shift(:, :, n_multi)) .* measurement_mask(:, :, n_multi), ...
                               angles_rad(n_multi), 'bilinear', 'crop');
            overlap_crop_all(:, :, n_multi) = rotated(row1:row2, col1:col2);
        end
    
        % Find an estimated value of defocus (Physical + Virtual) from the defocus aberration (Without high-order compensation)
        z_d = findDefocus_DAbI_FFT(overlap_crop_all, angles_sep, xsize, dpix_c, z_est, z_threshold_low, mag, na_illu,na_obj, sub_pixel_resolve, z_virtual,defocus_level);
        
    
        % From the expansion of defocus kernel 
        % sqrt(k0^2-kx^2-ky^2) = k0 - (kx^2+ky^2)/(2*k0) + high-order terms,
        % the defocus aberration have basic terms and high-order terms.    

        % Compensate for the higher-order terms
        z_0 = z_d;
        z_1 = z_d;
        count = 1;
        while (abs(z_0-z_1) > precision && count < max_iter) || count==1
            CTF_def = CTF .* exp(-1i*sqrt_term*z_1 - 1i*(0.5*(kx_coord.^2+ky_coord.^2)./k0.*z_1)); 
            for i = 1:nBrightField_multi
                idx1 = 2*i - 1;
                idx2 = 2*i;
            
                kxc1 = mycoord_mo(idx1, 1);
                kyc1 = mycoord_mo(idx1, 2);
                spectrum1 = circshift(CTF_def, [-kxc1, -kyc1]) .* circshift(CTF_vir, [-kxc1, -kyc1]) .* BF_fft(:,:,idx1);
            
                kxc2 = mycoord_mo(idx2, 1);
                kyc2 = mycoord_mo(idx2, 2);
                spectrum2 = circshift(CTF_def, [-kxc2, -kyc2]) .* circshift(CTF_vir, [-kxc2, -kyc2]) .* BF_fft(:,:,idx2);
                imStack_FFT_shift(:,:,i) = spectrum1 + spectrum2;
            end
    
    
            % Extract rotated spectrum from all image pairs
            for n_multi = 1:nBrightField_multi
                rotated = imrotate(abs(imStack_FFT_shift(:, :, n_multi)) .* measurement_mask(:, :, n_multi), ...
                                   angles_rad(n_multi), 'bilinear', 'crop');
                overlap_crop_all(:, :, n_multi) = rotated(row1:row2, col1:col2);
            end
            % Recalculate the defocus distance after compensating the higher-order aberration terms
            z_2 = findDefocus_DAbI_FFT(overlap_crop_all, angles_sep, xsize, dpix_c,z_est, z_threshold_low, mag, na_illu,na_obj, sub_pixel_resolve, z_virtual,defocus_level);
            z_2 = abs(z_2)*sign(z_d);
            if isnan(z_2)
                    if isnan(z_d)
                        z = 0;
                    else
                        z = z_d;
                    end
                    break
            else
                z_0 = z_1;
                z_1 = z_2;
                z = z_1;
            end
            count = count + 1;
        end
        % Subtract z_virtual to get the real defocus distance together with direction
        z = -z_virtual + z;
       
    end

end



