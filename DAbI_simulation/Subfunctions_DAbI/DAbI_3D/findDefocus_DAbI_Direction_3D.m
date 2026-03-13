function [z] = findDefocus_DAbI_Direction_3D(Nk_up, Nk_down, CTF_vir, mycoord, imStack, overlap_mask, angle_rad, z,row1, row2, col1,col2)
% findDefocus_DAbI_Direction determines the sign of the defocus value (z)
% by comparing the fringe density by adding known defocus aberration on 
% both side. The direction where the fringes are denser should be the
% direction of defocus.

% This method is only used for large defocus distance.

% Inputs:
% - Nk_up / Nk_down: Estimated fringe density for z +/- ∆z cases.
% - CTF_vir:         Virtual, known defocus aberration.
% - mycoord:         Illumination wavevectors (pixel shift).
% - imStack:         Two raw images (2-LED acquisition).
% - overlap_mask:    Binary mask for overlapping fringe region S.
% - angle_rad:       Rotation angle to align fringes horizontally.
% - z:               Estimated magnitude of defocus.
% - row1, row2, col1, col2: Indices to crop the aligned region of interest.

% Output:
% - z: Final signed defocus (positive or negative), based on fringe asymmetry.

%% Step 1: Apply conjugate and non-conjugate virtual CTFs to input images
    % Apply to image 1
    kxc1 = mycoord(1,1); kyc1 = mycoord(1,2);
    recon_FT1 = fftshift(fft2(imStack(:,:,1)));
    spectrum1_A = circshift(conj(CTF_vir), [-kxc1 -kyc1]).*recon_FT1;
    spectrum1_B = circshift(CTF_vir, [-kxc1 -kyc1]).*recon_FT1;
    % Apply to image 2
    kxc2 = mycoord(2,1); kyc2 = mycoord(2,2);
    recon_FT2 = fftshift(fft2(imStack(:,:,2)));
    spectrum2_A = circshift(conj(CTF_vir), [-kxc2 -kyc2]).*recon_FT2;
    spectrum2_B = circshift(CTF_vir, [-kxc2 -kyc2]).*recon_FT2;
    % Combine to calculate summed spectra
    overlap_A = (spectrum1_A + spectrum2_A).*overlap_mask;
    overlap_B = (spectrum1_B + spectrum2_B).*overlap_mask;
    % Rotate to align fringes horizontally and crop region of interest
    rotated_A = imrotate(abs(overlap_A), angle_rad, 'bilinear', 'crop');
    overlap_spectrum_A = rotated_A(row1:row2, col1:col2);
    rotated_B = imrotate(abs(overlap_B), angle_rad, 'nearest', 'crop');
    overlap_spectrum_B = rotated_B(row1:row2, col1:col2);

    % overlap_spectrum_A and overlap_spectrum_B contains fringes of
    % differnet density, apply FFT and locate the estimated peak regions to
    % decide the defocus direction
    spectrum_FFT_A = abs(fftshift(fft2(overlap_spectrum_A)));
    spectrum_FFT_B = abs(fftshift(fft2(overlap_spectrum_B)));

    f_k_up  = floor(1/Nk_up*size(spectrum_FFT_A,1));
    kxc  = floor(size(spectrum_FFT_A,1)/2);
    kx_1 = kxc - f_k_up;
    kx_2 = kxc + f_k_up;
    kyc  = floor(size(spectrum_FFT_A,2)/2);
    freq_area_A1 = sum(sum(spectrum_FFT_A(kx_1-3:kx_1+3,kyc-2:kyc+3)));
    freq_area_A2 = sum(sum(spectrum_FFT_A(kx_2-3:kx_2+3,kyc-2:kyc+3)));
    freq_area_B1 = sum(sum(spectrum_FFT_B(kx_1-3:kx_1+3,kyc-2:kyc+3)));
    freq_area_B2 = sum(sum(spectrum_FFT_B(kx_2-3:kx_2+3,kyc-2:kyc+3)));

    f_k_down  = floor(1/Nk_down*size(spectrum_FFT_A,1));
    kx_3 = kxc - f_k_down;
    kx_4 = kxc + f_k_down;
    if f_k_down < kxc -3
        freq_area_A3 = sum(sum(spectrum_FFT_A(kx_3-3:kx_3+3,kyc-2:kyc+3)));
        freq_area_A4 = sum(sum(spectrum_FFT_A(kx_4-3:kx_4+3,kyc-2:kyc+3)));
        freq_area_B3 = sum(sum(spectrum_FFT_B(kx_3-3:kx_3+3,kyc-2:kyc+3)));
        freq_area_B4 = sum(sum(spectrum_FFT_B(kx_4-3:kx_4+3,kyc-2:kyc+3)));
        if freq_area_A1+freq_area_A2+freq_area_B3+freq_area_B4 >= freq_area_B1+freq_area_B2+freq_area_A3+freq_area_A4 
            z = z;
        else
            z = -z;
        end
    else
        if freq_area_A1+freq_area_A2>= freq_area_B1+freq_area_B2
            z = z;
        else
            z = -z;
        end
    end



end