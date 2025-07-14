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

% To be released


end
