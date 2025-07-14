clc;clear all; close all;
warning off curvefit:fit:convertingY

%% Set up parameters
folder = 'Data'; % Calibration Data to be load
folder_name = dir(folder);
folder_name = folder_name(3:end);
z_defocus_matrix = zeros(length(folder_name), 67);

% load system parameters
load("system_parameters.mat")

sub_pixel_resolve = true;
patch_size = 1536;
vis_spec = false;

use_GPU = true;
if gpuDeviceCount("available") == 0
    useGPU = false;
end

%% Use DAbI to estimate defocus distance 
for sample_idx = 1:length(folder_name)
    sample_name = folder_name(sample_idx).name;
    filepath = [folder filesep sample_name filesep];

    xc = ceil(patch_size / 2) + 1;
    yc = ceil(patch_size / 2) + 1;
    na_rp_cal = 1 / wavelength * na_obj / (1 / (patch_size * dpix_c / mag));

    freqXY_calib_show = (dpix_c / mag * patch_size * na_calib) / wavelength +[xc,yc] ;
    freqXY_calib_show = freqXY_calib_show(:, [2, 1]);
    freqXY_calib = (dpix_c / mag * patch_size * na_calib) / wavelength;

    k0 = 2 * pi / wavelength;
    dk = 2 * pi * mag / (dpix_c * patch_size);
    kmax = pi * mag / dpix_c;
    [kx_coord, ky_coord] = meshgrid( linspace(-kmax, kmax - dk, patch_size), ...
                                 linspace(-kmax, kmax - dk, patch_size) );
    
    phi_aberration = add_aberration_zernike(kx_coord, ky_coord, zernike_abe,k0*na_obj);
    CTF = single((kx_coord.^2 + ky_coord.^2) < (k0 * na_obj)^2);
    CTF = max(CTF, circshift(rot90(CTF, 2), mod(patch_size, 2) * [1, 1]));
    CTF = CTF.*exp(-1i*phi_aberration);


    count = 1;
    zlist = [0:5:50,60:10:430,455:25:855];
    for z_pos = zlist            
        curr_z = num2str(floor((z_pos)*1000),"%06.0f");
        filename = [filepath sample_name '_z_' curr_z '_nm.mat'];
        load(filename)
        imlow_defocus = imlow_defocus(1024-patch_size/2+1:1024+patch_size/2,1024-patch_size/2+1:1024+patch_size/2,:);
    
        if vis_spec
            idx_r = na_obj / (1 / (patch_size * dpix_c / mag)) / wavelength;
            Fimlow = mat2gray(log(1+abs(fftshift(fft2(double(1*imlow_defocus(:,:,1))+1*double(imlow_defocus(:,:,2)))))));
            figure(1284);
            colormap(gray);
            clf;
            imshow(Fimlow,[])
            hold on;
            viscircles(freqXY_calib_show(1, :), idx_r, 'LineStyle', '--', 'EdgeColor', 'r');
            viscircles(freqXY_calib_show(2, :), idx_r, 'LineStyle', '--', 'EdgeColor', 'b');
            drawnow;
        end


        if use_GPU
            imlow_defocus = gpuArray(single(imlow_defocus));
            freqXY_calib = single(freqXY_calib);
            CTF = gpuArray(single(CTF));
        else
            imlow_defocus = single(imlow_defocus);
            freqXY_calib = single(freqXY_calib);
            CTF = single(CTF);
        end
        tic
        z_defocus = findDefocus_DAbI(imlow_defocus, freqXY_calib, CTF, na_illu, na_obj, wavelength, mag, dpix_c, sub_pixel_resolve,use_GPU);
        toc
        z_defocus_matrix(sample_idx, count) = z_defocus;
        count = count + 1;
        
        disp(['Real Z = ',num2str(z_pos)]);
        if z_defocus == 0
            disp('Calculated Z = 0');
        else
            disp(['Calculated Z = ',num2str(z_defocus)]);
        end
        
    
    end
end




