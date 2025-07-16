% =========================================================================

% Digital defocus Aberration Interference (DAbI) simulation code 
% Written by Shi Zhao and Haowen Zhou
% @ Caltech Biophotonics Laboratory | http://biophot.caltech.edu/
% The source code is licensed under GPL-3. 
% Version: March, 15th, 2025
%
% Reference: 
% @misc{zhou2025DAbI,
%       title={Digital defocus aberration interference for 
%               automated optical microscopy}, 
%       author={Haowen Zhou and Shi Zhao and Yujie Fan and Zhenyu Dong 
%               and Oumeng Zhang and Viviana Gradinaru and Changhuei Yang},
%       year={2025},
%       eprint={2507.10867},
%       archivePrefix={arXiv},
%       primaryClass={physics.optics},
%       url={https://arxiv.org/abs/2507.10867}, 
% }

% Project Page: https://hwzhou2020.github.io/DAbI-Web/
% GitHub Repo: https://github.com/hwzhou2020/DAbI

% =========================================================================

clearvars;  clc; close all;

%% Check dependancies 
if ~any(any(contains(struct2cell(ver), 'Image Processing Toolbox')))
    error('Image processing toolbox is required.');
end

%% Add dependent function paths
addpath(genpath('Subfunctions_Simulation'));
addpath(genpath('Subfunctions_DAbI'));

%% Function operators 
F  = @(x) fftshift(fft2(ifftshift(x))); % Fourier transform
Ft = @(x) ifftshift(ifft2(fftshift(x))); % Inverse Fourier transform
logamp = @(x) log10(abs(x)+1); % Log operation for spectrum visualization

%% Simulation parameters setting

defocus         = [200];   % Defocus distances (um) (Symbol: ∆z)
                            % Can have multiple values for a loop
                            
% General imaging parameters
mag             = 20;       % System magnification        
na_illu         = 0.39;     % Illumination NA （Symbol: sinβ<na_obj）   
na_obj          = 0.40;     % Objective lens NA
dpix_c          = 6.5;      % Camera pixel pitch (um) (Symbol: e)
wavelength      = 0.520;    % Wavelength (um)
downsampling    = 2;        % Downsampling rate for data generation
object_type     = "NSCLC";  % Options: "NSCLC" | "USAF"
obj_size        = 3072;     % Ground truth object size (<=3072)

% illumination settings
angle2use       = [pi/8, 3*pi/8];   
                            % The two azimuthal illumination angles 
                            % (Symbol: θ1 and θ2)

% Aberrations
zernike_coeffs  = zeros(1, 34);
                            % Coefficients of nomarlized Zernike polynomials 
                            % The sequence is based on Noll Index
                            % https://en.wikipedia.org/wiki/Zernike_polynomials
zernike_coeffs(4) = 0;      % Defocus aberration is already added in  
                            % "defocus". This term is set to 0

% GPU usage for reconstruction
use_GPU         = false;

% Visualization and saving
save_raw        = true;     % Save raw intensity image data
save_dir        = './Simulated_Data'; % Save directory
show_LEDs       = false;    % Visualize LED illumination angle in NA space
show_obj        = false;    % Visualize ground truth object
show_abe        = false;    % Visualize system aberration
show_raw        = false;    % Visualize the two intensity images (I1,I2)
show_fringe     = false;    % Visualize DAbI fringes (summed spectrum)

% Algorithm parameters (Only at the last section: Perform DAbI)
sub_pixel_resolve = true;   % Whether to use sub-pixel precision 
                            % to locate the valley positions
precision       = 0.01;    % Iteration precision (stopping criterion)
max_iter        = 5;        % Maximum iteration number

% Automatic settings based on above parameters
% Make save directory
if save_raw
    if ~exist(save_dir,"dir")
        mkdir(save_dir)
    end
end
% Check GPU
if gpuDeviceCount("available") == 0
    useGPU = false;
    disp('GPU not found. CPU is used instead.')
end
if ~any(any(contains(struct2cell(ver), 'Parallel Computing Toolbox')))
    useGPU = false;
    disp('Parallel Computing Toolbox is NOT installed. CPU is used instead.');
end

%% Calculate the illumination NA for each LED

% Illumination angles in NA space 
% illumination NA = sqrt( na_calib(:,1).^2 + na_calib(:,2).^2 )
na_calib = zeros(2,2);
na_calib(:,1) = na_illu .* cos(angle2use) ;
na_calib(:,2) = na_illu .* sin(angle2use) ;


% Visualization LED illumination angle in NA space
if show_LEDs
    figure(1001); hold on;
    plot(na_calib(:,1), na_calib(:,2), 'go', 'MarkerFaceColor', 'g', ...
        'MarkerSize', 12);
    title('LED illumination angle in NA space');
    axis square; grid on;
    xlabel('NA_x');
    ylabel('NA_y');
    xlim([-na_obj,na_obj])
    ylim([-na_obj,na_obj])
end


%% Object generation
if object_type == "USAF" % Pure Phase object
    im = double(imread('USAF-pc200nm.png'));
    % Crop center
    obj_phase = im(size(im,1)/2-obj_size/2+1:size(im,1)/2+obj_size/2, ...
        size(im,1)/2-obj_size/2+1:size(im,1)/2+obj_size/2,1)/255;
    obj_phase =obj_phase*pi/4;
    obj = ones(size(im(:,:,1))).*exp(1i*obj_phase);
elseif object_type == 'NSCLC'
    load("NSCLC.mat");
    % Crop center
    obj = obj(size(obj,1)/2-obj_size/2+1:size(obj,1)/2+obj_size/2, ...
        size(obj,1)/2-obj_size/2+1:size(obj,1)/2+obj_size/2,1);
end

% Visualize ground truth object
if show_obj
    figure(1002)
    subplot(1,3,1)
    imagesc(abs(obj));colormap gray; axis square;
    title('Designed Object Amplitude')
    axis off;
    subplot(1,3,2)
    imagesc(angle(obj));colormap gray; axis square;
    title('Designed Object Phase')
    axis off;
    subplot(1,3,3)
    imagesc(log(1+abs(fftshift(fft2(obj)))));colormap gray; axis square;
    title('Designed Object Spectrum')
    axis off;
end

%% Two raw intensity image generation (I1, I2)
spsize = (dpix_c)/mag;      % Sample plane pixel size
N = obj_size/downsampling;  % Size of the raw intensity image

for iz = 1:length(defocus)
    imlow_defocus = zeros(N,N,size(na_calib,1));
    if mod(N,2) == 1
        du = 1/(spsize*(N-1));
    else
        du = 1/(N*spsize);
    end
    % Illumination NA shift in pixel
    vled = na_calib(:,1)/wavelength;
    uled = na_calib(:,2)/wavelength;
    idx_u = round(uled/du);
    idx_v = round(vled/du);
    Ns = [];
    Ns(:,:,1) = -idx_v';
    Ns(:,:,2) = -idx_u'; 
    
    k0 = 2*pi/wavelength;
    kmax = pi/(spsize);
    dkx = kmax/(N/2);
    dky = kmax/(N/2);
    [kx, ky] = meshgrid(linspace(-kmax,kmax-dkx,N),linspace(-kmax,kmax-dky,N));
    % Amplitude of coherent transfer function (CTF)
    mask = (kx.^2+ky.^2<=(k0*na_obj)^2).*ones(N,N);
    % Added aberration (excluding defocus)
    phi_aberration = add_aberration_zernike(kx, ky, zernike_coeffs,k0*na_obj).*mask;
    
    % Visualize system aberration
    if show_abe
        figure(1003);imagesc(phi_aberration);
        axis off;axis image;
        title('Aberration at pupil phase');
        clim([-pi pi]);colormap(cNeoAlbedo);
    end
    
    % Coherent transfer function (CTF), including defocus
    H01 = mask.*(exp(-1i*real(sqrt(k0^2-kx.^2-ky.^2))*defocus(iz))).* ...
        (exp(-1i * phi_aberration));
    % Ground truth object spectrum
    OBJ = F(obj);
    
    % Generate two raw image
    for m = 1:size(na_calib,1)
        kxc = round((N*downsampling+1)/2-Ns(1,m,1));
        kyc = round((N*downsampling+1)/2-Ns(1,m,2));
        
        kxl = round(kxc-N/2);
        kyl = round(kyc-N/2);
        kxh = round(kxc+N/2-1);
        kyh = round(kyc+N/2-1);
        LOW_OBJ = OBJ(kxl:kxh,kyl:kyh) .* H01;
        low_obj = Ft(LOW_OBJ)/(downsampling^2);
        imlow_defocus(:,:,m)  = abs(low_obj).^2;
    end

    % Visualize the two intensity images (I1,I2)
    if show_raw
        figure(10041);
        imshow(imlow_defocus(:,:,1),[]);
        title('I1');colormap gray;axis equal;
        figure(10042);
        imshow(imlow_defocus(:,:,2),[]);
        title('I2');colormap gray;axis equal;
    end
    
    % Save raw intensity image data
    if save_raw
        file_name = object_type + sprintf('_z_%d_NA%03.0f_mag%d_DAbI.mat', defocus(iz), na_obj * 100, mag);
        saveName = fullfile(save_dir, file_name);
        save(saveName,'wavelength','na_illu','na_obj','mag','dpix_c','na_calib',...
                  'imlow_defocus','phi_aberration','-v7.3'); 
    end
    
    % Visualize fringes in the summed spectrum
    if show_fringe
        patch_size = size(imlow_defocus,1);
        xc = ceil(patch_size / 2) + 1;
        yc = ceil(patch_size / 2) + 1;
        freqXY_calib_show = (dpix_c / mag * patch_size * na_calib) / wavelength +[xc,yc] ;
        freqXY_calib_show = freqXY_calib_show(:, [2, 1]);
        idx_r = na_obj / (1 / (patch_size * dpix_c / mag)) / wavelength;
        Fimlow = mat2gray(log(1+abs(fftshift(fft2(double(1*imlow_defocus(:,:,1))+1*double(imlow_defocus(:,:,2)))))));
        figure(1005);
        imshow(Fimlow,[]);
        colormap(gray);hold on;
        viscircles(freqXY_calib_show(1, :), idx_r, 'LineStyle', '--', 'EdgeColor', 'r');
        drawnow;
        viscircles(freqXY_calib_show(2, :), idx_r, 'LineStyle', '--', 'EdgeColor', 'r');
        drawnow;
    end
    
    %% Perform DAbI
    patch_size = size(imlow_defocus,1);

    % Redefine spatial frequency paramters
    k0 = 2 * pi / wavelength;
    dk = 2 * pi * mag / (dpix_c * patch_size);
    kmax = pi * mag / dpix_c;
    [kx_coord, ky_coord] = meshgrid( linspace(-kmax, kmax - dk, patch_size), ...
                                 linspace(-kmax, kmax - dk, patch_size) );
    
    % Mimicing one-time calibration for system aberration
    phi_aberration = add_aberration_zernike(kx_coord, ky_coord, zernike_coeffs,k0*na_obj);
    
    % Compute CTF
    CTF = single((kx_coord.^2 + ky_coord.^2) < (k0 * na_obj)^2);
    CTF = max(CTF, circshift(rot90(CTF, 2), mod(patch_size, 2) * [1, 1]));
    CTF = CTF.*exp(-1i*phi_aberration);
    
    % Illumination angle in pixelated frequency space 
    freqXY_calib = (dpix_c / mag * patch_size * na_calib) / wavelength;

    % load variables to GPU
    if use_GPU
        imlow_defocus = gpuArray(single(imlow_defocus));
        freqXY_calib = single(freqXY_calib);
        CTF = gpuArray(single(CTF));
    else
        imlow_defocus = single(imlow_defocus);
        freqXY_calib = single(freqXY_calib);
        CTF = single(CTF);
    end
    

    % !! ---- DAbI start from here ---- !!
    z_defocus = findDefocus_DAbI(imlow_defocus, freqXY_calib, CTF, ...
        na_illu, na_obj, wavelength, mag, dpix_c, ...
        sub_pixel_resolve,use_GPU, precision, max_iter);

    fprintf('DAbI Estimated z: %.4f\n', z_defocus);
    fprintf('Ground Truth z:   %.4f\n', defocus(iz));

end


