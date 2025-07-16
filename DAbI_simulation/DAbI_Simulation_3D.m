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


%%
clear;close all;clc;

%% Check dependencies
% Package requirements:
% 1. Image Processing Toolbox
% 2. Parallel Computing Toolbox (Optional, for GPU computing)
if ~any(any(contains(struct2cell(ver), 'Image Processing Toolbox')))
    error('Image processing toolbox is required.');
end
if ~any(any(contains(struct2cell(ver), 'Parallel Computing Toolbox'))) && (gpuDeviceCount("available") ~= 0)
    error('Parallel Computing Toolbox is required.');
end

%% Add dependent function paths
addpath(genpath('Subfunctions_Simulation'));
addpath(genpath('Subfunctions_DAbI'));

%% Simulation parameters setting

defocus     = 200;        % unit: um. Defocus distances (Symbol: ∆z)

% General imaging parameters
mag         = 20;         % magnification of the objective lens
ps          = 6.5/mag;    % unit: um. Pixel size in object space
pixelSizeXY = ps;         % unit: um. Effective pixel size (lateral)
pixelSizeZ  = ps;         % unit: um. Allows uneven 3D grid, 
                          % i.e. pixelSizeZ can be different from pixelSizeXY
NA          = 0.40;       % numerical aperture (NA) of the imaging system
na_illu     = NA*0.98;       % Illumination NA （Symbol: na_illu = sinβ<NA）   
lambda      = 0.635;      % unit: um. wavelength of the illumination light

useGPU     = true;       % whether to use GPU for reconstruction 
if gpuDeviceCount("available") == 0
    useGPU = false;
end

% Simulated object parameters
xsize       = 2048;       % unit: pixel
ysize       = 2048;       % unit: pixel. Prefer to be the same as xsize


n_media     = 1.330;      % Background refractive index
dn          = 0.005;       % Refractive index variations 
samLayers   = 150;        % Number of sample layers to use 
                          % Maximum 340 for cell.mat
                          % samLayers * pixelSizeZ = sample thickness in um
zsize       = samLayers + 50;        
                          % 50 layers for media

% illumination settings
angle2use       = [0, pi/3];   
                            % The two azimuthal illumination angles 
                            % (Symbol: θ1 and θ2)


% Algorithm parameters (Only at the last section: Perform DAbI)
sub_pixel_resolve = true;   % Whether to use sub-pixel precision 
                            % to locate the valley positions
precision       = 0.2;     % Iteration precision (stopping criterion)
max_iter        = 5;        % Maximum iteration number

% Visualization and saving
save_raw        = false;     % Save raw intensity image data
save_dir        = './Simulated_Data_3D'; % Save directory
show_LEDs       = false;    % Visualize LED illumination angle in NA space
show_obj        = false;    % Visualize ground truth object
show_raw        = false;    % Visualize the two intensity images (I1,I2)
show_fringe     = false;    % Visualize DAbI fringes (summed spectrum)

% Make save directory
if save_raw
    if ~exist(save_dir,"dir")
        mkdir(save_dir)
    end
end

%% Design ground truth object
% Load sample data
sam = [];
disp('-------Sample Loading------')
tic
for i = 1:5
    filename = sprintf('sam_%d.mat', i);
    data = load(filename);
    sam = cat(3, sam, data.sam_part);
end
sam = double(sam);
toc

disp('------Sample Loaded------')
disp('------Image Generation------')
tic
% Ensure samLayers and zsize are positive integers
assert(samLayers > 0 && zsize >= samLayers, 'Invalid layer sizes');
% Get the center plane of the input stack
center_plane = floor(size(sam, 3) / 2) + 1;
% Compute cropping indices (make sure within bounds)
half_crop = floor(samLayers / 2);
start_crop = max(center_plane - half_crop, 1);
end_crop   = min(center_plane + (samLayers - half_crop - 1), size(sam, 3));
% Crop from original stack
sam_crop = sam(:, :, start_crop:end_crop);
% Initialize the padded output
sam_pad = ones(size(sam,1), size(sam,2), zsize, class(sam)).*n_media;
% Compute where to insert cropped volume
insert_start = floor((zsize - size(sam_crop, 3)) / 2) + 1;
insert_end   = insert_start + size(sam_crop, 3) - 1;
% Assign the cropped stack into the center of the padded volume
sam_pad(:, :, insert_start:insert_end) = (sam_crop-mean2(sam_crop))*dn + n_media;
% Designed object
obj = imresize(sam_pad,[xsize,ysize],'bilinear');
clear sam_pad sam_crop sam

% Visualize ground truth object
if show_obj
    implay(mat2gray(obj))
end

%% Generate Coherent Transfer Function (CTF)  
res_xy = 1/(xsize*pixelSizeXY);           % unit: 1/um. pixel size in x-y frequency domain
res_z = 1/(zsize*pixelSizeZ);             % unit: 1/um. pixel size in z frequency domain
maxCTF = round(NA/(lambda)/res_xy);       % unit: pixel. Calculate the maximal spatial frequency in Fourier domain
[Y,X] = meshgrid(1:ysize,1:xsize);
xc = floor(xsize/2+1);
yc = floor(ysize/2+1);
R_k4img = abs((X-xc) + 1i*(Y-yc));
CTF = R_k4img<= maxCTF;                 


%% Calculate the illumination NA for each LED
% Illumination angles in NA space 
% illumination NA = sqrt( na_calib(:,1).^2 + na_calib(:,2).^2 )
na_calib = zeros(2,2);
na_calib(:,1) = na_illu .* cos(angle2use) ;
na_calib(:,2) = na_illu .* sin(angle2use) ;
kIllu = (ps * xsize * na_calib) / lambda;

% Visualization LED illumination angle in NA space
if show_LEDs
    figure(1002); hold on;
    plot(na_calib(:,1), na_calib(:,2), 'go', 'MarkerFaceColor', 'g', ...
        'MarkerSize', 12);
    title('LED illumination angle in NA space');
    axis square; grid on;
    xlabel('NA_x');
    ylabel('NA_y');
    xlim([-NA,NA])
    ylim([-NA,NA])
end


%% Generate two raw intensity images
% Use multi-sclice beam propagation method to simulate raw intensity
% measurements on the camera  
imlow_defocus = imagingMultiSlice(obj,kIllu,CTF,lambda*1e3,pixelSizeXY,pixelSizeZ,defocus, ...
                            'n',n_media,'gpu',useGPU);

toc
disp('------Image Generation Done------')
% Visualize the two intensity images (I1,I2)
if show_raw
    figure(10031);
    imshow(imlow_defocus(:,:,1),[]);
    title('I1');colormap gray;axis equal;
    figure(10032);
    imshow(imlow_defocus(:,:,2),[]);
    title('I2');colormap gray;axis equal;
end

% Save raw intensity image data
if save_raw
    file_name = sprintf('Cell_3D_z_%d_NA%03.0f_mag%d_DAbI.mat', defocus, NA * 100, mag);
    saveName = fullfile(save_dir, file_name);
    save(saveName,'lambda','na_illu','NA','ps','mag','na_calib',...
              'imlow_defocus','-v7.3'); 
end
% Visualize fringes in the summed spectrum
if show_fringe
    figure(1004);
    imshow(log(1+abs(fftshift(fft2(imlow_defocus(:,:,1)+imlow_defocus(:,:,2))))),[]);
    colormap summer
end

%% Perform DAbI
% load variables to GPU
if useGPU
    imlow_defocus = gpuArray(single(imlow_defocus));
    kIllu = single(kIllu);
    CTF = gpuArray(single(CTF));
else
    imlow_defocus = single(imlow_defocus);
    kIllu = single(kIllu);
    CTF = single(CTF);
end

disp('------Start DAbI------')       
z_defocus = findDefocus_DAbI_3D(imlow_defocus, kIllu, CTF, na_illu, NA, ...
    lambda, mag, ps*mag, sub_pixel_resolve,useGPU,precision, max_iter);
fprintf('DAbI Estimated z: %.4f\n', z_defocus);
fprintf('Ground Truth z:   %.4f\n', defocus);
