

%% Set up parameters
folder = 'E:\MCIAC\Data\DualLED_Range_247_20x\Data_m'; % Calibration Data to be load
folder_name = dir(folder);
folder_name = folder_name(3:end);

% ---------- copy settings ----------
dst_root = 'E:\MCIAC\DAbI Code Formal\DAbI_experiments_2D\Data';

% ---------- your loop ----------
for sample_idx = [1,3,4,9,11,12,13,14,15,20]
    sample_name = folder_name(sample_idx).name;
    filepath = [folder filesep sample_name filesep];


    zlist = [0:10:30, 50:50:1200];
    for z_pos = zlist
        curr_z = num2str(floor((z_pos)*1000), "%06.0f");
        filename = [filepath sample_name '_z_' curr_z '_nm_DAbI.mat']; 

        dst_folder = fullfile(dst_root, sample_name);
        if ~exist(dst_folder)
            mkdir(dst_folder);
        end

        dst_file = fullfile(dst_folder, [sample_name '_z_' curr_z '_nm_DAbI.mat']);
        if exist(filename, 'file')
            copyfile(filename, dst_file, 'f'); 
        else
            warning('File not found, skip: %s', filename);
        end
    end
end



