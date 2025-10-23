zipFile = 'UKB_37k_connectome_data.zip';

% Unzip into current folder
unzip(zipFile, pwd);

weights = ["MD", "FA", "SC", "ISOVF", "ICVF"];

for w = weights
    weight = char(w);  % convert to char for compatibility with dir/fullfile
    
    % Build path to all matching CSVs
    connectome_pattern = fullfile(pwd, 'UKB_37k_connectome_data', 'matrices', weight, '*Raw*.csv');
    files = dir(connectome_pattern);
    L = numel(files);
    
    if L == 0
        warning('No files found for weight: %s', weight);
        continue;
    end
    
    % Initialize struct
    Connectome = struct();
    Connectome.(weight) = zeros(85, 85, L);
    
    % Initialize ID only once
    if ~exist('ID', 'var')
        ID = strings(L, 1);
    end
    
    filefolder = files(1).folder;
    
    for i = 1:L
        if ID(i) == ""
            thisfile = files(i).name;
            ID(i) = string(thisfile(1:7));
        else
            thisfile = ID(i) + "_Raw_" + weight + ".csv";
        end
        
        single_connectome_file = fullfile(filefolder, thisfile);
        Connectome.(weight)(:, :, i) = readmatrix(single_connectome_file);
    end
    
    % Save per-weight file
    writefilename = "UKB_37k_connectome_Raw_" + weight + ".mat";
    save(writefilename, "Connectome");
end

% Save the subject IDs once
save("UKB_37k_connectome_ID.mat", "ID");