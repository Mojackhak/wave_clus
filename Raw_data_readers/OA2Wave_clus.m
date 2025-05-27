function OA2Wave_clus(matFolder, saveFolder)
% OA2Wave_clus  Convert Alpha-Omega *.mat recordings to Wave_Clus format
%
%   OA2Wave_clus(matFolder, saveFolder)
%
%   • matFolder   – folder that contains the .mat files exported by
%                   NeuroOmega / Alpha Omega (one depth sweep per file).
%   • saveFolder  – destination folder for
%                     ①  *.mat files compatible with Wave_Clus
%                     ②  per-channel type CSV files (time / voltage)
%                     (optional)
%
%   Wave_Clus expects at least two variables inside each *.mat:
%       data   : [samples × channels] int16 or double   (raw or µV)
%       sr     : scalar                                 (Hz)
%
%   Author:  Mojack   Date: 2025.05.27

%% sanity checks
if nargin < 2
    error('Usage: OA2Wave_clus(<matFolder>, <saveFolder>)');
end
if ~isfolder(matFolder)
    error('Input folder "%s" not found.', matFolder);
end
if ~isfolder(saveFolder)
    mkdir(saveFolder);
end

% Alpha Omega channel types you care about – adapt if needed
chanTypes = {'CSPK', 'CLFP'};

save_csv = false;

% make subfolder for each chanTypes
for ct = 1:numel(chanTypes)
    channelType = chanTypes{ct};
    subfolder = fullfile(saveFolder, channelType);
    mkdir(subfolder)
end

% MER range in uV (used to guess gain if absent)
MER_range = [2, 13];

% get list of *.mat files
files = dir(fullfile(matFolder, '*.mat'));
if isempty(files)
    error('No *.mat files found in "%s".', matFolder);
end

% reserve output struct
% merStruct = struct([]);

%% loop over files ---------------------------------------------------------
for f = 1:numel(files)
    matFile = fullfile(files(f).folder, files(f).name);
    fprintf('[%2d/%2d]  %s\n', f, numel(files), files(f).name);

    % ---- load MAT --------------------------------------------------------
    D = load(matFile);                 % gives a struct with many fields
    keys = fieldnames(D);

    % Attempt to infer patient name from parent folder
    [~, parentFolder] = fileparts(files(f).folder);
    patientName = parentFolder;


    % container for this depth
    fileEntry = struct('filename', files(f).name, ...
        'channels',   struct);

    % Iterate over each type (MER, LFP, …)
    for ct = 1:numel(chanTypes)
        channelType = chanTypes{ct};

        % ---- extract the four keys we need ------------------------------
        adKey   = findFirstKey(keys, sprintf('^%s_\\d+$',              channelType));
        bitKey  = findFirstKey(keys, sprintf('^%s_\\d+_BitResolution$', channelType));
        gainKey = findFirstKey(keys, sprintf('^%s_\\d+_Gain$',          channelType));
        srKey   = findFirstKey(keys, sprintf('^%s_\\d+_KHz$',           channelType));

        if isempty(adKey) || isempty(bitKey) || isempty(srKey)
            fprintf('   ↳ %s: missing mandatory keys – skipping.\n', channelType);
            continue
        end

        adArray       = D.(adKey)(1,:);            % row vector
        bitResolution = double(D.(bitKey)(1));     % µV / bit
        if isempty(gainKey)
            fprintf('Gain key missing → Guess gain')
            gain_range = mean(abs(adArray)) * bitResolution./ MER_range;
            if max(gain_range) < 500  % Neuro Omega
                gain = 20;
            elseif min(gain_range) > 500 % NeuroNav
                gain = 1385;
            else
                gain = 1; % default
            end
            fprintf('   ↳ %s: Gain key missing → assumed %d.\n', channelType, gain);
        else
            gain = double(D.(gainKey)(1));
        end
        sr   = double(D.(srKey)(1)) * 1000;        % kHz → Hz

        % ---- scale to µV -------------------------------------------------
        voltage = (double(adArray) .* bitResolution ./ gain)';   % µV
        nSamp   = numel(voltage);
        time    = (0:nSamp-1).' / sr;                          % s

        % ---- write CSV (optional but keeps parity with Python) ----------
        base   = erase(files(f).name, '.mat');
        if save_csv
            csvOut = sprintf('%s_%s_%s.csv', patientName, base, channelType);
            csvOut = fullfile(saveFolder, channelType, csvOut);
            writetable(table(time, voltage), csvOut);
            fprintf('      saved  %s\n', csvOut);
        end

        % ---- populate output struct -------------------------------------
        ch.base   = base;
        ch.side   = base(1);
        ch.tract  = base(2:3);
        ch.depth  = base(5:end-4);
        ch.time           = time;
        ch.data        = voltage;
        ch.bitResolution  = bitResolution;
        ch.gain           = gain;
        ch.sr             = sr;
        ch.unit_sr        = 'Hz';
        ch.unit_time      = 's';
        ch.unit_voltage   = 'uV';
        % fileEntry.channels.(channelType) = ch;

        % ---- Wave_Clus MAT ----------------------------------------------
        wcMat = sprintf('%s_%s_%s_WC.mat', patientName, base, channelType);
        wcMat = fullfile(saveFolder, channelType, wcMat);
        save(wcMat, '-struct', 'ch', '-mat');
        fprintf('      saved  %s\n', wcMat);

    end

    % merStruct = [merStruct; fileEntry]; %#ok<AGROW>
end

fprintf('\nDone.  %d files processed.\n', numel(files));

% -------------------------------------------------------------------------
    function key = findFirstKey(list, regexPattern)
        % helper: return first fieldname that matches regex, or [] if none
        idx = find(~cellfun(@isempty, regexp(list, regexPattern, 'once')), 1);
        if isempty(idx)
            key = [];
        else
            key = list{idx};
        end
    end
% -------------------------------------------------------------------------

% save merStruct
% patientMat = sprintf('%s.mat', patientName);
% patientMat = fullfile(saveFolder, patientMat);
% save(patientMat, 'merStruct', '-mat');
end
