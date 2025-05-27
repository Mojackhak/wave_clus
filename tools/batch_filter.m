function batch_filter(matFolder, saveFolder, bandpass)
% batch_filter  Band-pass filter all .mat files in a folder.
%
%   batch_filter(matFolder, saveFolder, bandpass)
%
%   matFolder  – folder that contains .mat files.  Each file must store a
%                struct (or variables) with fields:
%                   data : [samples × channels] or [channels × samples]
%                   sr   : scalar sampling rate (Hz)
%   saveFolder – destination folder for filtered .mat files.
%   bandpass   – 1×2 vector [fLow fHigh]  (Hz), e.g. [300 3000]
%
%   The function designs a zero-phase 4th-order Butterworth band-pass and
%   writes out <originalFileName>_filt.mat in saveFolder.
%
%   Example
%   -------
%     batch_filter('raw_MER', 'filt_MER', [300 3000]);

    arguments
        matFolder  (1,:) char
        saveFolder (1,:) char
        bandpass   (1,2) double {mustBePositive, mustBeIncreasing(bandpass)}
    end

    if ~isfolder(matFolder)
        error('Input folder "%s" not found.', matFolder);
    end
    if ~isfolder(saveFolder)
        mkdir(saveFolder);
    end

    files = dir(fullfile(matFolder, '*.mat'));
    if isempty(files)
        warning('No .mat files found in "%s".', matFolder);
        return
    end

    for k = 1:numel(files)
        inFile  = fullfile(files(k).folder, files(k).name);
        fprintf('[%2d/%2d]  %s  ...', k, numel(files), files(k).name);

        S = load(inFile);                % load into struct S
        % --------------------------------------------------------------
        % Locate data / sr fields (supports arbitrary nested struct)
        % --------------------------------------------------------------
        [data, sr] = findDataSr(S);
        if isempty(data) || isempty(sr)
            fprintf(' skipped (no data/sr).\n');
            continue
        end

        % --------------------------------------------------------------
        % Filter design
        % --------------------------------------------------------------
        nyq = sr/2;
        if bandpass(2) >= nyq
            error('High-cut %.1f Hz ≥ Nyquist (%.1f Hz) in file %s', ...
                  bandpass(2), nyq, files(k).name);
        end
        Wn = bandpass / nyq;                       % normalised
        [b, a] = butter(4, Wn, 'bandpass');

        % Ensure data is (samples × channels)
        transpose_back = false;
        if size(data,1) < size(data,2)           % row-vector or channels×samples
            data = data.';                       % transpose
            transpose_back = true;
        end

        % --------------------------------------------------------------
        % Apply zero-phase filter
        % --------------------------------------------------------------
        data_filt = filtfilt(b, a, double(data));

        if transpose_back
            data_filt = data_filt.';             % restore original orientation
        end

        % --------------------------------------------------------------
        % Save
        % --------------------------------------------------------------
        [~, name] = fileparts(files(k).name);
        outFile = fullfile(saveFolder, sprintf('%s_filt.mat', name));
        data = data_filt; % keep same variable names
        save(outFile, 'data', 'sr', '-mat');

        fprintf(' saved → %s\n', outFile);
    end
end
% ======================================================================
function [data, sr] = findDataSr(S)
% Recursively search struct S for fields named 'data' and 'sr'.
    data = [];
    sr   = [];
    queue = {S};
    while ~isempty(queue) && (isempty(data) || isempty(sr))
        T = queue{1};  queue(1) = [];
        if isstruct(T)
            fn = fieldnames(T);
            for i = 1:numel(fn)
                if strcmp(fn{i}, 'data'), data = T.(fn{i}); end
                if strcmp(fn{i}, 'sr'),   sr   = T.(fn{i}); end
                queue{end+1} = T.(fn{i});               %#ok<AGROW>
            end
        end
    end
end
% ======================================================================
function mustBeIncreasing(vec)
% Custom validation — ensures vec(1) < vec(2).
    if vec(1) >= vec(2)
        eid = 'bandpass:notIncreasing';
        error(eid, 'bandpass must be [fLow fHigh] with fLow < fHigh.');
    end
end
