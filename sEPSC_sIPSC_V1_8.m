%% Batch Processing of ABF Files for EPSC/IPSC Analysis
clc;
clear;
close all;

%% Set up parameters (adjust as needed)
fs = 20000;              % Sampling rate (Hz)
dt = 1/fs;               % Sampling interval (s)
% Time segments (in seconds)
EPSC_timeRange = [1 14];    % EPSC segment from 1 to 14 s
IPSC_timeRange = [16 29];   % IPSC segment from 16 to 29 s

% Filter and smoothing parameters
cutoff = 100;            % Low-pass filter cutoff (Hz)
order = 4;               % Butterworth filter order
smoothingWindow = 50;    % Smoothing window (samples)
polyDegree = 4;          % Degree for baseline correction

% Detection parameters
minEventAmplitude = 15;                % Minimum event amplitude (pA)
minPeakDistanceSamples = round(0.1 * fs);% Minimum interval of 0.1 s in samples

%% Get list of ABF files in the current folder
abfFiles = dir('*.abf');

% Initialize master summary and raw event storage across all files
masterFileSummary = table();
masterRawEPSC = table();
masterRawIPSC = table();

%% Loop over each ABF file
for f = 1:length(abfFiles)
    % Get file name and path
    abfFileName = abfFiles(f).name;
    disp(['Processing file: ', abfFileName]);
    
    % Load ABF data (requires abfload to be available in path)
    [d, si, h] = abfload(abfFileName);
    [numSamples, numChannels, numSweeps] = size(d);
    t = (0:numSamples-1) * dt;  % Time vector for each sweep
    
    % Define indices for EPSC and IPSC segments based on time
    idxEPSC = (t >= EPSC_timeRange(1)) & (t <= EPSC_timeRange(2));
    idxIPSC = (t >= IPSC_timeRange(1)) & (t <= IPSC_timeRange(2));
    
    % Design the low-pass Butterworth filter
    [b, a] = butter(order, cutoff/(fs/2), 'low');
    
    % Initialize per-file storage for summary (per sweep) and raw events
    summaryData = [];
    fileRawEPSC = table();
    fileRawIPSC = table();
    
    %% Loop through each sweep in the current file
    for sweep = 1:numSweeps
        % Extract current trace (channel 1)
        currentTrace = d(:,1,sweep);
        
        % EPSC and IPSC segments and time vectors
        tEPSC = t(idxEPSC).';  % ensure column vector
        tIPSC = t(idxIPSC).';
        dataEPSC = currentTrace(idxEPSC);
        dataIPSC = currentTrace(idxIPSC);
        
        % Apply low-pass filtering
        dataEPSC_filt = filtfilt(b, a, dataEPSC);
        dataIPSC_filt = filtfilt(b, a, dataIPSC);
        
        % Baseline correction via polynomial fit
        p_EPSC = polyfit(tEPSC, dataEPSC_filt, polyDegree);
        baseline_EPSC = polyval(p_EPSC, tEPSC);
        dataEPSC_corrected = dataEPSC_filt - baseline_EPSC;
        
        p_IPSC = polyfit(tIPSC, dataIPSC_filt, polyDegree);
        baseline_IPSC = polyval(p_IPSC, tIPSC);
        dataIPSC_corrected = dataIPSC_filt - baseline_IPSC;
        
        % Smooth the corrected signals
        dataEPSC_smooth = smoothdata(dataEPSC_corrected, 'movmean', smoothingWindow);
        dataIPSC_smooth = smoothdata(dataIPSC_corrected, 'movmean', smoothingWindow);
        
        %% Detect EPSC events (negative deflections)
        % Invert signal so that downward deflections become positive peaks
        [epscPeakVals, epscPeakLocs] = findpeaks(-dataEPSC_smooth, ...
            'MinPeakHeight', minEventAmplitude, 'MinPeakDistance', minPeakDistanceSamples);
        % Convert inverted values back to negative amplitudes
        EPSC_amplitudes = -epscPeakVals;
        EPSC_eventTimes = tEPSC(epscPeakLocs);
        if ~isempty(EPSC_eventTimes)
            epscISI = [NaN; diff(EPSC_eventTimes)];
        else
            epscISI = [];
        end
        
        %% Detect IPSC events (positive deflections)
        [ipscPeakVals, ipscPeakLocs] = findpeaks(dataIPSC_smooth, ...
            'MinPeakHeight', minEventAmplitude, 'MinPeakDistance', minPeakDistanceSamples);
        IPSC_amplitudes = ipscPeakVals;
        IPSC_eventTimes = tIPSC(ipscPeakLocs);
        if ~isempty(IPSC_eventTimes)
            ipscISI = [NaN; diff(IPSC_eventTimes)];
        else
            ipscISI = [];
        end
        
        %% Compute summary metrics for this sweep
        EPSC_count = numel(EPSC_amplitudes);
        if EPSC_count > 0
            EPSC_meanAmp = mean(EPSC_amplitudes);
        else
            EPSC_meanAmp = NaN;
        end
        if numel(EPSC_eventTimes) > 1
            EPSC_meanISI = mean(diff(EPSC_eventTimes));
        else
            EPSC_meanISI = NaN;
        end
        
        IPSC_count = numel(IPSC_amplitudes);
        if IPSC_count > 0
            IPSC_meanAmp = mean(IPSC_amplitudes);
        else
            IPSC_meanAmp = NaN;
        end
        if numel(IPSC_eventTimes) > 1
            IPSC_meanISI = mean(diff(IPSC_eventTimes));
        else
            IPSC_meanISI = NaN;
        end
        
        summaryData = [summaryData; sweep, EPSC_count, EPSC_meanAmp, EPSC_meanISI, IPSC_count, IPSC_meanAmp, IPSC_meanISI];
        
        %% Append raw event data with sweep info
        if ~isempty(EPSC_amplitudes)
            nEPSC = numel(EPSC_amplitudes);
            tempEPSC = table(repmat(sweep, nEPSC, 1), repmat("EPSC", nEPSC, 1), (1:nEPSC)', ...
                EPSC_amplitudes, epscISI, 'VariableNames', {'Sweep','EventType','EventNumber','Amplitude_pA','ISI_s'});
            fileRawEPSC = [fileRawEPSC; tempEPSC];
        end
        
        if ~isempty(IPSC_amplitudes)
            nIPSC = numel(IPSC_amplitudes);
            tempIPSC = table(repmat(sweep, nIPSC, 1), repmat("IPSC", nIPSC, 1), (1:nIPSC)', ...
                IPSC_amplitudes, ipscISI, 'VariableNames', {'Sweep','EventType','EventNumber','Amplitude_pA','ISI_s'});
            fileRawIPSC = [fileRawIPSC; tempIPSC];
        end
        
    end  % end sweep loop
    
    %% Create per-file summary table (for sweeps)
    summaryTable = array2table(summaryData, ...
        'VariableNames', {'Sweep', 'EPSC_Count', 'EPSC_MeanAmplitude_pA', 'EPSC_MeanISI_s', ...
                           'IPSC_Count', 'IPSC_MeanAmplitude_pA', 'IPSC_MeanISI_s'});
    
    % Write per-file Excel summary with three worksheets:
    %   Sheet 1: Summary, Sheet "Raw_EPSC": raw EPSC data, Sheet "Raw_IPSC": raw IPSC data.
    [filepath, name, ext] = fileparts(abfFileName);
    fileSummaryName = fullfile(filepath, [name, '_EPSC_IPSC_summary.xlsx']);
    writetable(summaryTable, fileSummaryName, 'Sheet', 1);
    writetable(fileRawEPSC, fileSummaryName, 'Sheet', 'Raw_EPSC');
    writetable(fileRawIPSC, fileSummaryName, 'Sheet', 'Raw_IPSC');
    disp(['Generated summary Excel for file: ', fileSummaryName]);
    
    %% Save an example PNG for the file (using sweep 1)
    exampleSweep = 1;
    currentTrace_ex = d(:,1,exampleSweep);
    tEPSC_ex = t(idxEPSC).';
    tIPSC_ex = t(idxIPSC).';
    dataEPSC_ex = currentTrace_ex(idxEPSC);
    dataIPSC_ex = currentTrace_ex(idxIPSC);
    
    dataEPSC_filt_ex = filtfilt(b, a, dataEPSC_ex);
    dataIPSC_filt_ex = filtfilt(b, a, dataIPSC_ex);
    
    p_EPSC_ex = polyfit(tEPSC_ex, dataEPSC_filt_ex, polyDegree);
    baseline_EPSC_ex = polyval(p_EPSC_ex, tEPSC_ex);
    dataEPSC_corrected_ex = dataEPSC_filt_ex - baseline_EPSC_ex;
    
    p_IPSC_ex = polyfit(tIPSC_ex, dataIPSC_filt_ex, polyDegree);
    baseline_IPSC_ex = polyval(p_IPSC_ex, tIPSC_ex);
    dataIPSC_corrected_ex = dataIPSC_filt_ex - baseline_IPSC_ex;
    
    dataEPSC_smooth_ex = smoothdata(dataEPSC_corrected_ex, 'movmean', smoothingWindow);
    dataIPSC_smooth_ex = smoothdata(dataIPSC_corrected_ex, 'movmean', smoothingWindow);
    
    [epscPeakVals_ex, epscPeakLocs_ex] = findpeaks(-dataEPSC_smooth_ex, 'MinPeakHeight', minEventAmplitude, 'MinPeakDistance', minPeakDistanceSamples);
    [ipscPeakVals_ex, ipscPeakLocs_ex] = findpeaks(dataIPSC_smooth_ex, 'MinPeakHeight', minEventAmplitude, 'MinPeakDistance', minPeakDistanceSamples);
    
    figure;
    subplot(2,1,1);
    plot(tEPSC_ex, dataEPSC_corrected_ex, 'b'); hold on;
    plot(tEPSC_ex(epscPeakLocs_ex), dataEPSC_corrected_ex(epscPeakLocs_ex), 'ko','MarkerFaceColor','g');
    xlabel('Time (s)'); ylabel('EPSC Corrected Current (pA)');
    title(['Example Sweep ', num2str(exampleSweep), ' - EPSC (', abfFileName,')']);
    grid on;
    
    subplot(2,1,2);
    plot(tIPSC_ex, dataIPSC_corrected_ex, 'r'); hold on;
    plot(tIPSC_ex(ipscPeakLocs_ex), dataIPSC_corrected_ex(ipscPeakLocs_ex), 'ko','MarkerFaceColor','y');
    xlabel('Time (s)'); ylabel('IPSC Corrected Current (pA)');
    title(['Example Sweep ', num2str(exampleSweep), ' - IPSC (', abfFileName,')']);
    grid on;
    
    examplePNGFileName = fullfile(filepath, [name, '_example_sweep.png']);
    saveas(gcf, examplePNGFileName);
    disp(['Saved example PNG: ', examplePNGFileName]);
    
    %% Compute per-file averages (across sweeps) for master summary
    avgEPSC_amp = mean(summaryTable.EPSC_MeanAmplitude_pA, 'omitnan');
    avgEPSC_ISI = mean(summaryTable.EPSC_MeanISI_s, 'omitnan');
    avgIPSC_amp = mean(summaryTable.IPSC_MeanAmplitude_pA, 'omitnan');
    avgIPSC_ISI = mean(summaryTable.IPSC_MeanISI_s, 'omitnan');
    
    fileSummaryRow = table({abfFileName}, avgEPSC_amp, avgEPSC_ISI, avgIPSC_amp, avgIPSC_ISI, ...
        'VariableNames', {'FileName','Avg_EPSC_Amplitude_pA','Avg_EPSC_ISI_s','Avg_IPSC_Amplitude_pA','Avg_IPSC_ISI_s'});
    masterFileSummary = [masterFileSummary; fileSummaryRow];
    
    % Append raw event data to master raw tables (adding a FileName column)
    if ~isempty(fileRawEPSC)
        fileRawEPSC.FileName = repmat({abfFileName}, height(fileRawEPSC), 1);
        masterRawEPSC = [masterRawEPSC; fileRawEPSC];
    end
    if ~isempty(fileRawIPSC)
        fileRawIPSC.FileName = repmat({abfFileName}, height(fileRawIPSC), 1);
        masterRawIPSC = [masterRawIPSC; fileRawIPSC];
    end
end  % end file loop

%% Write master summary Excel file for all ABF files
masterSummaryFileName = 'Summary_EPSC_IPSC.xlsx';
writetable(masterFileSummary, masterSummaryFileName, 'Sheet', 1);
disp(['Master summary file generated: ', masterSummaryFileName]);

%% Generate cumulative probability plots (separate figures, square shape)

% --- Cumulative probability for EPSC Amplitude ---
% Convert negative EPSC amplitudes to positive using abs()
figure('Position', [100 100 500 500]);
[fEPSC_amp, xEPSC_amp] = ecdf(abs(masterRawEPSC.Amplitude_pA));
plot(xEPSC_amp, fEPSC_amp, 'b','LineWidth',2);
xlabel('EPSC Amplitude (pA)');
ylabel('Cumulative Probability');
title('Cumulative Probability for EPSC Amplitude');
axis square;
grid on;
saveas(gcf, 'Cumulative_Probability_EPSC_Amplitude.png');
disp('Saved EPSC amplitude cumulative probability plot.');

% --- Cumulative probability for IPSC Amplitude ---
figure('Position', [100 100 500 500]);
[fIPSC_amp, xIPSC_amp] = ecdf(masterRawIPSC.Amplitude_pA);
plot(xIPSC_amp, fIPSC_amp, 'r','LineWidth',2);
xlabel('IPSC Amplitude (pA)');
ylabel('Cumulative Probability');
title('Cumulative Probability for IPSC Amplitude');
axis square;
grid on;
saveas(gcf, 'Cumulative_Probability_IPSC_Amplitude.png');
disp('Saved IPSC amplitude cumulative probability plot.');

% --- Cumulative probability for EPSC ISI ---
% Remove NaN values from EPSC ISI
EPSC_ISI_all = masterRawEPSC.ISI_s(~isnan(masterRawEPSC.ISI_s));
figure('Position', [100 100 500 500]);
[fEPSC_ISI, xEPSC_ISI] = ecdf(EPSC_ISI_all);
plot(xEPSC_ISI, fEPSC_ISI, 'b','LineWidth',2);
xlabel('EPSC ISI (s)');
ylabel('Cumulative Probability');
title('Cumulative Probability for EPSC ISI');
axis square;
grid on;
saveas(gcf, 'Cumulative_Probability_EPSC_ISI.png');
disp('Saved EPSC ISI cumulative probability plot.');

% --- Cumulative probability for IPSC ISI ---
IPSC_ISI_all = masterRawIPSC.ISI_s(~isnan(masterRawIPSC.ISI_s));
figure('Position', [100 100 500 500]);
[fIPSC_ISI, xIPSC_ISI] = ecdf(IPSC_ISI_all);
plot(xIPSC_ISI, fIPSC_ISI, 'r','LineWidth',2);
xlabel('IPSC ISI (s)');
ylabel('Cumulative Probability');
title('Cumulative Probability for IPSC ISI');
axis square;
grid on;
saveas(gcf, 'Cumulative_Probability_IPSC_ISI.png');
disp('Saved IPSC ISI cumulative probability plot.');
