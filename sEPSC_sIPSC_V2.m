%% Batch Processing of ABF Files from Two Groups for EPSC/IPSC Analysis
clc;
clear;
close all;

%% Let the user choose folders for Group 1 and Group 2
groupFolders = cell(2,1);
groupFolders{1} = uigetdir(pwd, 'Select folder for Group 1 ABF files');
groupFolders{2} = uigetdir(pwd, 'Select folder for Group 2 ABF files');

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
minEventAmplitude = 10;                % Minimum event amplitude (pA)
minPeakDistanceSamples = round(0.1 * fs);% Minimum interval of 0.1 s in samples

%% Initialize master storage (across all groups)
masterFileSummary_all = table();
masterRawEPSC_all = table();
masterRawIPSC_all = table();

%% Loop over groups (Group 1 and Group 2)
for grp = 1:2
    disp(['Processing Group ', num2str(grp)]);
    % Get list of ABF files in the current group folder
    abfFiles = dir(fullfile(groupFolders{grp}, '*.abf'));
    
    % Loop over each ABF file in this group
    for f = 1:length(abfFiles)
        % Get file name and full path
        abfFileName = fullfile(abfFiles(f).folder, abfFiles(f).name);
        disp(['Processing file: ', abfFileName]);
        
        % Load ABF data (requires abfload)
        [d, si, h] = abfload(abfFileName);
        [numSamples, numChannels, numSweeps] = size(d);
        t = (0:numSamples-1) * dt;  % time vector for each sweep
        
        % Define indices for EPSC and IPSC segments based on time
        idxEPSC = (t >= EPSC_timeRange(1)) & (t <= EPSC_timeRange(2));
        idxIPSC = (t >= IPSC_timeRange(1)) & (t <= IPSC_timeRange(2));
        
        % Design the low-pass Butterworth filter
        [b, a] = butter(order, cutoff/(fs/2), 'low');
        
        % Initialize per-file storage for summary and raw event data
        summaryData = [];
        fileRawEPSC = table();
        fileRawIPSC = table();
        
        %% Loop through each sweep in the current file
        for sweep = 1:numSweeps
            % Extract current trace (channel 1)
            currentTrace = d(:,1,sweep);
            
            % Extract EPSC and IPSC segments and corresponding time vectors
            tEPSC = t(idxEPSC).';  % as column vector
            tIPSC = t(idxIPSC).';
            dataEPSC = currentTrace(idxEPSC);
            dataIPSC = currentTrace(idxIPSC);
            
            % Apply low-pass filtering
            dataEPSC_filt = filtfilt(b, a, dataEPSC);
            dataIPSC_filt = filtfilt(b, a, dataIPSC);
            
            % Baseline correction using polynomial fit
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
            
            %% Append raw event data for this sweep
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
            
        end  % End sweep loop
        
        %% Create per-file summary table (sweep-level)
        summaryTable = array2table(summaryData, ...
            'VariableNames', {'Sweep', 'EPSC_Count', 'EPSC_MeanAmplitude_pA', 'EPSC_MeanISI_s', ...
                               'IPSC_Count', 'IPSC_MeanAmplitude_pA', 'IPSC_MeanISI_s'});
        
        % Write per-file Excel summary with three worksheets:
        % Sheet 1: Summary; Sheet "Raw_EPSC": raw EPSC data; Sheet "Raw_IPSC": raw IPSC data.
        [~, baseName, ~] = fileparts(abfFileName);
        fileSummaryName = fullfile(abfFiles(f).folder, [baseName, '_EPSC_IPSC_summary.xlsx']);
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
        title(['Example Sweep ', num2str(exampleSweep), ' - EPSC (', baseName,')']);
        grid on;
        
        subplot(2,1,2);
        plot(tIPSC_ex, dataIPSC_corrected_ex, 'r'); hold on;
        plot(tIPSC_ex(ipscPeakLocs_ex), dataIPSC_corrected_ex(ipscPeakLocs_ex), 'ko','MarkerFaceColor','y');
        xlabel('Time (s)'); ylabel('IPSC Corrected Current (pA)');
        title(['Example Sweep ', num2str(exampleSweep), ' - IPSC (', baseName,')']);
        grid on;
        
        examplePNGFileName = fullfile(abfFiles(f).folder, [baseName, '_example_sweep.png']);
        saveas(gcf, examplePNGFileName);
        disp(['Saved example PNG: ', examplePNGFileName]);
        
        %% Compute per-file averages (across sweeps) for master summary
        avgEPSC_amp = mean(summaryTable.EPSC_MeanAmplitude_pA, 'omitnan');
        avgEPSC_ISI = mean(summaryTable.EPSC_MeanISI_s, 'omitnan');
        avgIPSC_amp = mean(summaryTable.IPSC_MeanAmplitude_pA, 'omitnan');
        avgIPSC_ISI = mean(summaryTable.IPSC_MeanISI_s, 'omitnan');
        
        fileSummaryRow = table({abfFileName}, grp, avgEPSC_amp, avgEPSC_ISI, avgIPSC_amp, avgIPSC_ISI, ...
            'VariableNames', {'FileName','Group','Avg_EPSC_Amplitude_pA','Avg_EPSC_ISI_s','Avg_IPSC_Amplitude_pA','Avg_IPSC_ISI_s'});
        masterFileSummary_all = [masterFileSummary_all; fileSummaryRow];
        
        % Append raw event data to master raw tables (add a Group column)
        if ~isempty(fileRawEPSC)
            fileRawEPSC.Group = repmat(grp, height(fileRawEPSC), 1);
            masterRawEPSC_all = [masterRawEPSC_all; fileRawEPSC];
        end
        if ~isempty(fileRawIPSC)
            fileRawIPSC.Group = repmat(grp, height(fileRawIPSC), 1);
            masterRawIPSC_all = [masterRawIPSC_all; fileRawIPSC];
        end
    end  % End file loop for current group
end  % End group loop

%% Write master summary Excel file for all ABF files
masterSummaryFileName = 'Summary_EPSC_IPSC.xlsx';
writetable(masterFileSummary_all, masterSummaryFileName, 'Sheet', 1);
disp(['Master summary file generated: ', masterSummaryFileName]);

%% Final cumulative probability plots (Group 1 vs Group 2)

% EPSC Amplitude (convert to positive using abs)
group1_EPSC_amp = abs(masterRawEPSC_all.Amplitude_pA(masterRawEPSC_all.Group==1));
group2_EPSC_amp = abs(masterRawEPSC_all.Amplitude_pA(masterRawEPSC_all.Group==2));

figure('Position', [100 100 500 500]);
% Compute empirical CDF for each group
[f1, x1] = ecdf(group1_EPSC_amp);
[f2, x2] = ecdf(group2_EPSC_amp);
plot(x1, f1, 'b','LineWidth',2); hold on;
plot(x2, f2, 'r','LineWidth',2);
xlabel('EPSC Amplitude (pA)');
ylabel('Cumulative Probability');
title('EPSC Amplitude CDF');
legend('Group 1','Group 2','Location','best');
axis square; grid on;
saveas(gcf, 'Cumulative_Probability_EPSC_Amplitude.png');
disp('Saved cumulative probability plot for EPSC amplitude.');

% IPSC Amplitude
group1_IPSC_amp = masterRawIPSC_all.Amplitude_pA(masterRawIPSC_all.Group==1);
group2_IPSC_amp = masterRawIPSC_all.Amplitude_pA(masterRawIPSC_all.Group==2);

figure('Position', [100 100 500 500]);
[f1, x1] = ecdf(group1_IPSC_amp);
[f2, x2] = ecdf(group2_IPSC_amp);
plot(x1, f1, 'b','LineWidth',2); hold on;
plot(x2, f2, 'r','LineWidth',2);
xlabel('IPSC Amplitude (pA)');
ylabel('Cumulative Probability');
title('IPSC Amplitude CDF');
legend('Group 1','Group 2','Location','best');
axis square; grid on;
saveas(gcf, 'Cumulative_Probability_IPSC_Amplitude.png');
disp('Saved cumulative probability plot for IPSC amplitude.');

% EPSC ISI
group1_EPSC_ISI = masterRawEPSC_all.ISI_s(masterRawEPSC_all.Group==1);
group2_EPSC_ISI = masterRawEPSC_all.ISI_s(masterRawEPSC_all.Group==2);
% Remove NaNs
group1_EPSC_ISI = group1_EPSC_ISI(~isnan(group1_EPSC_ISI));
group2_EPSC_ISI = group2_EPSC_ISI(~isnan(group2_EPSC_ISI));

figure('Position', [100 100 500 500]);
[f1, x1] = ecdf(group1_EPSC_ISI);
[f2, x2] = ecdf(group2_EPSC_ISI);
plot(x1, f1, 'b','LineWidth',2); hold on;
plot(x2, f2, 'r','LineWidth',2);
xlabel('EPSC ISI (s)');
ylabel('Cumulative Probability');
title('EPSC ISI CDF');
legend('Group 1','Group 2','Location','best');
axis square; grid on;
saveas(gcf, 'Cumulative_Probability_EPSC_ISI.png');
disp('Saved cumulative probability plot for EPSC ISI.');

% IPSC ISI
group1_IPSC_ISI = masterRawIPSC_all.ISI_s(masterRawIPSC_all.Group==1);
group2_IPSC_ISI = masterRawIPSC_all.ISI_s(masterRawIPSC_all.Group==2);
group1_IPSC_ISI = group1_IPSC_ISI(~isnan(group1_IPSC_ISI));
group2_IPSC_ISI = group2_IPSC_ISI(~isnan(group2_IPSC_ISI));

figure('Position', [100 100 500 500]);
[f1, x1] = ecdf(group1_IPSC_ISI);
[f2, x2] = ecdf(group2_IPSC_ISI);
plot(x1, f1, 'b','LineWidth',2); hold on;
plot(x2, f2, 'r','LineWidth',2);
xlabel('IPSC ISI (s)');
ylabel('Cumulative Probability');
title('IPSC ISI CDF');
legend('Group 1','Group 2','Location','best');
axis square; grid on;
saveas(gcf, 'Cumulative_Probability_IPSC_ISI.png');
disp('Saved cumulative probability plot for IPSC ISI.');
