clear all; close all; clc
tic


%% Load initial signals and sampling rate
ecg1 = load('signal_13.mat');  
ecg2 = load('signal_50.mat');
fs = 1000;
signal1 = ecg1.signal;
signal2 = ecg2.signal;



%% PART A - מבוסס נגזרות

% PART A - function for 1 best threshold
window = 30;
[Rwaves_A.R1, best_threshold_partA1, ~, ~] = Rwave_detection_D(signal1, fs, window);
[Rwaves_A.R2, best_threshold_partA2, ~, ~] = Rwave_detection_D(signal2, fs, window);

%% Part B – Matched Filter R-peak Detection
template_lengths = [50, 100, 150, 200, 300];

% Signal 1
[Rwaves_B.R1, corr_vec1, chosen1] = Rwave_detection_M(signal1, fs, template_lengths, 1.2);
most_common_L1 = template_lengths(mode(chosen1));
fprintf('Signal 1 - Selected Template Length: %d samples\n', most_common_L1);

% Signal 2
[Rwaves_B.R2, corr_vec2, chosen2] = Rwave_detection_M(signal2, fs, template_lengths, 1.2);
most_common_L2 = template_lengths(mode(chosen2));
fprintf('Signal 2 - Selected Template Length: %d samples\n', most_common_L2);



%% Plotting results


% PART A - plot Best Rwave_detection_D
duration = 20; % Duration is 20 seconds
samples = duration*fs;
time = (1:samples); % time in mili seconds
%[b, a] = butter(4, [5 15] / (fs/2), 'bandpass');
%ecg1_samples = filtfilt(b, a, ecg1(1:samples));
%ecg2_samples = filtfilt(b, a, ecg2(1:samples));
ecg1_samples = signal1(1:samples);
ecg2_samples = signal2(1:samples);

% Get only the R-peaks within the 20-second sample window
Rwaves_R1_ms = Rwaves_A.R1(Rwaves_A.R1 <= samples);
Rwaves_R2_ms = Rwaves_A.R2(Rwaves_A.R2 <= samples);
% Plot for Signal 1 (Signal 13)
figure;
plot(time, ecg1_samples); 
hold on;
plot(time(Rwaves_R1_ms), ecg1_samples(Rwaves_R1_ms), 'ro', 'MarkerFaceColor', 'r');
title(['Signal 13 - Best Threshold = ' num2str(best_threshold_partA1)]);
xlabel('Time [ms]');
ylabel('Signal [mV]');
grid on;
% Plot for Signal 2 (Signal 50)
figure;
plot(time, ecg2_samples); 
hold on;
plot(time(Rwaves_R2_ms), ecg2_samples(Rwaves_R2_ms), 'ro', 'MarkerFaceColor', 'r');
title(['Signal 50 - Best Threshold = ' num2str(best_threshold_partA2)]);
xlabel('Time [ms]');
ylabel('Signal [mV]');
grid on;


% Part B – Matched Filter R-peak Detection
duration_sec = 20;
samples = duration_sec * fs;
t20 = (0:samples-1) / fs;

% Plot raw ECG signals and detected R-peaks for first 20 seconds
figure;
plot(t20, signal1(1:samples));
hold on;
plot(Rwaves_B.R1(Rwaves_B.R1 <= samples) / fs, signal1(Rwaves_B.R1(Rwaves_B.R1 <= samples)), 'ro');
xlabel('Time [s]'); ylabel('Amplitude [mV]');
title('ECG1 – First 20 sec (Raw) with Detected R-peaks');
grid on;

figure;
plot(t20, signal2(1:samples));
hold on;
plot(Rwaves_B.R2(Rwaves_B.R2 <= samples) / fs, signal2(Rwaves_B.R2(Rwaves_B.R2 <= samples)), 'ro');
xlabel('Time [s]'); ylabel('Amplitude [mV]');
title('ECG2 – First 20 sec (Raw) with Detected R-peaks');
grid on;

% Plot correlation vectors for both signals (first 20 seconds) 
figure;
plot(0:samples-1, corr_vec1(1:samples));
xlabel('Sample Index'); ylabel('Correlation Value');
title('ECG1 – Correlation Vector (First 20 Seconds)');
grid on;

figure;
plot(0:samples-1, corr_vec2(1:samples));
xlabel('Sample Index'); ylabel('Correlation Value');
title('ECG2 – Correlation Vector (First 20 Seconds)');
grid on;


%% Part C – ROC Curve 

% PART A

% Load both signals and their ground-truth peaks
load('signal_03.mat');  
signal3 = sig;
load('peaks_03.mat');
true_peaks3 = peaks;
load('signal_04.mat'); 
signal4 = sig;
load('peaks_04.mat');
true_peaks4 = peaks;

% Parameters
fs = 1000;
duration = 20;
theta = 50;
window = 30;
% R-peak detection
[~, ~, R_peaks_all_3, thresholds3] = Rwave_detection_D(signal3, fs, window);
[~, ~, R_peaks_all_4, thresholds4] = Rwave_detection_D(signal4, fs, window);

% Performance evaluation
Pd_3 = zeros(1,5);
Pf_3 = zeros(1,5);
Pd_4 = zeros(1,5);
Pf_4 = zeros(1,5);

for i = 1:5
    [Pd_3(i), Pf_3(i)] = PE(true_peaks3, R_peaks_all_3{i}, theta);
    [Pd_4(i), Pf_4(i)] = PE(true_peaks4, R_peaks_all_4{i}, theta);
end

% שלב 1: ערכי הסף כמו בפועל
thresholds3_real = thresholds3;  % כבר מהפונקציה
thresholds4_real = thresholds4;

% שלב 2: תוויות (sig, thresh)
labels3 = arrayfun(@(t) sprintf('(03, %.3f)', t), thresholds3_real, 'UniformOutput', false);
labels4 = arrayfun(@(t) sprintf('(04, %.3f)', t), thresholds4_real, 'UniformOutput', false);

% שלב 3: מיזוג נתונים
Pf_all = [Pf_3, Pf_4];
Pd_all = [Pd_3, Pd_4];
labels_all = [labels3, labels4];
sig_id = [repmat(3, 1, 5), repmat(4, 1, 5)];  % לזהות מי הסיגנל
marker_style = {'o', 's'};  % עיגול ל־03, ריבוע ל־04
colors = {'b', 'r'};        % כחול ל־03, אדום ל־04

% מיון לפי Pf
[Pf_sorted, sort_idx] = sort(Pf_all);
Pd_sorted = Pd_all(sort_idx);
labels_sorted = labels_all(sort_idx);
sig_sorted = sig_id(sort_idx);

% ציור ROC Curve
figure;
hold on;
% צייר קו מקשר בין הנקודות (לפני ציור הנקודות עצמן)
plot(Pf_sorted, Pd_sorted, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);  % קו אפור עדין
for i = 1:length(Pf_sorted)
    sig = sig_sorted(i);
    color = colors{sig - 2};           % 3 → index 1, 4 → index 2
    marker = marker_style{sig - 2};
    % צייר נקודה
    plot(Pf_sorted(i), Pd_sorted(i), marker, 'MarkerSize', 8,'MarkerSize', 10,'MarkerEdgeColor', color, 'MarkerFaceColor', color);
    % צייר תווית צמודה לנקודה
    text(Pf_sorted(i), Pd_sorted(i), ['  ' labels_sorted{i}], ...
        'FontSize', 12,'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 12, 'Color', color);
end

xlabel('Pf (False Positive Rate)','FontSize', 12);
ylabel('Pd (True Positive Rate)','FontSize', 12);
title('Unified ROC Curve with Threshold Labels','FontSize', 14);

% יצירת ידיות (handles) מלאכותיות ללג'נד
h_line   = plot(NaN, NaN, '-',  'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
h_circle = plot(NaN, NaN, 'o',  'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
h_square = plot(NaN, NaN, 's',  'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
legend([h_line, h_circle, h_square], ...
       {'ROC Curve', 'Signal 03 (Blue Circle)', 'Signal 04 (Red Square)'},'Location', 'best','FontSize', 12);
grid on;


% PART B

% Load signals and peaks
% Signal A = signal_03 (peaks_A), Signal B = signal_04 (peaks_B)
load('signal_03.mat'); signal_A = sig;
load('peaks_03.mat');  peaks_A = peaks;
load('signal_04.mat'); signal_B = sig;
load('peaks_04.mat');  peaks_B = peaks;

% ROC Curve for Matched Filter
fs = 1000;
theta = 50;
template_lengths = [50, 100, 150, 200, 300];
scale = 1.2;

pd_A = zeros(1, length(template_lengths));
pf_A = zeros(1, length(template_lengths));
pd_B = zeros(1, length(template_lengths));
pf_B = zeros(1, length(template_lengths));

% Run detection on both signals with multiple template lengths
for i = 1:length(template_lengths)
    L = template_lengths(i);

    [rpeaks_A, ~, ~] = Rwave_detection_M(signal_A, fs, L, scale);
    [pd_A(i), pf_A(i)] = PE(peaks_A, rpeaks_A, theta);

    [rpeaks_B, ~, ~] = Rwave_detection_M(signal_B, fs, L, scale);
    [pd_B(i), pf_B(i)] = PE(peaks_B, rpeaks_B, theta);
end

% Combine all results
pd_all = [pd_A, pd_B];
pf_all = [pf_A, pf_B];
L_all = [template_lengths, template_lengths];
labels = [repmat({'#03'}, 1, length(template_lengths)), repmat({'#04'}, 1, length(template_lengths))];

% Sort by Pf 
[pf_sorted, sort_idx] = sort(pf_all);
pd_sorted = pd_all(sort_idx);
L_sorted = L_all(sort_idx);
label_sorted = labels(sort_idx);

% Generate colors based on signal
colors = cellfun(@(x) strcmp(x, '#03'), label_sorted); % Assigns 1 to points from signal #03 (A), and 0 to signal #04 (B)
color_map = [0 0 1; 1 0 0]; % blue, red
rgb_colors = color_map(colors+1, :);

% Plot
figure;
hold on; grid on;
xlabel('Pf (False Positive Rate)');
ylabel('Pd (True Positive Rate)');
title(sprintf('ROC Curve – Matched Filter (\\theta = %d )', theta));

% Add artificial start and end points to curve (for smoother ROC shape)
Pf_plot = [0, pf_sorted, 1];
Pd_plot = [0, pd_sorted, 1];

% Connect all points including edges
plot(Pf_plot, Pd_plot, 'k-', 'LineWidth', 1.2);

% Plot each point with appropriate color
for i = 1:length(pd_sorted)
    scatter(pf_sorted(i), pd_sorted(i), 60, rgb_colors(i,:), 'filled');
    offset_y = 0.0007 * (-1)^(i);  % alternate y-offset
    text(pf_sorted(i) + 0.0002, pd_sorted(i) + offset_y, ...
        sprintf('L=%d (%s)', L_sorted(i), label_sorted{i}), ...
        'Color', rgb_colors(i,:), 'FontSize', 8);
end

legend('ROC Curve', 'Location', 'SouthWest');

%% Save Results

% Convert R1 and R2 to row vectors
Rwaves_B.R1 = Rwaves_B.R1(:).';
Rwaves_B.R2 = Rwaves_B.R2(:).';

% Rwaves_B is used to store results from Part B (Matched Filter detection)
save('207298696_207814757.mat', 'Rwaves_B');
save('207298696_207814757.mat', 'Rwaves_A', '-append');
toc

%% function PE

function [Pd, Pf] = PE(true_peaks, detected_peaks, theta)
    TP = 0;
    for i = 1:length(true_peaks)
        if any(abs(detected_peaks - true_peaks(i)) <= theta)
            TP = TP + 1;
        end
    end
    Pd = TP / max(length(true_peaks), 1);  

    FP = 0;
    for i = 1:length(detected_peaks)
        if all(abs(true_peaks - detected_peaks(i)) > theta)
            FP = FP + 1;
        end
    end
    Pf = FP / max(length(detected_peaks), 1);  
end

