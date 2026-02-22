function [Rwaves_A, best_threshold_A, Rwaves_all_A, all_thresholds_A] = Rwave_detection_D(signal, fs, window)

% About the function:
%
% R-peak detection using 5 different threshold values
%
% Inputs:
%       signal - the input ECG signal
%       fs - sampling frequency
%       win - half window size around QRS_idx to search for local R peak
%
% Outputs:
%       Rwaves_A - the best chosen vector contains detected R-peak indices for the best threshold
%       best_threshold_A - the best threshold for the best chosen vector
%       R_peaks_all - 1x5 cell array, each cell contains detected R-peak indices for a specific threshold
%       thresholds  - list of 5 threshold values used for detection

n = length(signal);

% Step 1: Bandpass filter (5–15 Hz)
[b, a] = butter(4, [5 15] / (fs/2), 'bandpass');
x_filt = filtfilt(b, a, signal);

% Step 2: FS1 derivative
Y1 = zeros(n,1);
for i = 3:n-3
    Y1(i) = abs(x_filt(i+1)-x_filt(i-1));
end
Y2 = zeros(n,1);
for i = 3:n-3
    Y2(i) = abs(x_filt(i+2)-2*x_filt(i)+x_filt(i-2));
end
Y3 = zeros(n,1);
for i = 3:n-3
    Y3(i) = 1.3*Y1(i)+1.1*Y2(i);
end

%%%% Plot for y1, y2, y3
%time = (1:n) / fs * 1000; 
%figure;
%subplot(3,1,1);
%plot(time, Y1, 'b');
%title('Y1 = |x(i+1) - x(i-1)|');
%xlabel('Time [ms]');
%ylabel('Y1');
%grid on;
%subplot(3,1,2);
%plot(time, Y2, 'g');
%title('Y2 = |x(i+2) - 2x(i) + x(i-2)|');
%xlabel('Time [ms]');
%ylabel('Y2');
%grid on;
%subplot(3,1,3);
%plot(time, Y3, 'r');
%title('Y3 = 1.3·Y1 + 1.1·Y2');
%xlabel('Time [ms]');
%ylabel('Y3');
%grid on;

Y = Y3;

% Step 3: Define 5 thresholds as fractions of max(Y)
base_thresh = max(Y);
all_thresholds_A = [0.2, 0.375, 0.55, 0.6, 0.65] * base_thresh;

Rwaves_all_A = cell(1,5);  % will store R-peaks for each threshold
Rwaves_A = []; % best R-peaks

for t = 1:5
    thresh = all_thresholds_A(t);
    QRS_idx = [];

    % Step 4: Detect QRS candidates using current threshold
    for i = 3:n-3
        if Y(i) > thresh && Y(i+1) > thresh && Y(i+2) > thresh
            if Y(i+1)*x_filt(i+1) > 0 && Y(i+2)*x_filt(i+2) > 0
                QRS_idx(end+1) = i;
            end
        end
    end

    % Step 5: For each QRS candidate, find the local max (R peak) in ±win window
    R_peaks = [];
    for idx = QRS_idx
        start_i = max(1, idx - window);
        end_i = min(n, idx + window);
        [~, rel_idx] = max(x_filt(start_i:end_i));
        R_peaks(end+1) = start_i + rel_idx - 1;
    end

    % Store result for this threshold
    Rwaves_all_A{t} = unique(R_peaks); % all ThreshHolds
end

% best threshold chosen
Rwaves_A = Rwaves_all_A{3};
best_threshold_A = all_thresholds_A(3);
end
