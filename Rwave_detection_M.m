function [R_locs, corr_vec, chosen_template_indices] = Rwave_detection_M(ecg, fs, template_lengths, scale)

    if nargin < 4
        scale = 1;
    end

    %% Step 1: Filter the signal using a 4th-order Butterworth bandpass filter (5-15 Hz)
    ecg = ecg(:);  % Ensure column vector
    Wn = [5 15] / (fs/2);
    [b, a] = butter(4, Wn, 'bandpass');
    ecg_filtered = filtfilt(b, a, ecg);  % Zero-phase filtering

    %% Step 2: Normalize the signal to the range [-1, 1]
    ecg_norm = -1 + 2 * (ecg_filtered - min(ecg_filtered)) / (max(ecg_filtered) - min(ecg_filtered));

    %% Step 3: Select template center from strongest peak between 1s and 3s
    search_range = round(1*fs):round(3*fs);
    [~, rel_idx] = max(ecg_filtered(search_range));
    center_idx = rel_idx + search_range(1) - 1;

    %% Step 4: Extract templates of different lengths centered on the chosen peak
    templates = {};
    valid_lengths = [];  % Keep track of which lengths are valid
    for i = 1:length(template_lengths)
        L = template_lengths(i);
        start_idx = center_idx - floor((L-1)/2);
        end_idx = center_idx + floor(L/2);

        if start_idx >= 1 && end_idx <= length(ecg_norm)
            temp_norm = ecg_norm(start_idx:end_idx);
            if std(temp_norm) >= 0.01
                templates{end+1} = temp_norm;
                valid_lengths(end+1) = L;
            end
        end
    end
    % Check if any valid templates were found
    if isempty(valid_lengths)
        error('No valid templates found. Check signal quality or template lengths.');
    end

    %% Step 5: Compute  correlation between signal and each template
    num_templates = length(templates);
    corr_mat = zeros(length(ecg_norm), num_templates);
    for i = 1:num_templates
        temp = templates{i};
        corr_raw = sliding_corr(ecg_norm, temp);
        corr_mat(:, i) = movmean(corr_raw, round(0.01 * fs));  % smooth with 10ms window
    end
   
    [~, best_template_idx] = max(sum(corr_mat, 1));
    corr_vec = corr_mat(:, best_template_idx);
    best_template_idx_vec = best_template_idx * ones(size(corr_vec));

    %% Step 6: Find candidate peaks using adaptive thresholding in 10s windows
    win_size = 10 * fs;
    raw_peaks = [];
    min_corr_threshold = 0.1;

    for start_idx = 1:win_size:length(corr_vec)
        end_idx = min(start_idx + win_size - 1, length(corr_vec));
        local_corr = corr_vec(start_idx:end_idx);
        if std(local_corr) < 0.01
            continue;
        end
        local_thresh = max(mean(local_corr) + scale * std(local_corr), min_corr_threshold);
        local_peaks = find(local_corr >= local_thresh) + start_idx - 1;
        raw_peaks = [raw_peaks; local_peaks];
    end

    %% Step 7: Keep one strong peak per 250ms window and filter by amplitude
    min_dist = round(0.25 * fs);
    temp_peaks = [];
    chosen_template_indices = [];
    amp_thresh_global = prctile(ecg_filtered, 90);
    i = 1;
    while i <= length(raw_peaks)
        curr = raw_peaks(i);
        window_mask = raw_peaks >= curr & raw_peaks <= curr + min_dist;
        window_peaks = raw_peaks(window_mask);
        [~, imax] = max(corr_vec(window_peaks));
        peak_idx = window_peaks(imax);
        if corr_vec(peak_idx) > 0.25 && ecg_filtered(peak_idx) >= amp_thresh_global
            temp_peaks(end+1,1) = peak_idx;
            chosen_template_indices(end+1,1) = best_template_idx_vec(peak_idx);
        end
        i = find(raw_peaks > curr + min_dist, 1);
        if isempty(i), break; end
    end

    %% Step 8: Remove low amplitude peaks (below 70% of median)
    amps = ecg_filtered(temp_peaks);
    amp_thresh = median(amps) * 0.7;
    valid_idx = amps >= amp_thresh;
    temp_peaks = temp_peaks(valid_idx);
    chosen_template_indices = chosen_template_indices(valid_idx);

    %% Step 9: Refine each peak using both filtered and raw signal
    window_size = round(0.1 * fs);  % 100 ms
    R_locs = zeros(size(temp_peaks));
    for j = 1:length(temp_peaks)
        idx = temp_peaks(j);
        s = max(1, idx - window_size);
        e = min(length(ecg), idx + window_size);

        % Find local max in both filtered and raw ECG
        [~, i1] = max(ecg_filtered(s:e));
        [~, i2] = max(ecg(s:e));

        % Choose the one with higher absolute amplitude
        if abs(ecg_filtered(s + i1 - 1)) >= abs(ecg(s + i2 - 1))
            R_locs(j) = s + i1 - 1;
        else
            R_locs(j) = s + i2 - 1;
        end
    end

    %% Step 10: Final check to ensure peaks are at least 250ms apart
    min_final_dist = round(0.25 * fs);
    R_locs = sort(R_locs);
    final_peaks = [];
    last_accepted = -Inf;
    for i = 1:length(R_locs)
        if R_locs(i) - last_accepted >= min_final_dist
            final_peaks(end+1,1) = R_locs(i);
            last_accepted = R_locs(i);
        end
    end
    R_locs = final_peaks;
end

%% function sliding_corr 
% Compute sliding correlation between signal and template
% Based on the correlation formula taught in the tutorial 

function c = sliding_corr(x, y)
    N = length(x);
    L = length(y);
    c = zeros(N, 1);
    y_mean = mean(y);
    y_std = std(y);
    for i = 1:(N - L + 1)
        x_seg = x(i:i+L-1);
        x_mean = mean(x_seg);
        x_std = std(x_seg);
        if x_std == 0 || y_std == 0
            c(i + floor(L/2)) = 0;
        else
            numerator = sum((x_seg - x_mean) .* (y - y_mean));
            c(i + floor(L/2)) = numerator / ((L - 1) * x_std * y_std);
        end
    end
end
