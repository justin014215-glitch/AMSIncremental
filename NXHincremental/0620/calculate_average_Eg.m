function [Eg_avg, sample_ranges] = calculate_average_Eg(Eg, interval_count, samples_per_interval)
% 自動等量分區平均 Eg
% 輸入：Eg 張量、區間數量、每區樣本數
% 輸出：平均 Eg（每區一個）與起訖樣本索引

    total_samples = size(Eg, 3);
    Eg_avg = zeros(3, 3, interval_count);
    sample_ranges = zeros(interval_count, 2);

    for i = 1:interval_count
        start_idx = (i-1)*samples_per_interval + 1;
        end_idx = min(i*samples_per_interval, total_samples);
        Eg_avg(:,:,i) = mean(Eg(:,:,start_idx:end_idx), 3);
        sample_ranges(i,:) = [start_idx, end_idx];
    end
end
