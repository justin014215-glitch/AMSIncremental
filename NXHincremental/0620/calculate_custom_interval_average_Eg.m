function [Eg_avg_custom, custom_ranges] = calculate_custom_interval_average_Eg(Eg)
% 使用者手動輸入每個區間的起訖樣本，計算 Eg 的平均
% 輸入：Eg 張量（3×3×N）
% 輸出：平均 Eg（每區一個）、起訖索引紀錄

    total_samples = size(Eg, 3);
    num_intervals = input('請輸入自訂區間數量： ');

    Eg_avg_custom = zeros(3, 3, num_intervals);
    custom_ranges = zeros(num_intervals, 2);

    for i = 1:num_intervals
        fprintf('\n--- 區間 %d ---\n', i);
        start_idx = input('起始樣本編號：');
        end_idx = input('結束樣本編號：');

        if start_idx < 1 || end_idx > total_samples || start_idx > end_idx
            error('樣本範圍錯誤');
        end

        custom_ranges(i,:) = [start_idx, end_idx];
        Eg_avg_custom(:,:,i) = mean(Eg(:,:,start_idx:end_idx), 3);
    end
end

