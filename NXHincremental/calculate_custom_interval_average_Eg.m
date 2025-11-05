function [Eg_avg_custom, custom_ranges] = calculate_custom_interval_average_Eg(Eg)
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
        Eg_sum = zeros(3,3);
        for j = start_idx:end_idx
            Eg_sum = Eg_sum + Eg(:,:,j);
        end
        Eg_avg_custom(:,:,i) = Eg_sum / (end_idx - start_idx + 1);
    end

    fprintf('\n=== Eg 區間平均結果 ===\n');
    for i = 1:num_intervals
        fprintf('區間 %d：樣本 %d - %d\n', i, custom_ranges(i,1), custom_ranges(i,2));
        disp(Eg_avg_custom(:,:,i));
    end
end
