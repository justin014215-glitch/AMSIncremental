function [Eg_avg, sample_ranges] = calculate_average_Eg(Eg, interval_count, samples_per_interval)
% CALCULATE_AVERAGE_EG 計算地理座標系下有限應變橢球Eg的區間平均值
%
% 輸入參數:
%   Eg - 3x3xN 陣列，包含N個樣本的有限應變橢球
%   interval_count - 要劃分的區間數量
%   samples_per_interval - 每個區間內的樣本數量
%
% 輸出參數:
%   Eg_avg - 3x3xM 陣列，包含M個區間的平均有限應變橢球
%   sample_ranges - Mx2 陣列，每行包含一個區間的起始和結束樣本索引

    % 獲取樣本總數
    total_samples = size(Eg, 3);
    
    % 驗證參數合理性
    if interval_count * samples_per_interval > total_samples
        error('區間數量和每個區間的樣本數超過了總樣本數。');
    end
    
    % 初始化平均Eg矩陣和樣本範圍
    Eg_avg = zeros(3, 3, interval_count);
    sample_ranges = zeros(interval_count, 2);
    
    % 對每個區間計算平均值
    for i = 1:interval_count
        % 計算當前區間的樣本範圍
        start_idx = (i-1) * samples_per_interval + 1;
        end_idx = min(i * samples_per_interval, total_samples);
        sample_ranges(i, :) = [start_idx, end_idx];
        
        % 初始化當前區間的累加矩陣
        Eg_sum = zeros(3, 3);
        
        % 累加當前區間內的所有Eg矩陣
        for j = start_idx:end_idx
            Eg_sum = Eg_sum + Eg(:,:,j);
        end
        
        % 計算平均值
        Eg_avg(:,:,i) = Eg_sum / (end_idx - start_idx + 1);
    end
    
    % 顯示區間劃分結果
    fprintf('區間劃分結果：\n');
    for i = 1:interval_count
        fprintf('區間 %d: 樣本 %d 到 %d\n', i, sample_ranges(i, 1), sample_ranges(i, 2));
    end
end


