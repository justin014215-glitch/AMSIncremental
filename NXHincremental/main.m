%% 北橫磁感率異向性分析與應變增量主程式
clear; close all; clc;

% 輸入資料檔案
filename = 'NXHAMStest.xlsx';

% 設定預設區間劃分參數（如需固定區間平均）
interval_count = 4;           % 區間數量
samples_per_interval = 10;    % 每區樣本數

%% 讀取並計算 Eg（地理座標系統下有限應變橢球）
fprintf('執行原始 Eg 橢球計算...\n');
Eg = compute_Eg_from_AMS(filename);

%% 使用者自訂區間平均 Eg
fprintf('\n=== 自訂區間 Eg 平均分析 ===\n');
[Eg_avg_custom, custom_ranges] = calculate_custom_interval_average_Eg(Eg);

%% 使用固定區間平均進行增量應變分析
fprintf('\n=== 固定區間增量應變分析 ===\n');
[Einc_intervals, U_intervals, R_intervals, ...
 U_eigenvalues_intervals, U_eigenvectors_intervals] = ...
incremental_strain_interval_analysis(Eg, interval_count, samples_per_interval);
