%% 北橫磁感率異向性應變增量主程式

clear; clc; close all;

%% 參數設定與檔案讀取
filename = 'NXHAMS.xlsx';   % 原始 AMS 資料
interval_count = 4;             % 固定區間數
samples_per_interval = 10;      % 每區樣本數

% 讀取與轉換 Eg（地理座標系有限應變張量）
fprintf('==> 開始讀取並計算 Eg...\n');
Eg = compute_Eg_from_AMS(filename);

%% 使用者自訂區間平均 Eg
fprintf('\n==> 進行自訂區間 Eg 平均分析...\n');
[Eg_avg_custom, custom_ranges] = calculate_custom_interval_average_Eg(Eg);

%% 固定區間增量應變分析
fprintf('\n==> 進行固定區間增量應變分析...\n');
[Einc, U, R, eigvals_U, eigvecs_U, Eg_avg, sample_ranges] = ...
    incremental_strain_interval_analysis(Eg, interval_count, samples_per_interval);

%% 輸出與視覺化
fprintf('\n==> 輸出分析結果與圖形...\n');

% 顯示區間劃分資訊
fprintf('\n--- 區間劃分資訊 ---\n');
for i = 1:interval_count
    fprintf('區間 %d: 樣本 %d 到 %d\n', i, sample_ranges(i,1), sample_ranges(i,2));
end

% 顯示 Eg 平均結果
fprintf('\n--- 區間平均 Eg ---\n');
for i = 1:interval_count
    fprintf('區間 %d:\n'); disp(Eg_avg(:,:,i));
end

% 顯示增量應變與張量結果
for i = 1:interval_count-1
    fprintf('\n--- 區間 %d → %d 應變增量分析 ---\n', i, i+1);
    fprintf('Einc:\n'); disp(Einc(:,:,i));
    fprintf('U (stretch tensor):\n'); disp(U(:,:,i));
    fprintf('R (rotation tensor):\n'); disp(R(:,:,i));
    fprintf('主應變值: [%.4f, %.4f, %.4f]\n', eigvals_U(1,i), eigvals_U(2,i), eigvals_U(3,i));
    fprintf('主方向:\n'); disp(eigvecs_U(:,:,i));
end

%% 圖形化主值變化
figure;
subplot(1,1,1);
bar(1:interval_count-1, eigvals_U', 'grouped');
xlabel('區間對應段');
ylabel('主應變值');
legend({'U1','U2','U3'});
title('各區間之形變率主值');
grid on;


%% 匯出分析結果

fprintf('\n==> 儲存結果中...\n');

% 匯出主值為表格
T = table(sample_ranges(1:end-1,1), sample_ranges(1:end-1,2), ...
          sample_ranges(2:end,1), sample_ranges(2:end,2), ...
          eigvals_U(1,:)', eigvals_U(2,:)', eigvals_U(3,:)', ...
          'VariableNames', {'區間起始','區間結束','下一區起始','下一區結束','U1','U2','U3'});

[~, name, ~] = fileparts(filename);
writetable(T, [name '_incremental_result.xlsx']);

% 儲存變數至 .mat
save([name '_analysis_workspace.mat'], 'Eg', 'Eg_avg', 'Einc', 'U', 'R', ...
     'eigvals_U', 'eigvecs_U', 'sample_ranges', 'custom_ranges');

fprintf('✅ 全部完成！\n');
