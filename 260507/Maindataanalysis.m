%% AMS_MainAnalysis.m
% 主程式：依 Number 分組，對每組磁化率資料進行：
%   1. 計算各樣本的 Kg (地理座標系張量)
%   2. 標準化張量 (standardizeTensor)
%   3. 計算平均張量 (calculateMeanTensor)
%   4. 特徵值分析 (performEigenAnalysis)
%   5. Jelinek (1978) 統計誤差
%   6. Bootstrap 統計誤差
%
% 依賴：AMSFormulas.m (需與本檔案在同一目錄或 MATLAB Path 中)
%
% 輸入資料欄位 (Excel)：
%   Number, K1, K2, K3, dK1geo(trend), iK1geo(plunge),
%                        dK2geo, iK2geo, dK3geo, iK3geo
% =========================================================

clearvars; clc;

%% -------- 0. 讀取資料 ----------------------------------------
[file, path] = uigetfile('*.xlsx', '選擇原始資料 Excel 檔案');
if isequal(file, 0)
    error('未選擇檔案，程式結束。');
end
rawData = readtable(fullfile(path, file), 'VariableNamingRule', 'preserve');

fprintf('已讀取 %d 筆資料\n', height(rawData));

% 取出欄位（用 {} 存取，避免欄位名稱特殊字元問題）
Numbers  = rawData{:, 'Number'};
Names    = rawData{:, 'Name'};    % 樣本編號，如 2023020601A
K1_all   = rawData{:, 'K1'};
K2_all   = rawData{:, 'K2'};
K3_all   = rawData{:, 'K3'};
dK1_all  = rawData{:, 'dK1geo'}; % K1 Trend (度)
iK1_all  = rawData{:, 'iK1geo'}; % K1 Plunge (度)
dK2_all  = rawData{:, 'dK2geo'};
iK2_all  = rawData{:, 'iK2geo'};
dK3_all  = rawData{:, 'dK3geo'};
iK3_all  = rawData{:, 'iK3geo'};

%% -------- 1. 依 Number 分組處理 ------------------------------
uniqueNums = unique(Numbers, 'stable');
nGroups    = numel(uniqueNums);

% 預先分配結果表格欄位
results = struct();
results.Number      = zeros(nGroups, 1);
results.N_samples   = zeros(nGroups, 1);
results.SampleNames = cell(nGroups, 1);   % 各組樣本編號清單（逗號分隔字串）

% 平均張量特徵值與方向
results.mean_K1    = zeros(nGroups, 1);
results.mean_K2    = zeros(nGroups, 1);
results.mean_K3    = zeros(nGroups, 1);
results.mean_K1_trend  = zeros(nGroups, 1);
results.mean_K1_plunge = zeros(nGroups, 1);
results.mean_K2_trend  = zeros(nGroups, 1);
results.mean_K2_plunge = zeros(nGroups, 1);
results.mean_K3_trend  = zeros(nGroups, 1);
results.mean_K3_plunge = zeros(nGroups, 1);

% Jelinek 誤差
results.jelinek_E12 = zeros(nGroups, 1);
results.jelinek_E23 = zeros(nGroups, 1);
results.jelinek_E13 = zeros(nGroups, 1);
results.jelinek_K1_err = zeros(nGroups, 1);
results.jelinek_K2_err = zeros(nGroups, 1);
results.jelinek_K3_err = zeros(nGroups, 1);

% Bootstrap 誤差
results.boot_K1_std  = zeros(nGroups, 1);
results.boot_K2_std  = zeros(nGroups, 1);
results.boot_K3_std  = zeros(nGroups, 1);
results.boot_V1_conf = zeros(nGroups, 1);   % K1 方向 95% 信心角
results.boot_V3_conf = zeros(nGroups, 1);   % K3 方向 95% 信心角

fprintf('\n開始逐組計算...\n');
fprintf('%-8s %-8s %-10s %-10s %-10s  Jel_E12  Jel_E23  Boot_V1conf  Boot_V3conf\n', ...
    'Number','N','mean_K1','mean_K2','mean_K3');
fprintf('%s\n', repmat('-',1,80));

for g = 1:nGroups
    num = uniqueNums(g);
    idx = (Numbers == num);
    N   = sum(idx);

    results.Number(g)    = num;
    results.N_samples(g) = N;

    % ------- 2. 計算各樣本 Kg 並標準化 ----------------------
    tensorStack = zeros(3, 3, N);   % 3×3×N
    sampleIdx   = find(idx);

    % 收集本組所有樣本編號（含字母後綴，如 A/B/C/D）
    groupNames = Names(sampleIdx);
    results.SampleNames{g} = strjoin(groupNames, ', ');

    for s = 1:N
        si = sampleIdx(s);

        % 計算特徵向量矩陣 V
        V = AMSFormulas.calculateV( ...
            dK1_all(si), iK1_all(si), ...
            dK2_all(si), iK2_all(si), ...
            dK3_all(si), iK3_all(si));

        % 計算地理座標系磁化率張量
        %   Kg = V * diag([K1,K2,K3]) * V'
        Kdiag = diag([K1_all(si), K2_all(si), K3_all(si)]);
        Kg    = V * Kdiag * V';

        % 標準化：消除強度差異，保留形狀與方向
        Kg_norm = AMSFormulas.standardizeTensor(Kg);

        % 對稱化（防止浮點誤差累積）
        tensorStack(:,:,s) = AMSFormulas.symmetrize(Kg_norm);
    end

    % ------- 3. 計算平均張量與特徵值分析 ---------------------
    M_mean   = AMSFormulas.calculateMeanTensor(tensorStack);
    eigenRes = AMSFormulas.performEigenAnalysis(M_mean);

    results.mean_K1(g) = eigenRes.K1_val;
    results.mean_K2(g) = eigenRes.K2_val;
    results.mean_K3(g) = eigenRes.K3_val;
    results.mean_K1_trend(g)  = eigenRes.K1_dir(1);
    results.mean_K1_plunge(g) = eigenRes.K1_dir(2);
    results.mean_K2_trend(g)  = eigenRes.K2_dir(1);
    results.mean_K2_plunge(g) = eigenRes.K2_dir(2);
    results.mean_K3_trend(g)  = eigenRes.K3_dir(1);
    results.mean_K3_plunge(g) = eigenRes.K3_dir(2);

    % ------- 4. Jelinek (1978) 統計誤差 ----------------------
    jelinekStats = AMSFormulas.calculateJelinekStats(tensorStack);

    results.jelinek_E12(g)    = jelinekStats.E12;
    results.jelinek_E23(g)    = jelinekStats.E23;
    results.jelinek_E13(g)    = jelinekStats.E13;
    results.jelinek_K1_err(g) = jelinekStats.K1_err;
    results.jelinek_K2_err(g) = jelinekStats.K2_err;
    results.jelinek_K3_err(g) = jelinekStats.K3_err;

    % ------- 5. Bootstrap 統計誤差 ---------------------------
    bootStats = AMSFormulas.calculateBootstrapStats(tensorStack, 1000);

    % 安全讀取（AMSFormulas early-return 時可能缺 K2_std）
    results.boot_K1_std(g)  = getFieldSafe(bootStats, 'K1_std');
    results.boot_K2_std(g)  = getFieldSafe(bootStats, 'K2_std');
    results.boot_K3_std(g)  = getFieldSafe(bootStats, 'K3_std');
    results.boot_V1_conf(g) = getFieldSafe(bootStats, 'V1_conf');
    results.boot_V3_conf(g) = getFieldSafe(bootStats, 'V3_conf');

    % 進度輸出
    fprintf('%-8d %-8d %-10.4f %-10.4f %-10.4f  %7.2f  %7.2f  %11.2f  %11.2f\n', ...
        num, N, ...
        eigenRes.K1_val, eigenRes.K2_val, eigenRes.K3_val, ...
        jelinekStats.E12, jelinekStats.E23, ...
        bootStats.V1_conf, bootStats.V3_conf);
end

fprintf('\n計算完成！共處理 %d 組\n', nGroups);

%% -------- 6. 輸出結果到 Excel --------------------------------
outputFile = fullfile(path, 'AMS_Results.xlsx');

T = table( ...
    results.Number,    results.N_samples, ...
    results.SampleNames, ...
    results.mean_K1,   results.mean_K2,   results.mean_K3, ...
    results.mean_K1_trend,  results.mean_K1_plunge, ...
    results.mean_K2_trend,  results.mean_K2_plunge, ...
    results.mean_K3_trend,  results.mean_K3_plunge, ...
    results.jelinek_E12, results.jelinek_E23, results.jelinek_E13, ...
    results.jelinek_K1_err, results.jelinek_K2_err, results.jelinek_K3_err, ...
    results.boot_K1_std, results.boot_K2_std, results.boot_K3_std, ...
    results.boot_V1_conf, results.boot_V3_conf, ...
    'VariableNames', { ...
        'Number', 'N_Samples', ...
        'Sample_Names', ...
        'Mean_K1', 'Mean_K2', 'Mean_K3', ...
        'K1_Trend', 'K1_Plunge', ...
        'K2_Trend', 'K2_Plunge', ...
        'K3_Trend', 'K3_Plunge', ...
        'Jelinek_E12_deg', 'Jelinek_E23_deg', 'Jelinek_E13_deg', ...
        'Jelinek_K1_SE',   'Jelinek_K2_SE',   'Jelinek_K3_SE', ...
        'Boot_K1_Std', 'Boot_K2_Std', 'Boot_K3_Std', ...
        'Boot_V1_Conf95_deg', 'Boot_V3_Conf95_deg' ...
    });

writetable(T, outputFile, 'Sheet', 'AMS_Results');
fprintf('結果已輸出至：%s\n', outputFile);

%% -------- 7. 視覺化（修正版）---------------------------------
% 圖1：平均主磁化率（K1/K2/K3）＋ Jelinek 標準誤差（SE，相同單位）
% 圖2：Jelinek 95% 信心角（E12/E23/E13，單位：度）
% 圖3：Bootstrap 95% 信心角（V1/V3 方向，單位：度）
%
% 修正說明：
%   原版把角度誤差（E12，最大90°）當作磁化率的誤差棒，單位完全不匹配。
%   現在：
%     圖1 的誤差棒改用 Jelinek_K1_SE / K2_SE / K3_SE（與磁化率同單位）
%     Jelinek 信心角另開圖2 單獨顯示
% =========================================================

% ---- 圖1：平均主磁化率 ± Jelinek SE ----------------------------
figure('Name','Mean Susceptibility ± Jelinek SE', ...
       'NumberTitle','off', 'Color','white', 'Position',[100 100 900 420]);
hold on;
errorbar(results.Number, results.mean_K1, results.jelinek_K1_err, ...
    'o-', 'Color',[0.12 0.47 0.71], 'MarkerFaceColor',[0.12 0.47 0.71], ...
    'DisplayName','K_1 ± SE', 'LineWidth',1.2, 'CapSize',4);
errorbar(results.Number, results.mean_K2, results.jelinek_K2_err, ...
    's-', 'Color',[0.17 0.63 0.17], 'MarkerFaceColor',[0.17 0.63 0.17], ...
    'DisplayName','K_2 ± SE', 'LineWidth',1.2, 'CapSize',4);
errorbar(results.Number, results.mean_K3, results.jelinek_K3_err, ...
    '^-', 'Color',[0.84 0.15 0.16], 'MarkerFaceColor',[0.84 0.15 0.16], ...
    'DisplayName','K_3 ± SE', 'LineWidth',1.2, 'CapSize',4);
xlabel('Number'); ylabel('Normalized Susceptibility');
title('Mean Susceptibility with Jelinek Standard Error');
legend('Location','best'); grid on; box on;
xlim([min(results.Number)-1, max(results.Number)+1]);

% ---- 圖2：Jelinek 95% 信心角（角度）----------------------------
figure('Name','Jelinek 95% Confidence Angles', ...
       'NumberTitle','off', 'Color','white', 'Position',[100 560 900 380]);
hold on;
plot(results.Number, results.jelinek_E12, 'o-', ...
    'Color',[0.12 0.47 0.71], 'MarkerFaceColor',[0.12 0.47 0.71], ...
    'DisplayName','E_{12}  (K_1–K_2)', 'LineWidth',1.2, 'MarkerSize',5);
plot(results.Number, results.jelinek_E23, 's-', ...
    'Color',[0.17 0.63 0.17], 'MarkerFaceColor',[0.17 0.63 0.17], ...
    'DisplayName','E_{23}  (K_2–K_3)', 'LineWidth',1.2, 'MarkerSize',5);
plot(results.Number, results.jelinek_E13, '^-', ...
    'Color',[0.84 0.15 0.16], 'MarkerFaceColor',[0.84 0.15 0.16], ...
    'DisplayName','E_{13}  (K_1–K_3)', 'LineWidth',1.2, 'MarkerSize',5);
% 標示 45° 警戒線（重疊區，信心橢圓失去意義）
yline(45,'--','Color',[0.6 0.6 0.6],'LineWidth',1,'Label','45° (overlap warning)');
xlabel('Number'); ylabel('Confidence Angle (°)');
title('Jelinek 95% Confidence Angles (E_{12}, E_{23}, E_{13})');
legend('Location','best'); grid on; box on;
ylim([0 90]);
xlim([min(results.Number)-1, max(results.Number)+1]);

% ---- 圖3：Bootstrap 95% 信心角（角度）--------------------------
figure('Name','Bootstrap 95% Confidence Angles', ...
       'NumberTitle','off', 'Color','white', 'Position',[100 100 900 380]);
hold on;
plot(results.Number, results.boot_V1_conf, 'o-', ...
    'Color',[0.12 0.47 0.71], 'MarkerFaceColor',[0.12 0.47 0.71], ...
    'DisplayName','K_1 direction (95° cone)', 'LineWidth',1.2, 'MarkerSize',5);
plot(results.Number, results.boot_V3_conf, '^-', ...
    'Color',[0.84 0.15 0.16], 'MarkerFaceColor',[0.84 0.15 0.16], ...
    'DisplayName','K_3 direction (95° cone)', 'LineWidth',1.2, 'MarkerSize',5);
yline(45,'--','Color',[0.6 0.6 0.6],'LineWidth',1,'Label','45° (overlap warning)');
xlabel('Number'); ylabel('Confidence Angle (°)');
title('Bootstrap 95% Confidence Angles for K_1 and K_3 Directions');
legend('Location','best'); grid on; box on;
ylim([0 90]);
xlim([min(results.Number)-1, max(results.Number)+1]);
%% -------- 本地輔助函式 ----------------------------------------
function val = getFieldSafe(s, fieldName)
    % 安全讀取 struct 欄位，欄位不存在時回傳 NaN
    if isfield(s, fieldName)
        val = s.(fieldName);
    else
        val = NaN;
    end
end