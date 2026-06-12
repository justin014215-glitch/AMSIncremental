%% AMS_StrainPrep.m
% 主程式：依 Number 分組，計算各組平均應變張量 Eg
% 【不做標準化】，保留絕對磁化率強度，供後續應變增量計算使用
%
% 流程：
%   1. 各樣本：計算 V（方向餘弦）→ Kg = V * diag([K1,K2,K3]) * V'（地理座標）
%   2. 同組所有 Kg 取分量平均 → 平均 Kg（地理座標）
%   3. 對平均 Kg 做特徵值分解 → 平均主磁化率 [K1m,K2m,K3m] 和平均主軸方向 V_mean
%   4. Slate 模型：[K1m,K2m,K3m] → Er（主軸座標系應變張量，對角矩陣）
%   5. 座標轉換：Eg = V_mean * Er * V_mean'（回到地理座標系）
%   6. Jelinek & Bootstrap 誤差
%
% 資料欄位說明：
%   Number    - 站點編號（數值，用於分組）
%   Group     - 分區標籤（字元，如 A/B/C/D，供 Stereonet 篩選用）
%   Name      - 樣本名稱
%   Spec Name - 方磚名稱
%
% 輸出：
%   - AMS_StrainPrep_Results.xlsx  (各組 Eg 特徵值/方向/誤差)
%   - groupEg.mat                  (各組 3×3 Eg 張量，供應變增量主程式使用)
%
% 依賴：AMSFormulas.m, plotStereonet.m
% =========================================================

clearvars; clc;

%% -------- 0. 參數設定 ----------------------------------------
% Slate 模型係數 (Borradaile & Jackson 2004)
slateA = 6.897;
slateB = 0.007;

%% -------- 1. 讀取資料 ----------------------------------------
[file, path] = uigetfile('*.xlsx', '選擇原始資料 Excel 檔案');
if isequal(file, 0)
    error('未選擇檔案，程式結束。');
end
rawData = readtable(fullfile(path, file), 'VariableNamingRule', 'preserve');

fprintf('已讀取 %d 筆資料\n', height(rawData));

% 取出欄位
Numbers   = rawData{:, 'Number'};
Names     = rawData{:, 'Name'};
K1_all    = rawData{:, 'K1'};
K2_all    = rawData{:, 'K2'};
K3_all    = rawData{:, 'K3'};
dK1_all   = rawData{:, 'dK1geo'};
iK1_all   = rawData{:, 'iK1geo'};
dK2_all   = rawData{:, 'dK2geo'};
iK2_all   = rawData{:, 'iK2geo'};
dK3_all   = rawData{:, 'dK3geo'};
iK3_all   = rawData{:, 'iK3geo'};

% 讀取 Group 分區欄位（A/B/C/D），清理可能的引號與空白
if ismember('Group', rawData.Properties.VariableNames)
    Groups_raw = rawData{:, 'Group'};
    % 統一清理：去除引號與前後空白
    Groups_all = cellfun(@(s) strtrim(strrep(s, '''', '')), ...
                         Groups_raw, 'UniformOutput', false);
    hasGroupCol = true;
    fprintf('已讀取 Group 分區欄位（%s）\n', ...
        strjoin(unique(Groups_all, 'stable'), ', '));
else
    Groups_all  = repmat({''}, height(rawData), 1);
    hasGroupCol = false;
    fprintf('警告：未找到 Group 欄位，分區篩選功能將無法使用\n');
end

%% -------- 2. 依 Number 分組處理 ------------------------------
uniqueNums = unique(Numbers, 'stable');
nGroups    = numel(uniqueNums);

% 結果儲存
results = struct();
results.Number      = zeros(nGroups, 1);
results.N_samples   = zeros(nGroups, 1);
results.SampleNames = cell(nGroups, 1);
results.GroupLabel  = cell(nGroups, 1);   % ← 新增：各站點的 Group 標籤

% 平均 Kg 特徵值與方向 (磁化率空間)
results.mean_K1        = zeros(nGroups, 1);
results.mean_K2        = zeros(nGroups, 1);
results.mean_K3        = zeros(nGroups, 1);
results.mean_K1_trend  = zeros(nGroups, 1);
results.mean_K1_plunge = zeros(nGroups, 1);
results.mean_K2_trend  = zeros(nGroups, 1);
results.mean_K2_plunge = zeros(nGroups, 1);
results.mean_K3_trend  = zeros(nGroups, 1);
results.mean_K3_plunge = zeros(nGroups, 1);

% 應變張量 Eg 特徵值與方向 (應變空間)
results.Eg_e1        = zeros(nGroups, 1);
results.Eg_e2        = zeros(nGroups, 1);
results.Eg_e3        = zeros(nGroups, 1);
results.Eg_e1_trend  = zeros(nGroups, 1);
results.Eg_e1_plunge = zeros(nGroups, 1);
results.Eg_e2_trend  = zeros(nGroups, 1);
results.Eg_e2_plunge = zeros(nGroups, 1);
results.Eg_e3_trend  = zeros(nGroups, 1);
results.Eg_e3_plunge = zeros(nGroups, 1);

% Jelinek 誤差
results.jelinek_E12    = zeros(nGroups, 1);
results.jelinek_E23    = zeros(nGroups, 1);
results.jelinek_E13    = zeros(nGroups, 1);
results.jelinek_K1_err = zeros(nGroups, 1);
results.jelinek_K2_err = zeros(nGroups, 1);
results.jelinek_K3_err = zeros(nGroups, 1);

% Bootstrap 誤差
results.boot_K1_std  = zeros(nGroups, 1);
results.boot_K2_std  = zeros(nGroups, 1);
results.boot_K3_std  = zeros(nGroups, 1);
results.boot_V1_conf = zeros(nGroups, 1);
results.boot_V3_conf = zeros(nGroups, 1);

% 儲存各組完整 Eg 張量
groupEg   = zeros(3, 3, nGroups);
groupNums = zeros(nGroups, 1);

fprintf('\n開始逐組計算...\n');
fprintf('%-6s %-5s %-5s %-12s %-12s %-12s  %-10s %-10s  Jel_E12  Boot_V1conf\n', ...
    'Number','Grp','N','mean_K1','mean_K2','mean_K3','Eg_e1','Eg_e3');
fprintf('%s\n', repmat('-',1,95));

for g = 1:nGroups
    num = uniqueNums(g);
    idx = (Numbers == num);
    N   = sum(idx);

    results.Number(g)    = num;
    results.N_samples(g) = N;

    % 取得本組行號與樣本名稱
    sampleIdx  = find(idx);
    groupNames = Names(sampleIdx);
    results.SampleNames{g} = strjoin(groupNames, ', ');

    % 取得本組 Group 標籤（同一 Number 下通常同一 Group，取第一筆）
    if hasGroupCol
        grpLabel = Groups_all{sampleIdx(1)};
    else
        grpLabel = '';
    end
    results.GroupLabel{g} = grpLabel;

    % ------- (a) 建構各樣本 Kg（不標準化）---------------------
    tensorStack = zeros(3, 3, N);

    for s = 1:N
        si = sampleIdx(s);
        V = AMSFormulas.calculateV( ...
            dK1_all(si), iK1_all(si), ...
            dK2_all(si), iK2_all(si), ...
            dK3_all(si), iK3_all(si));
        Kdiag = diag([K1_all(si), K2_all(si), K3_all(si)]);
        Kg    = V * Kdiag * V';
        tensorStack(:,:,s) = AMSFormulas.symmetrize(Kg);
    end

    % ------- (b) 平均 Kg → 特徵值分解 -------------------------
    M_mean   = AMSFormulas.calculateMeanTensor(tensorStack);
    eigenRes = AMSFormulas.performEigenAnalysis(M_mean);

    results.mean_K1(g)        = eigenRes.K1_val;
    results.mean_K2(g)        = eigenRes.K2_val;
    results.mean_K3(g)        = eigenRes.K3_val;
    results.mean_K1_trend(g)  = eigenRes.K1_dir(1);
    results.mean_K1_plunge(g) = eigenRes.K1_dir(2);
    results.mean_K2_trend(g)  = eigenRes.K2_dir(1);
    results.mean_K2_plunge(g) = eigenRes.K2_dir(2);
    results.mean_K3_trend(g)  = eigenRes.K3_dir(1);
    results.mean_K3_plunge(g) = eigenRes.K3_dir(2);

    % ------- (c) 平均主磁化率 → Er（主軸座標系）--------------
    V_mean = eigenRes.V;
    Er = AMSFormulas.calculateEr( ...
        eigenRes.K1_val, eigenRes.K2_val, eigenRes.K3_val, ...
        slateA, slateB);

    % ------- (d) Er → Eg（地理座標）--------------------------
    Eg = AMSFormulas.calculateEg(Er, V_mean);
    Eg = AMSFormulas.symmetrize(Eg);

    groupEg(:,:,g) = Eg;
    groupNums(g)   = num;

    % Eg 特徵值分析
    EgRes = AMSFormulas.performEigenAnalysis(Eg);
    results.Eg_e1(g)        = EgRes.K1_val;
    results.Eg_e2(g)        = EgRes.K2_val;
    results.Eg_e3(g)        = EgRes.K3_val;
    results.Eg_e1_trend(g)  = EgRes.K1_dir(1);
    results.Eg_e1_plunge(g) = EgRes.K1_dir(2);
    results.Eg_e2_trend(g)  = EgRes.K2_dir(1);
    results.Eg_e2_plunge(g) = EgRes.K2_dir(2);
    results.Eg_e3_trend(g)  = EgRes.K3_dir(1);
    results.Eg_e3_plunge(g) = EgRes.K3_dir(2);

    % ------- (e) Jelinek 誤差 ---------------------------------
    jelinekStats = AMSFormulas.calculateJelinekStats(tensorStack);
    results.jelinek_E12(g)    = jelinekStats.E12;
    results.jelinek_E23(g)    = jelinekStats.E23;
    results.jelinek_E13(g)    = jelinekStats.E13;
    results.jelinek_K1_err(g) = jelinekStats.K1_err;
    results.jelinek_K2_err(g) = jelinekStats.K2_err;
    results.jelinek_K3_err(g) = jelinekStats.K3_err;

    % ------- (f) Bootstrap 誤差 ------------------------------
    bootStats = AMSFormulas.calculateBootstrapStats(tensorStack, 1000);
    results.boot_K1_std(g)  = getFieldSafe(bootStats, 'K1_std');
    results.boot_K2_std(g)  = getFieldSafe(bootStats, 'K2_std');
    results.boot_K3_std(g)  = getFieldSafe(bootStats, 'K3_std');
    results.boot_V1_conf(g) = getFieldSafe(bootStats, 'V1_conf');
    results.boot_V3_conf(g) = getFieldSafe(bootStats, 'V3_conf');

    % 進度輸出
    fprintf('%-6d %-5s %-5d %-12.4e %-12.4e %-12.4e  %-10.6f %-10.6f  %7.2f  %11.2f\n', ...
        num, grpLabel, N, ...
        eigenRes.K1_val, eigenRes.K2_val, eigenRes.K3_val, ...
        EgRes.K1_val, EgRes.K3_val, ...
        jelinekStats.E12, getFieldSafe(bootStats, 'V1_conf'));
end

fprintf('\n計算完成！共處理 %d 組\n', nGroups);

%% -------- 3. 輸出結果到 Excel --------------------------------
outputFile = fullfile(path, 'AMS_StrainPrep_Results.xlsx');

T = table( ...
    results.Number,      results.N_samples,   results.GroupLabel,  results.SampleNames, ...
    results.mean_K1,     results.mean_K2,     results.mean_K3, ...
    results.mean_K1_trend,  results.mean_K1_plunge, ...
    results.mean_K2_trend,  results.mean_K2_plunge, ...
    results.mean_K3_trend,  results.mean_K3_plunge, ...
    results.Eg_e1,       results.Eg_e2,       results.Eg_e3, ...
    results.Eg_e1_trend, results.Eg_e1_plunge, ...
    results.Eg_e2_trend, results.Eg_e2_plunge, ...
    results.Eg_e3_trend, results.Eg_e3_plunge, ...
    results.jelinek_E12, results.jelinek_E23,  results.jelinek_E13, ...
    results.jelinek_K1_err, results.jelinek_K2_err, results.jelinek_K3_err, ...
    results.boot_K1_std, results.boot_K2_std,  results.boot_K3_std, ...
    results.boot_V1_conf, results.boot_V3_conf, ...
    'VariableNames', { ...
        'Number', 'N_Samples', 'Group', 'Sample_Names', ...
        'Mean_K1', 'Mean_K2', 'Mean_K3', ...
        'K1_Trend', 'K1_Plunge', ...
        'K2_Trend', 'K2_Plunge', ...
        'K3_Trend', 'K3_Plunge', ...
        'Eg_e1', 'Eg_e2', 'Eg_e3', ...
        'Eg_e1_Trend', 'Eg_e1_Plunge', ...
        'Eg_e2_Trend', 'Eg_e2_Plunge', ...
        'Eg_e3_Trend', 'Eg_e3_Plunge', ...
        'Jelinek_E12_deg', 'Jelinek_E23_deg', 'Jelinek_E13_deg', ...
        'Jelinek_K1_SE',   'Jelinek_K2_SE',   'Jelinek_K3_SE', ...
        'Boot_K1_Std', 'Boot_K2_Std', 'Boot_K3_Std', ...
        'Boot_V1_Conf95_deg', 'Boot_V3_Conf95_deg' ...
    });

writetable(T, outputFile, 'Sheet', 'Results');
fprintf('Excel 結果已輸出至：%s\n', outputFile);

%% -------- 4. 儲存 Eg 張量供應變增量計算 ----------------------
matFile = fullfile(path, 'groupEg.mat');
save(matFile, 'groupEg', 'groupNums', 'results');
fprintf('Eg 張量已儲存至：%s\n', matFile);
fprintf('  → 變數 groupEg     : 3×3×%d 張量陣列\n', nGroups);
fprintf('  → 變數 groupNums   : 對應 Number 編號\n');
fprintf('  → results.GroupLabel: 各站點 A/B/C/D 分區標籤\n');
fprintf('\n後續應變增量計算請載入 groupEg.mat，\n');
fprintf('並呼叫 AMSFormulas.calculateF(Eg_initial, Eg_final)\n');

disp('程式執行完畢。');

%% -------- 5. 繪製 Stereonet ----------------------------------------
% 修正：plotStereonet 需要同時傳入 results（平均主軸）與 rawData（原始方磚）
% ● Mean 模式：只畫平均主軸符號（K1/K2/K3），不畫橢圓
% ● Raw  模式：畫選定站點的個別方磚 + 平均主軸大符號 + Bootstrap/Jelinek 信心橢圓
% =========================================================

fprintf('\n========== Stereonet 輸出 ==========\n');
fprintf('可用 Number 範圍：%d ~ %d（共 %d 組）\n', ...
    min(results.Number), max(results.Number), nGroups);
if hasGroupCol
    fprintf('可用 Group 分區：%s\n', ...
        strjoin(unique(results.GroupLabel, 'stable'), ', '));
end

plotChoice = questdlg( ...
    '要繪製 Stereonet 嗎？', 'Stereonet 輸出', ...
    '互動式選擇', '畫全部', '不畫', '互動式選擇');

if ~strcmp(plotChoice, '不畫') && ~isempty(plotChoice)

    % 儲存路徑
    stereoDir = fullfile(path, 'Stereonets');
    if ~exist(stereoDir, 'dir'), mkdir(stereoDir); end

    switch plotChoice

        % ---- 互動式：彈出選單讓使用者選擇資料類型、站點、橢圓類型 ----
        case '互動式選擇'
            % results 與 rawData 都傳入，讓使用者在函式內選擇要用哪種
            plotStereonet(results, rawData, ...
                'SaveFig',  true, ...
                'SavePath', stereoDir);

        % ---- 畫全部站點（仍讓使用者選資料類型與橢圓）----
        case '畫全部'
            plotStereonet(results, rawData, ...
                'SaveFig',  true, ...
                'SavePath', stereoDir);
    end

    fprintf('Stereonet 圖檔已儲存至：%s\n', stereoDir);
end
%% -------- 批次輸出所有站點 Stereonet ----------------------
batchDir = fullfile(path, 'Stereonets_Batch');
meanDir  = fullfile(batchDir, 'Mean');
rawDir   = fullfile(batchDir, 'Raw');
strainDir = fullfile(batchDir, 'Strain');
if ~exist(meanDir,   'dir'), mkdir(meanDir);   end
if ~exist(rawDir,    'dir'), mkdir(rawDir);    end
if ~exist(strainDir, 'dir'), mkdir(strainDir); end

% 快速 renderer
set(0, 'DefaultFigureRenderer', 'painters');

fprintf('\n========== 批次輸出所有站點 Stereonet ==========\n');

% 先把需要的資料擷取成一般陣列，parfor 內不能存取 struct 的 cell 欄位
nG         = nGroups;
numList    = results.Number;
grpLabels  = results.GroupLabel;
Eg_e1_t    = results.Eg_e1_trend;
Eg_e1_p    = results.Eg_e1_plunge;
Eg_e2_t    = results.Eg_e2_trend;
Eg_e2_p    = results.Eg_e2_plunge;
Eg_e3_t    = results.Eg_e3_trend;
Eg_e3_p    = results.Eg_e3_plunge;
Eg_e1_vals = results.Eg_e1;
Eg_e3_vals = results.Eg_e3;

parfor g = 1:nG
    num = numList(g);

    % --- Mean ---
    plotStereonetBatch(results, rawData, num, 'mean', 'none', meanDir, 200);

    % --- Raw ---
    plotStereonetBatch(results, rawData, num, 'raw', 'both', rawDir, 200);

    % --- Strain ---
    plotStrainStereonetBatch( ...
        num, grpLabels{g}, ...
        Eg_e1_t(g), Eg_e1_p(g), ...
        Eg_e2_t(g), Eg_e2_p(g), ...
        Eg_e3_t(g), Eg_e3_p(g), ...
        Eg_e1_vals(g), Eg_e3_vals(g), ...
        strainDir);
end

fprintf('批次輸出完成。\n');
%% -------- 本地輔助函式 ----------------------------------------
function val = getFieldSafe(s, fieldName)
    if isfield(s, fieldName)
        val = s.(fieldName);
    else
        val = NaN;
    end
end