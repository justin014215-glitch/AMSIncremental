%% AMS_IncrementalStrain.m
% 應變增量計算主程式
%
% 前置需求：先執行 AMS_StrainPrep.m → 產生 groupEg.mat
%
% 目前功能：
%   [模式 1] 選擇任意兩組 Group，計算應變增量 F = Ef2 * Ef1^(-1)
%   [模式 2] 依分區(A/B/C/D)計算應變增量 ← 分區代表值計算方式待定
%
% 依賴：AMSFormulas.m
% =========================================================

clearvars; clc;

fprintf('============================================================\n');
fprintf('   AMS 應變增量計算主程式\n');
fprintf('============================================================\n\n');

%% -------- 1. 載入 groupEg.mat -----------------------------------
[matFile, matPath] = uigetfile('*.mat', '選擇 groupEg.mat 檔案');
if isequal(matFile, 0), error('未選擇檔案，程式結束。'); end
load(fullfile(matPath, matFile), 'groupEg', 'groupNums', 'results');
nGroups = size(groupEg, 3);

fprintf('已載入 %d 組 Eg 張量\n', nGroups);

% 顯示各組資訊
fprintf('\n%-8s %-12s %-12s %-12s  %-18s\n', ...
    'Number', 'Eg_e1', 'Eg_e2', 'Eg_e3', 'Sample_Names');
fprintf('%s\n', repmat('-', 1, 70));
for g = 1:nGroups
    EgRes = AMSFormulas.performEigenAnalysis(groupEg(:,:,g));
    name_str = '';
    if isfield(results, 'SampleNames') && g <= numel(results.SampleNames)
        name_str = results.SampleNames{g};
    end
    fprintf('%-8d %-12.6f %-12.6f %-12.6f  %s\n', ...
        groupNums(g), EgRes.K1_val, EgRes.K2_val, EgRes.K3_val, name_str);
end
fprintf('\n');

%% -------- 2. 選擇計算模式 ----------------------------------------
modeChoice = questdlg('選擇計算模式', '應變增量模式', ...
    '模式1：選兩個樣本', '模式2：依分區(A/B/C/D)', '模式1：選兩個樣本');

if isempty(modeChoice), disp('已取消'); return; end

%% ================================================================
%  模式 1：選擇任意兩個 Group 計算應變增量
%% ================================================================
if strcmp(modeChoice, '模式1：選兩個樣本')

    keepGoing    = true;
    allPairResults = {};

    while keepGoing
        fprintf('\n--- 配對 #%d ---\n', numel(allPairResults)+1);

        % 選擇初始 Group（Ef1）
        numList = arrayfun(@(x) num2str(x), groupNums, 'UniformOutput', false);
        [idx1, ok] = listdlg( ...
            'PromptString', '選擇初始狀態 Group（Ef1，較早期/西側）：', ...
            'SelectionMode', 'single', ...
            'ListString', numList, ...
            'Name', '選擇 Ef1');
        if ~ok, fprintf('已取消\n'); break; end

        % 選擇最終 Group（Ef2）
        [idx2, ok] = listdlg( ...
            'PromptString', '選擇最終狀態 Group（Ef2，較晚期/東側）：', ...
            'SelectionMode', 'single', ...
            'ListString', numList, ...
            'Name', '選擇 Ef2');
        if ~ok, fprintf('已取消\n'); break; end

        num1 = groupNums(idx1);
        num2 = groupNums(idx2);
        Ef1  = groupEg(:,:,idx1);
        Ef2  = groupEg(:,:,idx2);

        pairRes = computeIncrementalStrain(Ef1, Ef2, num1, num2);
        allPairResults{end+1} = pairRes; %#ok<AGROW>

        cont = questdlg('是否再計算一組配對？', '繼續？', '繼續', '結束', '結束');
        if strcmp(cont, '結束'), keepGoing = false; end
    end

    if ~isempty(allPairResults)
        outputResults(allPairResults, matPath, 'PairResults');
    end

end % 模式 1

%% ================================================================
%  模式 2：依分區(A/B/C/D)計算應變增量
%  ⚠ 各分區代表 Eg 的計算方式尚未決定，目前以分量平均佔位
%% ================================================================
if strcmp(modeChoice, '模式2：依分區(A/B/C/D)')

    fprintf('\n============================\n');
    fprintf('  模式 2：依分區計算應變增量\n');
    fprintf('============================\n');
    fprintf('⚠  各分區代表 Eg 目前以「分量平均」計算，待確認後修改 calcSectionRepEg()\n\n');

    %% ---- 2a. 設定各分區包含的 Group Numbers --------------------
    % ⚠ 請依彭筱君(2015)分區結果填入對應 Number
    % 若保持空陣列 []，程式會彈出對話框讓你手動輸入
    sectionGroups.A = [];   % 例如 [1, 2, 3]
    sectionGroups.B = [];   % 例如 [4, 5, 6]
    sectionGroups.C = [];   % 例如 [7, 8, 9]
    sectionGroups.D = [];   % 例如 [10, 11, 12]

    sections = {'A','B','C','D'};
    for s = 1:4
        secName = sections{s};
        if isempty(sectionGroups.(secName))
            ans_str = inputdlg( ...
                {sprintf('區段 %s 包含的 Group Numbers（空格分隔，可留空跳過）：', secName)}, ...
                sprintf('設定區段 %s', secName), 1, {''});
            if ~isempty(ans_str) && ~isempty(strtrim(ans_str{1}))
                sectionGroups.(secName) = str2num(ans_str{1}); %#ok<ST2NM>
            end
        end
    end

    %% ---- 2b. 計算各分區代表 Eg ----------------------------------
    sectionEg    = struct();
    sectionValid = struct();

    for s = 1:4
        secName = sections{s};
        nums    = sectionGroups.(secName);

        if isempty(nums)
            fprintf('區段 %s：未設定，跳過\n', secName);
            sectionValid.(secName) = false;
            continue;
        end

        [Eg_rep, ok] = calcSectionRepEg(nums, groupNums, groupEg);

        if ok
            sectionEg.(secName)    = Eg_rep;
            sectionValid.(secName) = true;
            EgRes = AMSFormulas.performEigenAnalysis(Eg_rep);
            fprintf('區段 %s (Groups: %s)  Eg_e1=%.6f  e3=%.6f\n', ...
                secName, num2str(nums), EgRes.K1_val, EgRes.K3_val);
        else
            fprintf('區段 %s：Group 不存在，跳過\n', secName);
            sectionValid.(secName) = false;
        end
    end

    %% ---- 2c. 計算相鄰分區間應變增量 ----------------------------
    pairs = {'A','B'; 'B','C'; 'C','D'; 'A','D'};
    pairResults2 = {};

    for p = 1:size(pairs, 1)
        s1 = pairs{p,1};  s2 = pairs{p,2};
        if ~sectionValid.(s1) || ~sectionValid.(s2)
            fprintf('配對 %s→%s：資料不足，跳過\n', s1, s2);
            continue;
        end
        label   = sprintf('%s→%s', s1, s2);
        pairRes = computeIncrementalStrain( ...
            sectionEg.(s1), sectionEg.(s2), s1, s2, label);
        pairResults2{end+1} = pairRes; %#ok<AGROW>
    end

    if ~isempty(pairResults2)
        outputResults(pairResults2, matPath, 'SectionResults');
        % 額外記錄分區設定
        outputFile = fullfile(matPath, 'IncrementalStrain_SectionResults.xlsx');
        sNames = {}; sNums = {};
        for s = 1:4
            sn = sections{s};
            sNames{end+1} = sn; %#ok<AGROW>
            sNums{end+1}  = num2str(sectionGroups.(sn)); %#ok<AGROW>
        end
        writetable( ...
            table(sNames', sNums', 'VariableNames', {'Section','Group_Numbers'}), ...
            outputFile, 'Sheet', 'SectionDefinition');
    end

end % 模式 2

disp('程式執行完畢。');


%% ================================================================
%  核心計算函式
%% ================================================================

function res = computeIncrementalStrain(Ef1, Ef2, label1, label2, label)
% 計算應變增量 F = Ef2 * Ef1^(-1) 並做特徵值分析
%
% F 一般非對稱（含旋轉），此處取對稱部分做特徵值分析
% 若後續要分析旋轉成分可用 F 直接做極分解

    if nargin < 5
        label = sprintf('%s → %s', num2str(label1), num2str(label2));
    end

    fprintf('\n[%s]\n', label);

    % 計算增量
    F     = AMSFormulas.calculateF(Ef1, Ef2);   % F = Ef2 * Ef1^(-1)
    F_sym = AMSFormulas.symmetrize(F);           % 取對稱部分

    % 特徵值分析
    FRes  = AMSFormulas.performEigenAnalysis(F_sym);

    % 形狀參數 T
    T_val = NaN;
    if all([FRes.K1_val, FRes.K2_val, FRes.K3_val] > 0)
        lnR1 = log(FRes.K1_val / FRes.K2_val);
        lnR2 = log(FRes.K2_val / FRes.K3_val);
        if (lnR1 + lnR2) ~= 0
            T_val = (lnR2 - lnR1) / (lnR2 + lnR1);
        end
    end

    % 輸出
    fprintf('  det(F) = %.6f（體積守恆應 ≈ 1）\n', det(F));
    fprintf('  增量主應變：\n');
    fprintf('    e1 = %10.6f  Trend/Plunge: %6.1f / %5.1f\n', ...
        FRes.K1_val, FRes.K1_dir(1), FRes.K1_dir(2));
    fprintf('    e2 = %10.6f  Trend/Plunge: %6.1f / %5.1f\n', ...
        FRes.K2_val, FRes.K2_dir(1), FRes.K2_dir(2));
    fprintf('    e3 = %10.6f  Trend/Plunge: %6.1f / %5.1f\n', ...
        FRes.K3_val, FRes.K3_dir(1), FRes.K3_dir(2));
    if ~isnan(T_val)
        shape = '中性';
        if T_val >  0.1, shape = 'Oblate 扁平狀';
        elseif T_val < -0.1, shape = 'Prolate 雪茄狀'; end
        fprintf('  L=%.4f  F=%.4f  T=%.4f  (%s)\n', FRes.L, FRes.F, T_val, shape);
    end

    % 打包
    res.label  = label;
    res.label1 = label1;
    res.label2 = label2;
    res.F      = F;
    res.F_sym  = F_sym;
    res.FRes   = FRes;
    res.T_val  = T_val;
    res.detF   = det(F);
end


function [Eg_rep, ok] = calcSectionRepEg(nums, groupNums, groupEg)
% ⚠ 各分區代表 Eg 計算（目前：分量平均）
% 待確認後可改為其他方法：
%   - 選取代表性單一樣本
%   - 加權平均
%   - bootstrap 中位張量

    tensorStack = [];
    count = 0;
    for i = 1:numel(nums)
        gIdx = find(groupNums == nums(i));
        if isempty(gIdx)
            warning('Group %d 不存在，跳過', nums(i));
            continue;
        end
        count = count + 1;
        tensorStack(:,:,count) = groupEg(:,:,gIdx); %#ok<AGROW>
    end

    if count == 0
        Eg_rep = nan(3,3); ok = false;
    elseif count == 1
        Eg_rep = tensorStack(:,:,1); ok = true;
    else
        % ⚠ 目前：分量平均（待確認）
        Eg_rep = mean(tensorStack, 3);
        Eg_rep = AMSFormulas.symmetrize(Eg_rep);
        ok = true;
    end
end


%% ================================================================
%  輸出函式
%% ================================================================

function outputResults(pairResults, matPath, tag)
% 彙整輸出 Excel、.mat、Stereonet

    nP = numel(pairResults);

    % 組裝 Table
    Label     = cellfun(@(r) r.label,            pairResults, 'UniformOutput', false)';
    Ef1_ID    = cellfun(@(r) num2str(r.label1),  pairResults, 'UniformOutput', false)';
    Ef2_ID    = cellfun(@(r) num2str(r.label2),  pairResults, 'UniformOutput', false)';
    det_F     = cellfun(@(r) r.detF,             pairResults)';
    F_e1      = cellfun(@(r) r.FRes.K1_val,      pairResults)';
    F_e2      = cellfun(@(r) r.FRes.K2_val,      pairResults)';
    F_e3      = cellfun(@(r) r.FRes.K3_val,      pairResults)';
    e1_Trend  = cellfun(@(r) r.FRes.K1_dir(1),   pairResults)';
    e1_Plunge = cellfun(@(r) r.FRes.K1_dir(2),   pairResults)';
    e2_Trend  = cellfun(@(r) r.FRes.K2_dir(1),   pairResults)';
    e2_Plunge = cellfun(@(r) r.FRes.K2_dir(2),   pairResults)';
    e3_Trend  = cellfun(@(r) r.FRes.K3_dir(1),   pairResults)';
    e3_Plunge = cellfun(@(r) r.FRes.K3_dir(2),   pairResults)';
    L_val     = cellfun(@(r) r.FRes.L,           pairResults)';
    F_val     = cellfun(@(r) r.FRes.F,           pairResults)';
    T_shape   = cellfun(@(r) r.T_val,            pairResults)';

    T = table(Label, Ef1_ID, Ef2_ID, det_F, ...
        F_e1, F_e2, F_e3, ...
        e1_Trend, e1_Plunge, e2_Trend, e2_Plunge, e3_Trend, e3_Plunge, ...
        L_val, F_val, T_shape, ...
        'VariableNames', { ...
        'Label','Ef1_ID','Ef2_ID','det_F', ...
        'F_e1','F_e2','F_e3', ...
        'e1_Trend','e1_Plunge','e2_Trend','e2_Plunge','e3_Trend','e3_Plunge', ...
        'L_lineation','F_foliation','T_shape'});

    outputFile = fullfile(matPath, ['IncrementalStrain_' tag '.xlsx']);
    writetable(T, outputFile, 'Sheet', tag);
    fprintf('\nExcel 已輸出：%s\n', outputFile);

    % .mat
    save(fullfile(matPath, ['IncrementalStrain_' tag '.mat']), 'pairResults');

    % Stereonet
    stereoDir = fullfile(matPath, ['IncrementalStereonets_' tag]);
    if ~exist(stereoDir, 'dir'), mkdir(stereoDir); end

    for p = 1:nP
        res  = pairResults{p};
        FRes = res.FRes;

        figure('Name', ['應變增量 - ' res.label], ...
               'Color', 'white', 'Position', [100 100 650 650]);
        Stereonet(0, pi/2, 10*pi/180, 1);
        ax = gca; hold(ax, 'on');

        [x1,y1] = schmidtProject(FRes.K1_dir(1), FRes.K1_dir(2));
        plot(ax, x1, y1, 's', 'MarkerFaceColor', [0.12 0.47 0.71], ...
            'MarkerEdgeColor', [0.05 0.25 0.50], 'MarkerSize', 14, 'LineWidth', 1.5);
        text(ax, x1+0.05, y1+0.05, sprintf('e1=%.4f', FRes.K1_val), ...
            'FontSize', 9, 'Color', [0.12 0.47 0.71]);

        [x2,y2] = schmidtProject(FRes.K2_dir(1), FRes.K2_dir(2));
        plot(ax, x2, y2, '^', 'MarkerFaceColor', [0.17 0.63 0.17], ...
            'MarkerEdgeColor', [0.05 0.38 0.05], 'MarkerSize', 14, 'LineWidth', 1.5);

        [x3,y3] = schmidtProject(FRes.K3_dir(1), FRes.K3_dir(2));
        plot(ax, x3, y3, 'o', 'MarkerFaceColor', [0.84 0.15 0.16], ...
            'MarkerEdgeColor', [0.55 0.05 0.05], 'MarkerSize', 14, 'LineWidth', 1.5);
        text(ax, x3+0.05, y3+0.05, sprintf('e3=%.4f', FRes.K3_val), ...
            'FontSize', 9, 'Color', [0.84 0.15 0.16]);

        h1 = plot(ax, NaN, NaN, 's', 'MarkerFaceColor', [0.12 0.47 0.71], ...
            'MarkerEdgeColor',[0.05 0.25 0.50], 'MarkerSize', 12, ...
            'DisplayName', 'e_1 (max extension)');
        h2 = plot(ax, NaN, NaN, '^', 'MarkerFaceColor', [0.17 0.63 0.17], ...
            'MarkerEdgeColor',[0.05 0.38 0.05], 'MarkerSize', 12, ...
            'DisplayName', 'e_2');
        h3 = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor', [0.84 0.15 0.16], ...
            'MarkerEdgeColor',[0.55 0.05 0.05], 'MarkerSize', 12, ...
            'DisplayName', 'e_3 (max shortening)');
        legend(ax, [h1 h2 h3], 'Location', 'southoutside', ...
            'Orientation', 'horizontal', 'FontSize', 10);

        title(ax, sprintf('應變增量 — %s\nT=%.3f  det(F)=%.4f', ...
            res.label, res.T_val, res.detF), ...
            'FontSize', 12, 'FontWeight', 'bold');
        axis(ax, 'off'); hold(ax, 'off');

        fname = sprintf('IncrementalStrain_%s_%s.png', tag, ...
            strrep(strrep(res.label, '→', 'to'), ' ', ''));
        exportgraphics(gcf, fullfile(stereoDir, fname), 'Resolution', 300);
        fprintf('Stereonet saved: %s\n', fname);
    end
    fprintf('圖檔已儲存至：%s\n', stereoDir);
end


function [x, y] = schmidtProject(trend_deg, plunge_deg)
    r = sqrt(2) * sin((pi/2 - deg2rad(plunge_deg)) / 2);
    x = r * sin(deg2rad(trend_deg));
    y = r * cos(deg2rad(trend_deg));
end