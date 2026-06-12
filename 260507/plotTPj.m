%% plotTPj.m
% 從 groupEg.mat 讀取應變張量，計算 T 與 Pj 並畫 Jelinek T-Pj 圖
%
% 公式（Jelinek 1981）：
%   Km = (K1·K2·K3)^(1/3)
%   Pj = exp√(2·[(ln(K1/Km))² + (ln(K2/Km))² + (ln(K3/Km))²])
%   T  = (ln(K2/K3) - ln(K1/K2)) / (ln(K2/K3) + ln(K1/K2))
%
% 此處 K1≥K2≥K3 為 Eg 特徵值（應變主軸）
%
% 輸入：groupEg.mat（由 AMS_StrainPrep.m 產生）
% 輸出：T-Pj 散布圖（依 Group A/B/C/D 分色），可選擇儲存 PNG
% =========================================================

clearvars; clc;

%% -------- 1. 載入資料 ----------------------------------------
[matFile, matPath] = uigetfile('*.mat', '選擇 groupEg.mat');
if isequal(matFile, 0), error('未選擇檔案，程式結束。'); end
load(fullfile(matPath, matFile), 'groupEg', 'groupNums', 'results');

nGroups = size(groupEg, 3);
fprintf('已載入 %d 組 Eg 張量\n', nGroups);

%% -------- 2. 計算各組 T 與 Pj --------------------------------
T_vals  = nan(nGroups, 1);
Pj_vals = nan(nGroups, 1);
e1_vals = nan(nGroups, 1);
e2_vals = nan(nGroups, 1);
e3_vals = nan(nGroups, 1);

for g = 1:nGroups
    Eg  = groupEg(:,:,g);
    res = AMSFormulas.performEigenAnalysis(Eg);

    e1 = res.K1_val;   % 最大主應變
    e2 = res.K2_val;
    e3 = res.K3_val;   % 最小主應變

    e1_vals(g) = e1;
    e2_vals(g) = e2;
    e3_vals(g) = e3;

    % 只計算三個特徵值均為正的情況
    if e1 <= 0 || e2 <= 0 || e3 <= 0
        warning('Group %d 有非正特徵值，跳過', groupNums(g));
        continue;
    end

    % Km：幾何平均
    Km = (e1 + e2 + e3)/3;

    % Pj（公式 4.8）
    Pj_vals(g) = exp(sqrt(2 * ( log(e1/Km)^2 + log(e2/Km)^2 + log(e3/Km)^2 )));

    % T（公式 4.7）
    lnR1 = log(e1/e2);   % ln(K1/K2)
    lnR2 = log(e2/e3);   % ln(K2/K3)
    denom = lnR2 + lnR1;
    if abs(denom) > 1e-12
        T_vals(g) = (lnR2 - lnR1) / denom;
    else
        T_vals(g) = 0;   % 完全等軸：T=0
    end
end

%% -------- 3. 取得 Group 標籤 ---------------------------------
if isfield(results, 'GroupLabel')
    groupLabels = results.GroupLabel;
else
    groupLabels = repmat({''}, nGroups, 1);
end

% 整理所有出現的 Group，決定顏色
allGroups = unique(groupLabels, 'stable');
allGroups = allGroups(~cellfun(@isempty, allGroups));

% 修改後
colorMap = [ ...
    0.84 0.15 0.16;   % 紅   A
    0.17 0.63 0.17;   % 綠   B
    0.12 0.47 0.71;   % 藍   C
    0.93 0.69 0.13;   % 黃   D
    0.58 0.40 0.74;   % 紫   E
    0.09 0.75 0.81;   % 青   F
    0.77 0.43 0.20;   % 棕   G
    0.50 0.50 0.50];  % 灰   H

fprintf('\n%-8s %-8s %-10s %-10s %-10s %-8s %-8s\n', ...
    'Number','Group','e1','e2','e3','T','Pj');
fprintf('%s\n', repmat('-',1,62));

%% -------- 4. 繪圖 --------------------------------------------
fig = figure('Name','T-Pj Diagram (Jelinek 1981)', ...
    'NumberTitle','off', 'Color','white', 'Position',[100 100 680 620]);
ax = axes(fig);
hold(ax, 'on');

% 背景分區：先用大值佔位，畫完後用 xlim 裁切
xMax = max(1.05, max(Pj_vals(~isnan(Pj_vals)))*1.08);
patch(ax, [1 1 xMax xMax], [0 1 1 0], ...
    [0.94 0.97 1.00], 'EdgeColor','none', 'FaceAlpha',0.5);   % 上半：Oblate
patch(ax, [1 1 xMax xMax], [0 -1 -1 0], ...
    [1.00 0.96 0.94], 'EdgeColor','none', 'FaceAlpha',0.5);   % 下半：Prolate

% T=0 中線
yline(ax, 0, '-', 'Color',[0.6 0.6 0.6], 'LineWidth', 0.8);

% 散布點
legendHandles = gobjects(0);
legendLabels  = {};

if isempty(allGroups)
    % 無 Group 資訊：全部畫同一顏色
    valid = ~isnan(T_vals) & ~isnan(Pj_vals);
    h = scatter(ax, Pj_vals(valid), T_vals(valid), 60, ...
        'filled', 'MarkerFaceColor', colorMap(1,:), ...
        'MarkerEdgeColor', colorMap(1,:)*0.6, 'LineWidth', 0.8);
    legendHandles(end+1) = h;
    legendLabels{end+1}  = 'All';
    printRows(groupNums(valid), groupLabels(valid), ...
        e1_vals(valid), e2_vals(valid), e3_vals(valid), ...
        T_vals(valid), Pj_vals(valid));
else
    for gi = 1:numel(allGroups)
        grp  = allGroups{gi};
        mask = strcmp(groupLabels, grp) & ~isnan(T_vals) & ~isnan(Pj_vals);
        if ~any(mask), continue; end

        cIdx  = mod(gi-1, size(colorMap,1)) + 1;
        color = colorMap(cIdx,:);

        h = scatter(ax, Pj_vals(mask), T_vals(mask), 65, ...
            'filled', 'MarkerFaceColor', color, ...
            'MarkerEdgeColor', color*0.6, 'LineWidth', 0.8, ...
            'DisplayName', sprintf('Group %s', grp));
        legendHandles(end+1) = h; %#ok<AGROW>
        legendLabels{end+1}  = sprintf('Group %s', grp); %#ok<AGROW>

        % Number 標籤
        nums_g = groupNums(mask);
        Pj_g   = Pj_vals(mask);
        T_g    = T_vals(mask);
        for k = 1:numel(nums_g)
            text(ax, Pj_g(k)+0.002, T_g(k)+0.015, ...
                num2str(nums_g(k)), ...
                'FontSize', 7, 'Color', color*0.75, ...
                'HorizontalAlignment','left');
        end

        printRows(groupNums(mask), groupLabels(mask), ...
            e1_vals(mask), e2_vals(mask), e3_vals(mask), ...
            T_vals(mask), Pj_vals(mask));
    end
end

%% -------- 5. 軸設定與標註 ------------------------------------
% X 軸從 1 開始（Pj≥1 恆成立）
xlim(ax, [1, max(1.05, max(Pj_vals(~isnan(Pj_vals)))*1.08)]);
ylim(ax, [-1.05, 1.05]);

% 分區文字
text(ax, 1.005, 0.92, 'Oblate  (0 < T < 1)', ...
    'FontSize', 9, 'Color',[0.30 0.45 0.70], 'FontAngle','italic');
text(ax, 1.005, -0.92, 'Prolate  (-1 < T < 0)', ...
    'FontSize', 9, 'Color',[0.70 0.30 0.25], 'FontAngle','italic');

xlabel(ax, 'P_j  (Corrected degree of anisotropy)', 'FontSize', 12);
ylabel(ax, 'T  (Shape parameter)', 'FontSize', 12);
title(ax, 'T – P_j Diagram  (Jelinek 1981)', 'FontSize', 14, 'FontWeight','bold');
grid(ax, 'on'); box(ax, 'on');

lgd = legend(ax, legendHandles, legendLabels, ...
    'Location','best', 'FontSize', 10);

hold(ax, 'off');

%% -------- 6. 互動式刪點 --------------------------------------
% 方式一：點選圖上的點
% 方式二：輸入 Number 清單
fprintf('\n========== 互動式刪點 ==========\n');
fprintf('目前共 %d 個有效點\n', sum(~isnan(T_vals) & ~isnan(Pj_vals)));

% 記錄哪些點要隱藏（不重新計算，只是從圖上移除）
excludeNums = [];

keepGoing = true;
while keepGoing
    delCh = questdlg('刪除方式？', '互動式刪點', ...
        '點選圖上的點', '輸入 Number', '完成', '完成');
    if isempty(delCh) || strcmp(delCh, '完成'), break; end

    switch delCh

        % ---- 方式一：點選 ----
        case '點選圖上的點'
            fprintf('請在圖上點選要刪除的點（點完按 Enter 結束）\n');
            figure(fig);   % 確保焦點在圖上
            while true
                try
                    [xClick, yClick, btn] = ginput(1);
                catch
                    break;
                end
                if isempty(btn) || btn == 13   % Enter 結束
                    break;
                end

                % 找最近的點（在投影座標系內）
                valid = ~isnan(T_vals) & ~isnan(Pj_vals) & ...
                        ~ismember(groupNums, excludeNums);
                if ~any(valid), break; end

                % 正規化距離（X/Y 軸範圍不同，需各自縮放）
                xl = xlim(ax); yl = ylim(ax);
                dx = (Pj_vals(valid) - xClick) / (xl(2)-xl(1));
                dy = (T_vals(valid)  - yClick) / (yl(2)-yl(1));
                dist = sqrt(dx.^2 + dy.^2);

                validIdx = find(valid);
                [minDist, minI] = min(dist);

                if minDist < 0.05   % 點擊範圍閾值（圖寬的 5%）
                    hitNum = groupNums(validIdx(minI));
                    excludeNums(end+1) = hitNum; %#ok<AGROW>
                    fprintf('  已移除 No.%d  (T=%.4f, Pj=%.4f)\n', ...
                        hitNum, T_vals(validIdx(minI)), Pj_vals(validIdx(minI)));
                else
                    fprintf('  未點中任何點，請再試一次\n');
                end
                refreshPlot(ax, allGroups, groupLabels, groupNums, T_vals, Pj_vals, excludeNums, colorMap);
            end

        % ---- 方式二：輸入 Number ----
        case '輸入 Number'
            ans_str = inputdlg( ...
                {'輸入要刪除的 Number（空格分隔，例如：3 7 15）：'}, ...
                '刪除指定點', 1, {''});
            if isempty(ans_str) || isempty(strtrim(ans_str{1}))
                continue;
            end
            newNums = str2num(ans_str{1}); %#ok<ST2NM>
            for n = newNums
                if ismember(n, groupNums)
                    excludeNums(end+1) = n; %#ok<AGROW>
                    fprintf('  已移除 No.%d\n', n);
                else
                    fprintf('  警告：No.%d 不存在，跳過\n', n);
                end
            end
            refreshPlot(ax, allGroups, groupLabels, groupNums, T_vals, Pj_vals, excludeNums, colorMap);
    end

    % 問是否繼續
    cont = questdlg( ...
        sprintf('已排除 %d 個點，繼續刪點？', numel(excludeNums)), ...
        '繼續？', '繼續刪點', '完成', '復原全部', '完成');
    if strcmp(cont, '復原全部')
        excludeNums = [];
        refreshPlot(ax, allGroups, groupLabels, groupNums, T_vals, Pj_vals, excludeNums, colorMap);
        fprintf('  已復原所有刪除\n');
    elseif strcmp(cont, '完成')
        keepGoing = false;
    end
end

fprintf('最終排除 Number：%s\n', num2str(excludeNums));
fprintf('剩餘有效點：%d 個\n', ...
    sum(~isnan(T_vals) & ~isnan(Pj_vals) & ~ismember(groupNums, excludeNums)));

%% -------- 7. 儲存 --------------------------------------------
saveAns = questdlg('儲存圖片？', 'T-Pj 輸出', '儲存 PNG', '不儲存', '儲存 PNG');
if strcmp(saveAns, '儲存 PNG')
    outFile = fullfile(matPath, 'TPj_Diagram.png');
    exportgraphics(fig, outFile, 'Resolution', 300);
    fprintf('\n圖片已儲存至：%s\n', outFile);
end

disp('plotTPj 執行完畢。');


%% ================================================================
%  本地輔助函式
%% ================================================================
function printRows(nums, grps, e1, e2, e3, T, Pj)
    for i = 1:numel(nums)
        fprintf('%-8d %-8s %-10.5f %-10.5f %-10.5f %-8.4f %-8.4f\n', ...
            nums(i), grps{i}, e1(i), e2(i), e3(i), T(i), Pj(i));
    end
end

function refreshPlot(ax, allGroups, groupLabels, groupNums, ...
                     T_vals, Pj_vals, excludeNums, colorMap)
% 依據剩餘點的數值，自動動態重新縮放 X 軸，並修正背景與圖例

    % 1. 清除舊的點、標籤、背景區塊與圖例，避免疊加與名稱錯亂
    delete(findobj(ax, 'Type','scatter'));
    delete(findobj(ax, 'Type','text'));
    delete(findobj(ax, 'Type','patch'));
    legend(ax, 'off'); % 關閉舊圖例

    % 2. 動態計算「剩餘有效點」的最大 Pj 值，決定新的 X 軸上限
    validMask = ~isnan(T_vals) & ~isnan(Pj_vals) & ~ismember(groupNums, excludeNums);
    if any(validMask)
        currentMaxPj = max(Pj_vals(validMask));
        xMax = max(1.05, currentMaxPj * 1.08); % 留 8% 邊緣裕度
    else
        xMax = 1.05; % 若無剩餘點，預設最小寬度
    end

    % 3. 重新補上背景分區（完美延伸到新的動態 xMax）
    hold(ax, 'on');
    patch(ax, [1 1 xMax xMax], [0 1 1 0], ...
        [0.94 0.97 1.00], 'EdgeColor','none', 'FaceAlpha',0.5);   % 上半：Oblate
    patch(ax, [1 1 xMax xMax], [0 -1 -1 0], ...
        [1.00 0.96 0.94], 'EdgeColor','none', 'FaceAlpha',0.5);   % 下半：Prolate

    % 4. 重新繪製未被排除的散布點，並記錄控制代碼以重建精確圖例
    legendHandles = gobjects(0);
    legendLabels  = {};

    if isempty(allGroups)
        valid = validMask;
        if any(valid)
            h = scatter(ax, Pj_vals(valid), T_vals(valid), 60, ...
                'filled', 'MarkerFaceColor', colorMap(1,:), ...
                'MarkerEdgeColor', colorMap(1,:)*0.6, 'LineWidth', 0.8);
            legendHandles(end+1) = h;
            legendLabels{end+1}  = 'All';
            addNumberLabels(ax, find(valid), groupNums, T_vals, Pj_vals, colorMap(1,:));
        end
    else
        for gi = 1:numel(allGroups)
            grp  = allGroups{gi};
            mask = strcmp(groupLabels, grp) & validMask;
            if ~any(mask), continue; end
            
            cIdx  = mod(gi-1, size(colorMap,1)) + 1;
            color = colorMap(cIdx,:);
            h = scatter(ax, Pj_vals(mask), T_vals(mask), 65, ...
                'filled', 'MarkerFaceColor', color, ...
                'MarkerEdgeColor', color*0.6, 'LineWidth', 0.8, ...
                'DisplayName', sprintf('Group %s', grp));
            
            legendHandles(end+1) = h; %#ok<AGROW>
            legendLabels{end+1}  = sprintf('Group %s', grp); %#ok<AGROW>
            addNumberLabels(ax, find(mask), groupNums, T_vals, Pj_vals, color);
        end
    end
    
    % 5. 關鍵修正：將 X 軸緊貼重新縮放後的大小
    xlim(ax, [1, xMax]);
    ylim(ax, [-1.05, 1.05]);

    % 6. 分區文字與圖例補回
    text(ax, 1.005, 0.92,  'Oblate  (0 < T < 1)', ...
        'FontSize',9,'Color',[0.30 0.45 0.70],'FontAngle','italic');
    text(ax, 1.005, -0.92, 'Prolate  (-1 < T < 0)', ...
        'FontSize',9,'Color',[0.70 0.30 0.25],'FontAngle','italic');
    
    if ~isempty(legendHandles)
        legend(ax, legendHandles, legendLabels, 'Location','best', 'FontSize', 10);
    end
    
    hold(ax, 'off');
    drawnow;
end

function addNumberLabels(ax, idxList, groupNums, T_vals, Pj_vals, color)
    for k = 1:numel(idxList)
        i = idxList(k);
        text(ax, Pj_vals(i)+0.002, T_vals(i)+0.015, ...
            num2str(groupNums(i)), ...
            'FontSize',7, 'Color',color*0.75, 'HorizontalAlignment','left');
    end
end
