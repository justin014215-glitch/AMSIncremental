%% AMS_IncrementalStrain.m
% 應變增量計算主程式
%
% 前置需求：先執行 AMS_StrainPrep.m → 產生 groupEg.mat
%
% 模式 1：手選任意兩組，計算單一增量
% 模式 2：依分區 (A/B/C/D) 自動分組，可選擇計算「單一分區內部連續增量」或「相鄰分區對連續增量」
% 模式 3：從 No.1 連續計算到最後一個樣本（1→2, 2→3 … n-1→n）
%          所有增量主軸畫在同一張 stereonet，並輸出彙整 Excel
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
fprintf('\n%-8s %-6s %-12s %-12s %-12s  %-18s\n', ...
    'Number','Group','Eg_e1','Eg_e2','Eg_e3','Sample_Names');
fprintf('%s\n', repmat('-',1,78));
for g = 1:nGroups
    EgRes = AMSFormulas.performEigenAnalysis(groupEg(:,:,g));
    nm = '';
    if isfield(results,'SampleNames') && g<=numel(results.SampleNames)
        nm = results.SampleNames{g};
    end
    gl = '';
    if isfield(results,'GroupLabel') && g<=numel(results.GroupLabel)
        gl = results.GroupLabel{g};
    end
    fprintf('%-8d %-6s %-12.6f %-12.6f %-12.6f  %s\n', ...
        groupNums(g), gl, EgRes.K1_val, EgRes.K2_val, EgRes.K3_val, nm);
end
fprintf('\n');

%% -------- 2. 選擇計算模式 ----------------------------------------
modeChoice = questdlg('選擇計算模式','應變增量模式', ...
    '模式1：選兩個樣本', ...
    '模式2：依分區(A/B/C/D)', ...
    '模式3：連續增量(1→2→…→n)', ...
    '模式3：連續增量(1→2→…→n)');

if isempty(modeChoice), disp('已取消'); return; end

%% ================================================================
%  模式 1：手選兩組
%% ================================================================
if strcmp(modeChoice,'模式1：選兩個樣本')

    keepGoing = true;
    allPairResults = {};

    while keepGoing
        fprintf('\n--- 配對 #%d ---\n', numel(allPairResults)+1);
        numList = arrayfun(@(x) num2str(x), groupNums, 'UniformOutput',false);

        [idx1,ok] = listdlg('PromptString','選擇初始狀態 Ef1：', ...
            'SelectionMode','single','ListString',numList,'Name','選擇 Ef1');
        if ~ok, break; end
        [idx2,ok] = listdlg('PromptString','選擇最終狀態 Ef2：', ...
            'SelectionMode','single','ListString',numList,'Name','選擇 Ef2');
        if ~ok, break; end

        pairRes = computeIncrement(groupEg(:,:,idx1), groupEg(:,:,idx2), ...
            groupNums(idx1), groupNums(idx2));
        allPairResults{end+1} = pairRes;

        cont = questdlg('是否再計算一組？','繼續？','繼續','結束','結束');
        if strcmp(cont,'結束'), keepGoing=false; end
    end

    if ~isempty(allPairResults)
        outputResults(allPairResults, matPath, 'PairResults');
    end
end

%% ================================================================
%  模式 2：依分區 A/B/C/D — 連續增量（可選單一分區內部 或 跨分區）
%% ================================================================
if strcmp(modeChoice,'模式2：依分區(A/B/C/D)')

    fprintf('\n============================\n');
    fprintf('  模式 2：依分區連續應變增量\n');
    fprintf('============================\n');

    % ---- 從 results.GroupLabel 自動建立分區索引 -----------------
    if ~isfield(results,'GroupLabel')
        errordlg('results 中找不到 GroupLabel 欄位，請重新執行 AMS_StrainPrep','欄位缺失');
        return;
    end

    cleanLabels = cellfun(@(s) normalizeGroupStr(s), results.GroupLabel, ...
        'UniformOutput', false);

    % 找出實際存在的分區（依 A B C D 順序）
    validSections = {};
    for s = {'A','B','C','D'}
        sn = s{1};
        if any(strcmp(cleanLabels, sn))
            validSections{end+1} = sn; 
        end
    end

    if isempty(validSections)
        errordlg('未找到任何有效分區','分區不足');
        return;
    end

    % 建立每個分區的站點清單（依 Number 升序排列）
    sectionNums = struct();
    fprintf('\n各分區站點（依 Number 排序）：\n');
    for s = 1:numel(validSections)
        sn   = validSections{s};
        mask = strcmp(cleanLabels, sn);
        nums = sort(groupNums(mask), 'ascend');
        sectionNums.(sn) = nums;
        fprintf('  分區 %s：No.%s（共 %d 站）\n', sn, ...
            strjoin(arrayfun(@num2str,nums,'UniformOutput',false),', '), numel(nums));
    end

    % ---- 動態選單組合：同時提供「單一分區內部」與「跨分區」選項 -----
    menuOptions = {};
    
    % 1. 先加入單一分區內部連續增量選項
    for s = 1:numel(validSections)
        sn = validSections{s};
        if numel(sectionNums.(sn)) >= 2
            menuOptions{end+1} = sprintf('僅分區 %s 內部連續', sn); %#ok<AGROW>
        end
    end
    
    % 2. 再加入跨相鄰分區的選項
    for s = 1:numel(validSections)-1
        menuOptions{end+1} = sprintf('%s → %s', ...
            validSections{s}, validSections{s+1}); %#ok<AGROW>
    end

    [pairSel, ok] = listdlg( ...
        'PromptString', '選擇要計算的連續類型（可多選）：', ...
        'SelectionMode', 'multiple', ...
        'ListString',    menuOptions, ...
        'Name',          '選擇計算分區', ...
        'ListSize',      [280 220]);
    if ~ok, disp('已取消'); return; end

    selectedOptions = menuOptions(pairSel);

    % ---- 開始針對選定的選項進行計算與合圖 -------------------------
    for so = 1:numel(selectedOptions)
        currentOpt = selectedOptions{so};
        
        if startsWith(currentOpt, '僅分區')
            %% --- 狀況 A：只看單個分區內部的疊加 ---
            % 提取分區字母，例如從 '僅分區 A 內部連續' 拿掉字串取得 'A'
            sA = extractBetween(currentOpt, "僅分區 ", " 內部連續");
            sA = sA{1};
            sB = ''; % 單一分區沒有終點分區
            
            seqNums = sectionNums.(sA);
            nSeq    = numel(seqNums);
            tag     = sprintf('Sec_%s_Internal', sA);
            
            fprintf('\n====== 僅分區 %s 內部連續（共 %d 個增量）======\n', sA, nSeq-1);
            
        else
            %% --- 狀況 B：原本的 A → B 跨分區連續 ---
            parts = strsplit(currentOpt, ' → ');
            sA = parts{1};   
            sB = parts{2};   

            numsA = sectionNums.(sA);
            numsB = sectionNums.(sB);

            % 拼接 A+B
            seqNums = [numsA(:); numsB(:)];
            nSeq    = numel(seqNums);
            tag     = sprintf('Sec_%s_to_%s', sA, sB);
            
            fprintf('\n====== 分區對 %s → %s（共 %d 個增量）======\n', sA, sB, nSeq-1);
        end

        if nSeq < 2
            fprintf('站點不足，跳過：%s\n', currentOpt);
            continue;
        end

        % 逐對計算增量
        seqResults = cell(nSeq-1, 1);
        for p = 1:nSeq-1
            n1   = seqNums(p);
            n2   = seqNums(p+1);
            idx1 = find(groupNums == n1, 1);
            idx2 = find(groupNums == n2, 1);
            if isempty(idx1)||isempty(idx2)
                continue;
            end

            % 計算增量
            seqResults{p} = computeIncrement(groupEg(:,:,idx1), groupEg(:,:,idx2), n1, n2);
            
            % 判斷是否為跨分區（只有在兩者皆存在且標籤不同時才是跨分區）
            g1 = cleanLabels{groupNums==n1};
            g2 = cleanLabels{groupNums==n2};
            if ~strcmp(g1, g2)
                seqResults{p}.label   = [seqResults{p}.label ' *'];
                seqResults{p}.isCross = true;
            else
                seqResults{p}.isCross = false;
            end
        end

        % 移除空資料並輸出
        seqResults = seqResults(~cellfun(@isempty, seqResults));
        if ~isempty(seqResults)
            outputResults(seqResults, matPath, tag, sA, sB);
        end
    end
end

%% ================================================================
%  模式 3：連續增量 1→2→3→…→n
%% ================================================================
if strcmp(modeChoice,'模式3：連續增量(1→2→…→n)')

    fprintf('\n============================\n');
    fprintf('  模式 3：連續應變增量\n');
    fprintf('============================\n\n');

    if nGroups < 2
        errordlg('至少需要 2 組資料','資料不足'); return;
    end

    [sortedNums, sortIdx] = sort(groupNums, 'ascend');
    nPairs = nGroups - 1;

    fprintf('計算順序：');
    for g = 1:nGroups
        if g < nGroups, fprintf('No.%d → ', sortedNums(g));
        else,           fprintf('No.%d\n', sortedNums(g));
        end
    end
    fprintf('\n');

    seqResults = cell(nPairs,1);
    for p = 1:nPairs
        g1 = sortIdx(p);
        g2 = sortIdx(p+1);
        seqResults{p} = computeIncrement( ...
            groupEg(:,:,g1), groupEg(:,:,g2), ...
            groupNums(g1), groupNums(g2));
    end

    outputResults(seqResults, matPath, 'SeqResults');
end

disp('程式執行完畢。');


%% ================================================================
%  核心計算
%% ================================================================
function res = computeIncrement(Ef1, Ef2, label1, label2)
    label = sprintf('No.%s → No.%s', num2str(label1), num2str(label2));
    fprintf('\n[%s]\n', label);

    F     = AMSFormulas.calculateF(Ef1, Ef2);
    F_sym = AMSFormulas.symmetrize(F);
    FRes  = AMSFormulas.performEigenAnalysis(F_sym);

    T_val = NaN;
    if all([FRes.K1_val FRes.K2_val FRes.K3_val] > 0)
        lnR1 = log(FRes.K1_val/FRes.K2_val);
        lnR2 = log(FRes.K2_val/FRes.K3_val);
        if (lnR1+lnR2) ~= 0
            T_val = (lnR2-lnR1)/(lnR2+lnR1);
        end
    end

    fprintf('  det(F)=%.6f  e1=%.6f  e2=%.6f  e3=%.6f\n', ...
        det(F), FRes.K1_val, FRes.K2_val, FRes.K3_val);
    fprintf('  e1: %5.1f°/%4.1f°    e3: %5.1f°/%4.1f°\n', ...
        FRes.K1_dir(1),FRes.K1_dir(2), FRes.K3_dir(1),FRes.K3_dir(2));
    
    res.label   = label;
    res.label1  = label1;
    res.label2  = label2;
    res.F       = F;
    res.F_sym   = F_sym;
    res.FRes    = FRes;
    res.T_val   = T_val;
    res.detF    = det(F);
    res.isCross = false;
end


function s = normalizeGroupStr(s)
    s = strtrim(strrep(s, '''', ''));
    for c = 1:numel(s)
        code = double(s(c));
        if code >= 65313 && code <= 65338   
            s(c) = char(code - 65248);
        end
    end
end


%% ================================================================
%  輸出（Excel + 合圖 stereonet）
%% ================================================================
function outputResults(pairResults, matPath, tag, sectionA, sectionB)
    nP = numel(pairResults);

    C1=[0.12 0.47 0.71]; C1e=[0.05 0.25 0.50];
    C2=[0.17 0.63 0.17]; C2e=[0.05 0.38 0.05];
    C3=[0.84 0.15 0.16]; C3e=[0.55 0.05 0.05];

    if nargin < 4, sectionA=''; sectionB=''; end

    %% --- Excel 輸出 -------------------------------------------
    Label     = cellfun(@(r) r.label,            pairResults,'UniformOutput',false)';
    Ef1_ID    = cellfun(@(r) num2str(r.label1), pairResults,'UniformOutput',false)';
    Ef2_ID    = cellfun(@(r) num2str(r.label2), pairResults,'UniformOutput',false)';
    det_F     = cellfun(@(r) r.detF,            pairResults)';
    F_e1      = cellfun(@(r) r.FRes.K1_val,     pairResults)';
    F_e2      = cellfun(@(r) r.FRes.K2_val,     pairResults)';
    F_e3      = cellfun(@(r) r.FRes.K3_val,     pairResults)';
    e1_Trend  = cellfun(@(r) r.FRes.K1_dir(1),  pairResults)';
    e1_Plunge = cellfun(@(r) r.FRes.K1_dir(2),  pairResults)';
    e2_Trend  = cellfun(@(r) r.FRes.K2_dir(1),  pairResults)';
    e2_Plunge = cellfun(@(r) r.FRes.K2_dir(2),  pairResults)';
    e3_Trend  = cellfun(@(r) r.FRes.K3_dir(1),  pairResults)';
    e3_Plunge = cellfun(@(r) r.FRes.K3_dir(2),  pairResults)';
    L_val     = cellfun(@(r) r.FRes.L,          pairResults)';
    Fv_val    = cellfun(@(r) r.FRes.F,          pairResults)';
    T_shape   = cellfun(@(r) r.T_val,           pairResults)';

    T = table(Label,Ef1_ID,Ef2_ID,det_F, ...
        F_e1,F_e2,F_e3, ...
        e1_Trend,e1_Plunge,e2_Trend,e2_Plunge,e3_Trend,e3_Plunge, ...
        L_val,Fv_val,T_shape, ...
        'VariableNames',{ ...
        'Label','Ef1_ID','Ef2_ID','det_F', ...
        'F_e1','F_e2','F_e3', ...
        'e1_Trend','e1_Plunge','e2_Trend','e2_Plunge','e3_Trend','e3_Plunge', ...
        'L_lineation','F_foliation','T_shape'});

    xlsFile = fullfile(matPath, ['IncrementalStrain_' tag '.xlsx']);
    writetable(T, xlsFile, 'Sheet', tag);
    fprintf('\nExcel 輸出：%s\n', xlsFile);
    save(fullfile(matPath,['IncrementalStrain_' tag '.mat']), 'pairResults');

    %% --- 儲存路徑 ---------------------------------------------
    stereoDir = fullfile(matPath, ['IncrementalStereonets_' tag]);
    if ~exist(stereoDir,'dir'), mkdir(stereoDir); end

    %% --- 繪製立體投影網圖 ---------------------------------------
    if ~isempty(sectionA) && ~isempty(sectionB)
        figTitle = sprintf('應變增量  %s → %s（共 %d 個增量）', sectionA, sectionB, nP);
    elseif ~isempty(sectionA) && isempty(sectionB)
        figTitle = sprintf('應變增量  分區 %s 內部連續（共 %d 個增量）', sectionA, nP);
    else
        figTitle = sprintf('應變增量合圖（共 %d 對）', nP);
    end

    fig = figure('Name', figTitle, 'Color','white', 'Position',[150 150 700 720]);
    Stereonet(0, pi/2, 10*pi/180, 1);
    ax = gca; hold(ax,'on');

    for p = 1:nP
        res  = pairResults{p};
        FRes = res.FRes;
        lbl  = res.label;

        isCross = isfield(res,'isCross') && res.isCross;
        if isCross
            faceC1='none'; faceC2='none'; faceC3='none';
            lw = 1.8;
        else
            faceC1=C1; faceC2=C2; faceC3=C3;
            lw = 1.2;
        end

        [x1,y1] = schmidtProject(FRes.K1_dir(1), FRes.K1_dir(2));
        plot(ax,x1,y1,'s','MarkerFaceColor',faceC1,'MarkerEdgeColor',C1e, 'MarkerSize',13,'LineWidth',lw);
        text(ax,x1+0.04,y1+0.03, lbl, 'FontSize',7,'Color',C1e);

        [x2,y2] = schmidtProject(FRes.K2_dir(1), FRes.K2_dir(2));
        plot(ax,x2,y2,'^','MarkerFaceColor',faceC2,'MarkerEdgeColor',C2e, 'MarkerSize',13,'LineWidth',lw);

        [x3,y3] = schmidtProject(FRes.K3_dir(1), FRes.K3_dir(2));
        plot(ax,x3,y3,'o','MarkerFaceColor',faceC3,'MarkerEdgeColor',C3e, 'MarkerSize',13,'LineWidth',lw);
    end

    %% --- 圖例優化 (不寫顏色) -----------------------------------
    h1=plot(ax,NaN,NaN,'s','MarkerFaceColor',C1,'MarkerEdgeColor',C1e, 'MarkerSize',13,'DisplayName','e_1');
    h2=plot(ax,NaN,NaN,'^','MarkerFaceColor',C2,'MarkerEdgeColor',C2e, 'MarkerSize',13,'DisplayName','e_2');
    h3=plot(ax,NaN,NaN,'o','MarkerFaceColor',C3,'MarkerEdgeColor',C3e, 'MarkerSize',13,'DisplayName','e_3');

    hasCross = any(cellfun(@(r) isfield(r,'isCross')&&r.isCross, pairResults));
    if hasCross
        hx=plot(ax,NaN,NaN,'s','MarkerFaceColor','none','MarkerEdgeColor',C1e, 'MarkerSize',13,'LineWidth',1.8,'DisplayName','跨分區增量');
        legend(ax,[h1 h2 h3 hx],'Location','southoutside', 'Orientation','horizontal','FontSize',10);
    else
        legend(ax,[h1 h2 h3],'Location','southoutside', 'Orientation','horizontal','FontSize',10);
    end

    title(ax, figTitle, 'FontSize',12,'FontWeight','bold');
    annotation('textbox',[0.12 0.01 0.76 0.04], ...
        'String','Lower hemisphere, equal-area projection (Schmidt Net)', ...
        'EdgeColor','none','FontSize',7,'Color',[0.5 0.5 0.5], 'HorizontalAlignment','center');
    axis(ax,'off'); hold(ax,'off');

    fname = sprintf('IncrementAll_%s.png', tag);
    exportgraphics(fig, fullfile(stereoDir,fname), 'Resolution',300);
    fprintf('合圖 saved → %s\n', fullfile(stereoDir,fname));
end


%% ================================================================
%  幾何工具
%% ================================================================
function [x,y] = schmidtProject(trend_deg, plunge_deg)
    r = sqrt(2) * sin((pi/2 - deg2rad(plunge_deg)) / 2);
    x = r * sin(deg2rad(trend_deg));
    y = r * cos(deg2rad(trend_deg));
end