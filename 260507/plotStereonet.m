function plotStereonet(results, rawData, varargin)
% plotStereonet  AMS 主軸下半球等面積投影（Schmidt Net）
%
% 必要輸入：
%   results  - AMS_StrainPrep 輸出的 struct（含平均主軸方向與 Jelinek 誤差）
%   rawData  - 原始資料 table（含各方磚 K1/K2/K3 與 trend/plunge）
%
% 選用參數（直接帶參、跳過選單）：
%   'DataType'    - 'mean' | 'raw'
%   'GroupLabels' - {'A','B'} 等
%   'Numbers'     - [1 3 5] 等
%   'ShowLabels'  - true/false（預設 false）
%   'EllipseType' - 'bootstrap'|'jelinek'|'both'|'none'（預設 'none'）
%   'BootstrapN'  - Bootstrap 抽樣次數（預設 1000）
%   'SaveFig'     - true/false（預設 true）
%   'SavePath'    - 儲存路徑（預設 pwd）
%
% 視覺規則：
%   K1 □ 藍，K2 △ 綠，K3 ○ 紅
%   個別方磚：小符號半透明
%   平均主軸：大符號不透明
%   Bootstrap 橢圓：各軸同色實線（K1藍、K2綠、K3紅）
%   Jelinek  橢圓：各軸同色虛線（K1藍、K2綠、K3紅）
%
% 修改說明（2025）：
%   - Raw 模式下，三個主軸（K1/K2/K3）各自畫信心橢圓
%   - Bootstrap：PCA 擬合長短軸橢圓，三軸各自繪製
%   - Jelinek：以各主軸為中心、沿其與另兩軸所在平面方向展開的真正橢圓
%              K1 橢圓：長軸=E12（K1-K2方向）、短軸=E13（K1-K3方向）
%              K2 橢圓：長軸=E12（K2-K1方向）、短軸=E23（K2-K3方向）
%              K3 橢圓：長軸=E13（K3-K1方向）、短軸=E23（K3-K2方向）
% =========================================================

%% 解析參數
p = inputParser;
addParameter(p, 'DataType',    '',       @ischar);
addParameter(p, 'GroupLabels', {},       @iscell);
addParameter(p, 'Numbers',     [],       @isnumeric);
addParameter(p, 'ShowLabels',  false,    @islogical);
addParameter(p, 'EllipseType', 'none',   @ischar);
addParameter(p, 'BootstrapN',  1000,     @isnumeric);
addParameter(p, 'SaveFig',     true,     @islogical);
addParameter(p, 'SavePath',    pwd,      @ischar);
parse(p, varargin{:});
opt = p.Results;

% 顏色常數
C1=[0.12 0.47 0.71]; C1e=[0.05 0.25 0.50];  % K1 藍
C2=[0.17 0.63 0.17]; C2e=[0.05 0.38 0.05];  % K2 綠
C3=[0.84 0.15 0.16]; C3e=[0.55 0.05 0.05];  % K3 紅

%% 主選單
if isempty(opt.DataType)
    mainCh = questdlg('繪圖方式？','AMS Stereonet', ...
        '畫全部','互動選擇','取消','互動選擇');
    if isempty(mainCh) || strcmp(mainCh,'取消'), return; end

    if strcmp(mainCh,'互動選擇')
        opt = runInteractiveMenu(results, rawData, opt);
        if isempty(opt), return; end
    else
        opt = askAllMode(opt);
        if isempty(opt), return; end
    end
end

%% 驗證 DataType
if ~ismember(opt.DataType,{'mean','raw'})
    error('DataType 必須是 ''mean'' 或 ''raw''');
end

%% ============================================================
%  畫平均主軸（mean 模式：只畫符號，不畫橢圓）
%% ============================================================
if strcmp(opt.DataType,'mean')

    idx = filterMean(results, opt);
    if isempty(idx), errordlg('找不到符合條件的站點','錯誤'); return; end

    titleStr = makeTitle('mean', opt);
    createFigure(titleStr);
    ax = initAxes();

    for ii = 1:numel(idx)
        g = idx(ii);
        drawLarge(ax, ...
            results.mean_K1_trend(g), results.mean_K1_plunge(g), ...
            results.mean_K2_trend(g), results.mean_K2_plunge(g), ...
            results.mean_K3_trend(g), results.mean_K3_plunge(g), ...
            C1,C1e,C2,C2e,C3,C3e);

        if opt.ShowLabels
            lbl = meanLabel(results, g);
            addLabel(ax, results.mean_K1_trend(g), results.mean_K1_plunge(g), lbl, C1, 7);
            addLabel(ax, results.mean_K2_trend(g), results.mean_K2_plunge(g), lbl, C2, 7);
            addLabel(ax, results.mean_K3_trend(g), results.mean_K3_plunge(g), lbl, C3, 7);
        end
    end

    closeFigure(ax, titleStr, 'mean', opt, C1,C1e,C2,C2e,C3,C3e);

%% ============================================================
%  畫原始方磚（raw 模式：個別符號 + 平均符號 + 三軸橢圓）
%% ============================================================
else

    rowIdx = filterRaw(rawData, opt);
    if isempty(rowIdx), errordlg('找不到符合條件的方磚','錯誤'); return; end

    % 讀取欄位
    Numbers = rawData{:,'Number'};
    dK1=rawData{:,'dK1geo'}; iK1=rawData{:,'iK1geo'};
    dK2=rawData{:,'dK2geo'}; iK2=rawData{:,'iK2geo'};
    dK3=rawData{:,'dK3geo'}; iK3=rawData{:,'iK3geo'};
    K1v=rawData{:,'K1'}; K2v=rawData{:,'K2'}; K3v=rawData{:,'K3'};

    titleStr = makeTitle('raw', opt);
    createFigure(titleStr);
    ax = initAxes();

    numsInPlot = unique(Numbers(rowIdx),'stable');

    for ni = 1:numel(numsInPlot)
        thisNum = numsInPlot(ni);
        subIdx  = rowIdx(Numbers(rowIdx)==thisNum);
        nSub    = numel(subIdx);

        % (1) 個別方磚（小半透明符號）
        for ii = 1:nSub
            r = subIdx(ii);
            drawSmall(ax, dK1(r),iK1(r), dK2(r),iK2(r), dK3(r),iK3(r), C1,C2,C3);
            if opt.ShowLabels
                lbl = sprintf('No.%d', thisNum);
                addLabel(ax, dK1(r),iK1(r), lbl, C1*0.65, 5);
                addLabel(ax, dK3(r),iK3(r), lbl, C3*0.65, 5);
            end
        end

        % (2) 建構 tensorStack → 平均主軸
        tensorStack = zeros(3,3,nSub);
        for ii = 1:nSub
            r = subIdx(ii);
            V = AMSFormulas.calculateV(dK1(r),iK1(r),dK2(r),iK2(r),dK3(r),iK3(r));
            Kg = V * diag([K1v(r),K2v(r),K3v(r)]) * V';
            tensorStack(:,:,ii) = AMSFormulas.symmetrize(Kg);
        end
        M_mean = AMSFormulas.calculateMeanTensor(tensorStack);
        eigRes = AMSFormulas.performEigenAnalysis(M_mean);

        % (3) 平均主軸（大符號）
        drawLarge(ax, ...
            eigRes.K1_dir(1),eigRes.K1_dir(2), ...
            eigRes.K2_dir(1),eigRes.K2_dir(2), ...
            eigRes.K3_dir(1),eigRes.K3_dir(2), ...
            C1,C1e,C2,C2e,C3,C3e);
        if opt.ShowLabels
            addLabel(ax, eigRes.K1_dir(1),eigRes.K1_dir(2), sprintf('No.%d',thisNum), C1, 7);
            addLabel(ax, eigRes.K2_dir(1),eigRes.K2_dir(2), sprintf('No.%d',thisNum), C2, 7);
            addLabel(ax, eigRes.K3_dir(1),eigRes.K3_dir(2), sprintf('No.%d',thisNum), C3, 7);
        end

        % (4) 信心橢圓（三軸各自繪製）
        if ~strcmp(opt.EllipseType,'none') && nSub >= 3
            doB = ismember(opt.EllipseType,{'bootstrap','both'});
            doJ = ismember(opt.EllipseType,{'jelinek','both'});

            % --- Bootstrap：K1/K2/K3 三軸各自 PCA 橢圓（同色實線）---
            if doB
                bc = calcBootstrap(tensorStack, opt.BootstrapN);
                drawBootEllipse(ax, eigRes.K1_dir, eigRes.K2_dir, eigRes.K3_dir, bc, 'K1', C1);
                drawBootEllipse(ax, eigRes.K2_dir, eigRes.K1_dir, eigRes.K3_dir, bc, 'K2', C2);
                drawBootEllipse(ax, eigRes.K3_dir, eigRes.K1_dir, eigRes.K2_dir, bc, 'K3', C3);
                fprintf('  No.%d  Bootstrap 95%%  K1=%.1f°×%.1f°  K2=%.1f°×%.1f°  K3=%.1f°×%.1f°\n', ...
                    thisNum, ...
                    bc.K1_halfAngle, bc.K1_minorAngle, ...
                    bc.K2_halfAngle, bc.K2_minorAngle, ...
                    bc.K3_halfAngle, bc.K3_minorAngle);
            end

            % --- Jelinek：K1/K2/K3 三軸橢圓（同色虛線）---
            % K1 橢圓：長軸方向→K2，半角=E12；短軸方向→K3，半角=E13
            % K2 橢圓：長軸方向→K1，半角=E12；短軸方向→K3，半角=E23
            % K3 橢圓：長軸方向→K1，半角=E13；短軸方向→K2，半角=E23
            if doJ
                js = AMSFormulas.calculateJelinekStats(tensorStack);
                drawJelEllipse(ax, eigRes.K1_dir, eigRes.K2_dir, eigRes.K3_dir, js.E12, js.E13, C1);
                drawJelEllipse(ax, eigRes.K2_dir, eigRes.K1_dir, eigRes.K3_dir, js.E12, js.E23, C2);
                drawJelEllipse(ax, eigRes.K3_dir, eigRes.K1_dir, eigRes.K2_dir, js.E13, js.E23, C3);
                fprintf('  No.%d  Jelinek  E12=%.1f°  E13=%.1f°  E23=%.1f°\n', ...
                    thisNum, js.E12, js.E13, js.E23);
            end
        end
    end

    closeFigure(ax, titleStr, 'raw', opt, C1,C1e,C2,C2e,C3,C3e);
end

end % ==================== END MAIN ====================


%% ================================================================
%  選單
%% ================================================================
function opt = runInteractiveMenu(results, rawData, opt)
% 互動式選單

    % 步驟1：資料類型
    dtCh = questdlg('要畫哪種資料？','資料類型', ...
        '平均主軸','原始方磚（含誤差橢圓）','取消','原始方磚（含誤差橢圓）');
    if isempty(dtCh)||strcmp(dtCh,'取消'), opt=[]; return; end
    if strcmp(dtCh,'平均主軸'), opt.DataType='mean';
    else,                       opt.DataType='raw';
    end

    % 步驟2：篩選範圍
    fCh = questdlg('篩選範圍？','篩選', ...
        '依 Group 分區','依站點 Number','全部','依站點 Number');
    if isempty(fCh), opt=[]; return; end

    switch fCh
        case '依 Group 分區'
            if strcmp(opt.DataType,'mean')
                allGL = unique(normGroup(results.GroupLabel),'stable');
            else
                allGL = unique(normGroup(rawData{:,'Group'}),'stable');
            end
            allGL = allGL(~cellfun(@isempty,allGL));
            if isempty(allGL)
                errordlg('找不到 Group 資料','錯誤'); opt=[]; return;
            end
            [sel,ok] = listdlg('PromptString','選擇分區（可多選）：', ...
                'SelectionMode','multiple','ListString',allGL, ...
                'Name','Group 分區','ListSize',[220 200]);
            if ~ok, opt=[]; return; end
            opt.GroupLabels = allGL(sel);
            opt.Numbers     = [];

        case '依站點 Number'
            if strcmp(opt.DataType,'mean')
                allNums = results.Number;
                listStr = arrayfun(@(i) sprintf('[%s] No.%d', ...
                    results.GroupLabel{i}, allNums(i)), ...
                    (1:numel(allNums))','UniformOutput',false);
            else
                allNums = unique(rawData{:,'Number'},'stable');
                cG = normGroup(rawData{:,'Group'});
                listStr = arrayfun(@(n) sprintf('[%s] No.%d', ...
                    cG{find(rawData{:,'Number'}==n,1)}, n), ...
                    allNums,'UniformOutput',false);
            end
            [sel,ok] = listdlg('PromptString','選擇站點（可多選）：', ...
                'SelectionMode','multiple','ListString',listStr, ...
                'Name','站點 Number','ListSize',[280 380]);
            if ~ok, opt=[]; return; end
            opt.Numbers     = allNums(sel);
            opt.GroupLabels = {};

        otherwise  % 全部
            opt.GroupLabels = {};
            opt.Numbers     = [];
    end

    % 步驟3：標籤
    lCh = questdlg('顯示站點標籤？','標籤','顯示','不顯示','不顯示');
    opt.ShowLabels = strcmp(lCh,'顯示');

    % 步驟4：信心橢圓（原始模式才問）
    if strcmp(opt.DataType,'raw')
        eCh = questdlg('信心橢圓類型？（三軸各自繪製）','信心橢圓', ...
            'Bootstrap（實線）','Jelinek（虛線）','兩者都畫','兩者都畫');
        switch eCh
            case 'Bootstrap（實線）', opt.EllipseType='bootstrap';
            case 'Jelinek（虛線）',   opt.EllipseType='jelinek';
            case '兩者都畫',          opt.EllipseType='both';
            otherwise
                nCh = questdlg('','信心橢圓','不顯示','取消','不顯示');
                if strcmp(nCh,'不顯示'), opt.EllipseType='none';
                else, opt=[]; return;
                end
        end
    else
        opt.EllipseType = 'none';
    end
end

function opt = askAllMode(opt)
% 「畫全部」時問資料類型與橢圓類型
    dtCh = questdlg('畫全部：資料類型？','畫全部', ...
        '平均主軸','原始方磚（含誤差橢圓）','取消','原始方磚（含誤差橢圓）');
    if isempty(dtCh)||strcmp(dtCh,'取消'), opt=[]; return; end
    if strcmp(dtCh,'平均主軸'), opt.DataType='mean';
    else,                       opt.DataType='raw';
    end
    opt.GroupLabels = {};
    opt.Numbers     = [];

    if strcmp(opt.DataType,'raw')
        eCh = questdlg('信心橢圓類型？（三軸各自繪製）','信心橢圓', ...
            'Bootstrap（實線）','Jelinek（虛線）','兩者都畫','兩者都畫');
        switch eCh
            case 'Bootstrap（實線）', opt.EllipseType='bootstrap';
            case 'Jelinek（虛線）',   opt.EllipseType='jelinek';
            case '兩者都畫',          opt.EllipseType='both';
            otherwise,                opt.EllipseType='none';
        end
    else
        opt.EllipseType='none';
    end
end


%% ================================================================
%  篩選
%% ================================================================
function idx = filterMean(results, opt)
    if ~isempty(opt.GroupLabels)
        idx = find(ismember(normGroup(results.GroupLabel), normGroup(opt.GroupLabels)));
    elseif ~isempty(opt.Numbers)
        idx = find(ismember(results.Number, opt.Numbers));
    else
        idx = 1:numel(results.Number);
    end
end

function rowIdx = filterRaw(rawData, opt)
    Numbers = rawData{:,'Number'};
    cleanG  = normGroup(rawData{:,'Group'});
    if ~isempty(opt.GroupLabels)
        rowIdx = find(ismember(cleanG, normGroup(opt.GroupLabels)));
    elseif ~isempty(opt.Numbers)
        rowIdx = find(ismember(Numbers, opt.Numbers));
    else
        rowIdx = (1:numel(Numbers))';
    end
end


%% ================================================================
%  圖形建立與收尾
%% ================================================================
function createFigure(titleStr)
    figure('Name',titleStr,'NumberTitle','off', ...
        'Position',[100 100 700 740],'Color','white');
end

function ax = initAxes()
% 優先使用外部 Stereonet()，若不存在則 fallback 到內建底圖
    if exist('Stereonet','file') || exist('Stereonet','builtin')
        try
            Stereonet(0, pi/2, 10*pi/180, 1);
            ax = gca;
            hold(ax, 'on');
            return;
        catch
            % Stereonet() 存在但呼叫失敗，改用內建
        end
    end
    % ---- 內建 Schmidt Net 底圖 ----
    ax = gca;
    hold(ax, 'on');
    axis(ax, 'equal'); axis(ax, 'off');
    th = linspace(0, 2*pi, 361);
    plot(ax, cos(th), sin(th), 'k-', 'LineWidth', 1.5);
    plot(ax, [-1 1], [0 0], 'k-', 'LineWidth', 0.5);
    plot(ax, [0 0], [-1 1], 'k-', 'LineWidth', 0.5);
    text(ax,  0,    1.07, 'N','HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
    text(ax,  0,   -1.09, 'S','HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
    text(ax,  1.08, 0,    'E','HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
    text(ax, -1.08, 0,    'W','HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
    for pl = 10:10:80
        r = sqrt(2)*sin((pi/2-deg2rad(pl))/2);
        plot(ax, r*cos(th), r*sin(th), '-', 'Color',[0.78 0.78 0.78], 'LineWidth',0.35);
    end
    for az = 0:10:170
        pArr = linspace(0,90,91);
        xA = arrayfun(@(p) sqrt(2)*sin((pi/2-deg2rad(p))/2)*sin(deg2rad(az)),   pArr);
        yA = arrayfun(@(p) sqrt(2)*sin((pi/2-deg2rad(p))/2)*cos(deg2rad(az)),   pArr);
        xB = arrayfun(@(p) sqrt(2)*sin((pi/2-deg2rad(p))/2)*sin(deg2rad(az+180)),pArr);
        yB = arrayfun(@(p) sqrt(2)*sin((pi/2-deg2rad(p))/2)*cos(deg2rad(az+180)),pArr);
        plot(ax, [xA xB(end:-1:1)], [yA yB(end:-1:1)], '-', ...
            'Color',[0.78 0.78 0.78], 'LineWidth',0.35);
    end
    xlim(ax,[-1.15 1.15]); ylim(ax,[-1.15 1.15]);
end

function closeFigure(ax, titleStr, dataTag, opt, C1,C1e,C2,C2e,C3,C3e)
    % 基本圖例：K1 K2 K3 符號
    h1=plot(ax,NaN,NaN,'s','MarkerFaceColor',C1,'MarkerEdgeColor',C1e,...
        'MarkerSize',13,'DisplayName','K_1 (max)');
    h2=plot(ax,NaN,NaN,'^','MarkerFaceColor',C2,'MarkerEdgeColor',C2e,...
        'MarkerSize',13,'DisplayName','K_2 (int)');
    h3=plot(ax,NaN,NaN,'o','MarkerFaceColor',C3,'MarkerEdgeColor',C3e,...
        'MarkerSize',13,'DisplayName','K_3 (min)');
    handles=[h1 h2 h3];

    if strcmp(dataTag,'raw')
        % 個別方磚小符號（用 scatter 以支援 alpha）
        hs=scatter(ax,NaN,NaN,35,'s','MarkerFaceColor',C1,'MarkerEdgeColor',C1,...
            'DisplayName','個別方磚');
        hs.MarkerFaceAlpha=0.35; hs.MarkerEdgeAlpha=0.55;
        handles(end+1)=hs;

        doB=ismember(opt.EllipseType,{'bootstrap','both'});
        doJ=ismember(opt.EllipseType,{'jelinek','both'});
        if doB
            % Bootstrap 橢圓：三色各自說明
            hb1=plot(ax,NaN,NaN,'-','Color',C1,'LineWidth',1.8,'DisplayName','Bootstrap K_1');
            hb2=plot(ax,NaN,NaN,'-','Color',C2,'LineWidth',1.8,'DisplayName','Bootstrap K_2');
            hb3=plot(ax,NaN,NaN,'-','Color',C3,'LineWidth',1.8,'DisplayName','Bootstrap K_3');
            handles(end+1:end+3)=[hb1 hb2 hb3];
        end
        if doJ
            % Jelinek 橢圓：三色各自說明
            hj1=plot(ax,NaN,NaN,'--','Color',C1,'LineWidth',1.3,'DisplayName','Jelinek K_1');
            hj2=plot(ax,NaN,NaN,'--','Color',C2,'LineWidth',1.3,'DisplayName','Jelinek K_2');
            hj3=plot(ax,NaN,NaN,'--','Color',C3,'LineWidth',1.3,'DisplayName','Jelinek K_3');
            handles(end+1:end+3)=[hj1 hj2 hj3];
        end
    end

    legend(ax,handles,'Location','southoutside','Orientation','horizontal','FontSize',9);
    title(ax,titleStr,'FontSize',13,'FontWeight','bold');
    annotation('textbox',[0.12 0.01 0.76 0.04], ...
        'String','Lower hemisphere, equal-area projection (Schmidt Net)', ...
        'EdgeColor','none','FontSize',7,'Color',[0.5 0.5 0.5],'HorizontalAlignment','center');
    hold(ax,'off');
    axis(ax,'off');

    % 儲存
    if opt.SaveFig
        if ~isempty(opt.GroupLabels)
            ftag=['Group_' strjoin(normGroup(opt.GroupLabels),'')];
        elseif ~isempty(opt.Numbers)
            ftag=['No_' strjoin(arrayfun(@num2str,opt.Numbers,'UniformOutput',false),'_')];
        else
            ftag='All';
        end
        fname=fullfile(opt.SavePath, sprintf('Stereonet_%s_%s.png',dataTag,ftag));
        exportgraphics(gcf,fname,'Resolution',300);
        fprintf('Stereonet saved → %s\n', fname);
    end
end


%% ================================================================
%  繪圖輔助
%% ================================================================
function drawLarge(ax, d1,i1, d2,i2, d3,i3, C1,C1e,C2,C2e,C3,C3e)
    [x1,y1]=schmidtProject(d1,i1);
    [x2,y2]=schmidtProject(d2,i2);
    [x3,y3]=schmidtProject(d3,i3);
    plot(ax,x1,y1,'s','MarkerFaceColor',C1,'MarkerEdgeColor',C1e,'MarkerSize',13,'LineWidth',1.2);
    plot(ax,x2,y2,'^','MarkerFaceColor',C2,'MarkerEdgeColor',C2e,'MarkerSize',13,'LineWidth',1.2);
    plot(ax,x3,y3,'o','MarkerFaceColor',C3,'MarkerEdgeColor',C3e,'MarkerSize',13,'LineWidth',1.2);
end

function drawSmall(ax, d1,i1, d2,i2, d3,i3, C1,C2,C3)
    [x1,y1]=schmidtProject(d1,i1);
    [x2,y2]=schmidtProject(d2,i2);
    [x3,y3]=schmidtProject(d3,i3);
    h1=scatter(ax,x1,y1,40,'s','MarkerFaceColor',C1,'MarkerEdgeColor',C1,'LineWidth',0.5);
    h2=scatter(ax,x2,y2,40,'^','MarkerFaceColor',C2,'MarkerEdgeColor',C2,'LineWidth',0.5);
    h3=scatter(ax,x3,y3,40,'o','MarkerFaceColor',C3,'MarkerEdgeColor',C3,'LineWidth',0.5);
    h1.MarkerFaceAlpha=0.30; h1.MarkerEdgeAlpha=0.50;
    h2.MarkerFaceAlpha=0.30; h2.MarkerEdgeAlpha=0.50;
    h3.MarkerFaceAlpha=0.30; h3.MarkerEdgeAlpha=0.50;
end

function addLabel(ax, trend, plunge, lbl, color, fs)
    [x,y]=schmidtProject(trend,plunge);
    text(ax,x+0.03,y+0.03,lbl,'FontSize',fs,'Color',color);
end

function lbl = meanLabel(results, g)
    if isfield(results,'GroupLabel') && ~isempty(results.GroupLabel{g})
        lbl=sprintf('%s-No.%d', results.GroupLabel{g}, results.Number(g));
    else
        lbl=sprintf('No.%d', results.Number(g));
    end
end


%% ================================================================
%  Bootstrap 橢圓
%% ================================================================

function drawBootEllipse(ax, kDir, otherDir1, otherDir2, bc, which, color)
% 以 kDir 為中心畫 Bootstrap PCA 橢圓
% otherDir1/otherDir2 用於確認橢圓軸方向（目前由 bc 內部記錄，不需額外使用）
    halfMaj = bc.([which '_halfAngle']);
    halfMin = bc.([which '_minorAngle']);
    perpT   = bc.([which '_perpTrend']);
    perpP   = bc.([which '_perpPlunge']);
    if isnan(halfMaj) || halfMaj <= 0, return; end
    drawEllipse(ax, kDir(1),kDir(2), halfMaj,halfMin, perpT,perpP, color, '-', 1.8);
end

function bc = calcBootstrap(tensorStack, nBoot)
% 計算三個主軸各自的 Bootstrap PCA 橢圓參數
    N = size(tensorStack,3);
    bc = struct( ...
        'K1_halfAngle',NaN,'K1_minorAngle',NaN,'K1_perpTrend',NaN,'K1_perpPlunge',NaN, ...
        'K2_halfAngle',NaN,'K2_minorAngle',NaN,'K2_perpTrend',NaN,'K2_perpPlunge',NaN, ...
        'K3_halfAngle',NaN,'K3_minorAngle',NaN,'K3_perpTrend',NaN,'K3_perpPlunge',NaN);
    if N < 3, return; end

    M0 = AMSFormulas.calculateMeanTensor(tensorStack);
    e0 = AMSFormulas.performEigenAnalysis(M0);
    v0 = { tp2v(e0.K1_dir(1),e0.K1_dir(2)), ...
           tp2v(e0.K2_dir(1),e0.K2_dir(2)), ...
           tp2v(e0.K3_dir(1),e0.K3_dir(2)) };

    bV = zeros(3, nBoot, 3);   % [xyz, bootstrap樣本, 軸編號]
    for b = 1:nBoot
        si = randi(N, N, 1);
        Mb = AMSFormulas.calculateMeanTensor(tensorStack(:,:,si));
        eb = AMSFormulas.performEigenAnalysis(Mb);
        vb = { tp2v(eb.K1_dir(1),eb.K1_dir(2)), ...
               tp2v(eb.K2_dir(1),eb.K2_dir(2)), ...
               tp2v(eb.K3_dir(1),eb.K3_dir(2)) };
        for k = 1:3
            vbk = vb{k};
            if dot(vbk, v0{k}) < 0, vbk = -vbk; end
            bV(:,b,k) = vbk;
        end
    end

    keys = {'K1','K2','K3'};
    for k = 1:3
        [ha,ma,pT,pP] = fitBootEllipse(squeeze(bV(:,:,k)), v0{k});
        bc.([keys{k} '_halfAngle'])  = ha;
        bc.([keys{k} '_minorAngle']) = ma;
        bc.([keys{k} '_perpTrend'])  = pT;
        bc.([keys{k} '_perpPlunge']) = pP;
    end
end

function [halfAngle, minorAngle, perpTrend, perpPlunge] = fitBootEllipse(bVecs, v0)
% Bootstrap 散布橢圓
%
% 長短軸方向：用 Bingham 方向散布矩陣（T = bVecs*bVecs'/N）特徵值分解
%             完全在球面上定義，不依賴切平面近似
% 半角大小：  投影到各主軸方向後取球面角距的 95% 分位數
%
% 修改說明：
%   舊版用切平面 PCA 找長短軸方向，靠近赤道時方向有近似誤差。
%   新版改用 3×3 方向散布矩陣的特徵向量，與 Jelinek 橢圓的幾何定義一致。

    halfAngle = NaN; minorAngle = NaN; perpTrend = NaN; perpPlunge = NaN;
    nB = size(bVecs, 2);
    if nB < 10, return; end

    % --- Step 1：Bingham 方向散布矩陣特徵值分解 -------------------
    % T 的特徵向量即橢圓主軸方向（球面上精確定義）
    % 對應 v0 的特徵向量是最大特徵值（最集中方向）
    % 另外兩個特徵向量是橢圓的長軸（次大）和短軸（最小）
    T = (bVecs * bVecs') / nB;   % 3×3 方向散布矩陣
    [eV, eD] = eig(T);
    [~, si] = sort(diag(eD), 'descend');

    % si(1) 對應 v0（最集中方向），si(2) 長軸，si(3) 短軸
    % 確認 si(1) 確實對應 v0（與 v0 的夾角最小）
    dots = abs([dot(eV(:,si(1)),v0), dot(eV(:,si(2)),v0), dot(eV(:,si(3)),v0)]);
    [~, v0idx] = max(dots);
    axOrder = si([v0idx==1, v0idx==2, v0idx==3]);   % 重排：[v0軸, 長軸, 短軸]
    % 如果 v0idx==1，axOrder = si([1,2,3]) = si；其他情況重排
    allIdx = [1 2 3];
    otherIdx = allIdx(allIdx ~= v0idx);
    majVec = eV(:, si(otherIdx(1)));   % 長軸特徵向量（球面上）
    minVec = eV(:, si(otherIdx(2)));   % 短軸特徵向量（球面上）

    % 確保方向朝向下半球（與 v0 同側）
    if dot(majVec, v0) < 0, majVec = -majVec; end
    if dot(minVec, v0) < 0, minVec = -minVec; end

    % --- Step 2：計算各 bootstrap 樣本沿長/短軸方向的球面角距 -----
    aMaj = zeros(nB, 1);
    aMin = zeros(nB, 1);

    for b = 1:nB
        vb  = bVecs(:,b) / norm(bVecs(:,b));

        % 與平均方向 v0 的總球面角距
        ang = acosd(min(1, max(-1, dot(vb, v0))));

        % 切平面投影（只用來分解方向，不用來找軸）
        vp = vb - dot(vb, v0) * v0;
        if norm(vp) > 1e-10
            vp = vp / norm(vp);
            % 沿長軸方向的角距分量
            aMaj(b) = ang * dot(vp, majVec);
            % 沿短軸方向的角距分量
            aMin(b) = ang * dot(vp, minVec);
        end
    end

    % --- Step 3：95% 分位數當作半角 --------------------------------
    halfAngle  = prctile(abs(aMaj), 95);
    minorAngle = prctile(abs(aMin), 95);

    % --- Step 4：長軸方向轉為 Trend/Plunge（供 drawEllipse 使用）--
    pv = majVec;
    if pv(3) < 0, pv = -pv; end
    perpPlunge = asind(min(1, max(-1, pv(3))));
    perpTrend  = atan2d(pv(2), pv(1));
    if perpTrend < 0, perpTrend = perpTrend + 360; end
end


%% ================================================================
%  Jelinek 橢圓（真正的橢圓，非圓形）
%% ================================================================

function drawJelEllipse(ax, kDir, adjDir1, adjDir2, halfMaj_deg, halfMin_deg, color)
% 以 kDir 為中心畫 Jelinek 橢圓
%
% 橢圓定義：
%   長軸方向 = kDir 到 adjDir1 所在的大圓平面內、垂直於 kDir 的方向
%   短軸方向 = kDir 到 adjDir2 所在的大圓平面內、垂直於 kDir 的方向
%   長軸半角 = halfMaj_deg（對應 adjDir1 那對的 Jelinek 角，如 E12）
%   短軸半角 = halfMin_deg（對應 adjDir2 那對的 Jelinek 角，如 E13）
%
% 輸入：
%   kDir       - 本軸方向 [trend, plunge]（度）
%   adjDir1    - 長軸對應的相鄰軸方向 [trend, plunge]（度）
%   adjDir2    - 短軸對應的相鄰軸方向 [trend, plunge]（度）
%   halfMaj_deg - 長軸半角（Jelinek Eij，度）
%   halfMin_deg - 短軸半角（Jelinek Eik，度）
%   color      - 線條顏色

    if isnan(halfMaj_deg) || halfMaj_deg <= 0, return; end
    if isnan(halfMin_deg) || halfMin_deg <= 0, halfMin_deg = halfMaj_deg; end

    v0  = tp2v(kDir(1),    kDir(2));
    va1 = tp2v(adjDir1(1), adjDir1(2));
    va2 = tp2v(adjDir2(1), adjDir2(2));

    % 長軸方向：在 v0–va1 大圓平面內，垂直於 v0
    perpMaj = va1 - dot(va1, v0) * v0;
    if norm(perpMaj) < 1e-6
        % va1 與 v0 幾乎平行，退化為圓形
        halfMin_deg = halfMaj_deg;
        if abs(v0(3)) < 0.9, ref = [0;0;1]; else, ref = [1;0;0]; end
        perpMaj = cross(v0, ref);
    end
    perpMaj = perpMaj / norm(perpMaj);
    if perpMaj(3) < 0, perpMaj = -perpMaj; end

    % 短軸方向：在 v0–va2 大圓平面內，垂直於 v0
    perpMin = va2 - dot(va2, v0) * v0;
    if norm(perpMin) < 1e-6
        perpMin = cross(v0, perpMaj);
    end
    perpMin = perpMin / norm(perpMin);

    % 取得 perpMaj 的 trend/plunge 供 drawEllipse 使用
    if perpMaj(3) < 0, perpMaj = -perpMaj; end
    pP = asind(perpMaj(3));
    pT = atan2d(perpMaj(2), perpMaj(1));
    if pT < 0, pT = pT + 360; end

    % 確認長短軸方向正交；若不正交則以 perpMaj×v0 為短軸
    if abs(dot(perpMaj, perpMin)) > 0.1
        perpMin = cross(v0, perpMaj);
        perpMin = perpMin / norm(perpMin);
    end

    drawEllipse(ax, kDir(1),kDir(2), halfMaj_deg,halfMin_deg, pT,pP, color, '--', 1.3);
end


%% ================================================================
%  球面橢圓核心（Rodrigues 旋轉 + 正確跨半球處理）
%% ================================================================
function drawEllipse(ax, t0,p0, halfMaj,halfMin, perpT,perpP, color, ls, lw)
% 在球面上繪製橢圓（Schmidt 等面積投影，下半球）
%
% 跨半球處理說明：
%   橢圓上的點若跑到上半球（vo(3) < 0），取對徑點（-vo）投影到下半球。
%   對徑點在投影圖上對稱出現在另一側，需插入 NaN 斷開連線，
%   兩段弧分開繪製，視覺上呈現為兩條對稱弧線（標準下半球慣例）。

    if isnan(halfMaj) || halfMaj <= 0, return; end
    if isnan(halfMin) || halfMin <= 0, halfMin = halfMaj; end
    halfMaj = min(halfMaj, 89.9);
    halfMin = min(halfMin, 89.9);

    v0 = tp2v(t0, p0);

    % 建立長軸方向 vP（垂直於 v0）
    vP = tp2v(perpT, perpP);
    vP = vP - dot(vP, v0) * v0;
    if norm(vP) < 1e-6
        if abs(v0(3)) < 0.9, ref = [0;0;1]; else, ref = [1;0;0]; end
        vP = cross(v0, ref);
    end
    vP  = vP  / norm(vP);
    % 短軸方向：與 v0、vP 均垂直
    vP2 = cross(v0, vP);
    vP2 = vP2 / norm(vP2);

    % 生成橢圓上的 360 個點
    nPts = 361;
    phi  = linspace(0, 2*pi, nPts);
    vArr = zeros(3, nPts);   % 原始向量（可能在上半球）
    hArr = zeros(1, nPts);   % 1=下半球, -1=上半球（對徑反射）

    for k = 1:nPts
        vo = rotVec(v0,  vP,  deg2rad(halfMaj * cos(phi(k))));
        vo = rotVec(vo,  vP2, deg2rad(halfMin * sin(phi(k))));
        vo = vo / norm(vo);
        if vo(3) >= 0
            vArr(:,k) = vo;
            hArr(k)   = 1;
        else
            vArr(:,k) = -vo;   % 對徑反射
            hArr(k)   = -1;
        end
    end

    % 找半球切換點，在切換前後插入 NaN 斷開
    xe = nan(1, nPts);
    ye = nan(1, nPts);
    for k = 1:nPts
        % 若與前一點不同半球，保持 NaN（斷點已由初始化設定）
        if k > 1 && hArr(k) ~= hArr(k-1)
            % 斷點：此點開始新段，保留值但前一點後補 NaN 已由下一段起始處理
            % 實際做法：直接跳過，讓 NaN 保持
            continue;
        end
        vo = vArr(:,k);
        pk = asind(min(1, max(-1, vo(3))));
        tk = atan2d(vo(2), vo(1));
        if tk < 0, tk = tk + 360; end
        [xe(k), ye(k)] = schmidtProject(tk, pk);
    end

    % 分段繪製（下半球段 / 對徑段各自連續）
    % 重新整理：收集連續段
    segs = splitSegments(xe, ye, hArr);
    for s = 1:numel(segs)
        plot(ax, segs{s}(1,:), segs{s}(2,:), ls, 'Color', color, 'LineWidth', lw);
    end
end

function segs = splitSegments(xe, ye, hArr)
% 依半球標記把橢圓點序列切成連續段
% 若橢圓完全在同一半球，首尾相接形成閉合曲線
    nPts = numel(xe);
    segs = {};
    xSeg = []; ySeg = []; curH = hArr(1);

    for k = 1:nPts
        if hArr(k) == curH && ~isnan(xe(k))
            xSeg(end+1) = xe(k); %#ok<AGROW>
            ySeg(end+1) = ye(k); %#ok<AGROW>
        else
            if numel(xSeg) >= 2
                segs{end+1} = [xSeg; ySeg]; %#ok<AGROW>
            end
            xSeg = []; ySeg = []; curH = hArr(k);
            if ~isnan(xe(k))
                xSeg = xe(k); ySeg = ye(k);
            end
        end
    end
    % 收尾最後一段
    if numel(xSeg) >= 2
        segs{end+1} = [xSeg; ySeg]; %#ok<AGROW>
    end

    % 若只有一段且首尾半球相同（完整橢圓），閉合它
    if numel(segs) == 1
        segs{1}(:,end+1) = segs{1}(:,1);
    end

    % 若有兩段且首尾半球相同，把最後一段接到第一段前面（橢圓跨邊界後接合）
    if numel(segs) == 2 && hArr(1) == hArr(end)
        segs{1} = [segs{2}, segs{1}, segs{2}(:,1)];
        segs(2) = [];
    end
end


%% ================================================================
%  幾何工具
%% ================================================================
function [x,y] = schmidtProject(trend_deg, plunge_deg)
% Schmidt 等面積投影（下半球）
    r = sqrt(2) * sin((pi/2 - deg2rad(plunge_deg)) / 2);
    x = r * sin(deg2rad(trend_deg));
    y = r * cos(deg2rad(trend_deg));
end

function v = schmidtInv(x, y)
% Schmidt 投影逆運算（還原單位向量）
    r2 = x^2 + y^2;
    z  = max(0, 1 - r2/2);
    v  = [x * sqrt(1 - r2/4); y * sqrt(1 - r2/4); z];
    if norm(v) > 1e-10, v = v / norm(v); end
end

function v = tp2v(trend_deg, plunge_deg)
% Trend/Plunge → 單位方向向量（地理座標 X=北, Y=東, Z=向下）
    t  = deg2rad(trend_deg);
    pp = deg2rad(plunge_deg);
    v  = [cos(pp)*cos(t); cos(pp)*sin(t); sin(pp)];
end

function vr = rotVec(v, axis, angle_rad)
% Rodrigues 旋轉公式：v 繞 axis 旋轉 angle_rad
    vr = v * cos(angle_rad) + cross(axis, v) * sin(angle_rad) + ...
         axis * dot(axis, v) * (1 - cos(angle_rad));
end

function out = normGroup(cellIn)
    out = cellfun(@normStr, cellIn, 'UniformOutput', false);
end

function s = normStr(s)
    s = strtrim(strrep(s, '''', ''));
    for c = 1:numel(s)
        code = double(s(c));
        if code >= 65313 && code <= 65338
            s(c) = char(code - 65248);
        end
    end
end

function ts = makeTitle(dtype, opt)
    if ~isempty(opt.GroupLabels)
        tag = sprintf('Group %s', strjoin(opt.GroupLabels, '/'));
    elseif ~isempty(opt.Numbers)
        tag = sprintf('No.%s', strjoin(arrayfun(@num2str, opt.Numbers, 'UniformOutput', false), ','));
    else
        tag = 'All';
    end
    if strcmp(dtype,'mean'), src = '平均主軸'; else, src = '原始方磚'; end
    ts = sprintf('AMS Stereonet [%s] — %s', src, tag);
end