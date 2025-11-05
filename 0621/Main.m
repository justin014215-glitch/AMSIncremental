% ==========================
% 主程式：AMS 應變分析器
% ==========================

%config.CleanSmallValues = true;   % 關掉就設為 false
%config.Threshold = 1e-12;         % 調整你覺得多小要忽略

function main()
    % 主程式：AMS 磁感率異向性應變增量分析器
    clc; close all;

    % 檔案選擇
    [filename, pathname] = uigetfile({'*.xlsx;*.csv;*.txt'}, '請選擇 AMS 數據檔');
    if isequal(filename, 0)
        disp('使用者取消選擇。');
        return;
    end

    % 建立分析物件
    analyzer = AMSStrainAnalyzer(fullfile(pathname, filename), ...
        'SlateCoeffA', 6.897, ...
        'SlateCoeffB', 0.007, ...
        'Verbose', true);

    fprintf('\n=== 北橫磁感率異向性應變增量分析啟動 ===\n');
    
    % 計算原始磁化率張量與 Eg_raw
    analyzer.computeErRawAll();
    analyzer.computeEgFromErRaw();

    % 計算 3D Eg 並儲存變數
    analyzer.computeFiniteStrainTensors();
    assignSampleResultsToWorkspace(analyzer, '3D');
    % 設定是否清理小數值（如 1e-15）開關
    config.CleanSmallValues = true;
    config.Threshold = 1e-10;

    % 儲存分析物件供後續使用
    assignin('base', 'ams_analyzer', analyzer);

    fprintf('\n✅ 所有樣本應變結果已儲存至工作區\n');
    fprintf('🔍 範例：輸入 `Er_1`、`V_3` 或 `Eg_5` 存取對應結果\n');

    % 啟動互動選單（選擇性）
    ask = input('\n是否要進入互動選單模式？(y/n)：', 's');
    if strcmpi(ask, 'y')
        interactiveMenu(analyzer,config);
    else
        disp('✅ 分析結束。結果已可直接使用。');
    end


end
function A_clean = cleanMatrix(A, threshold)
    % 將絕對值小於 threshold 的元素設為 0
    if nargin < 2
        threshold = 1e-10;
    end
    A_clean = A;
    A_clean(abs(A) < threshold) = 0;
end

function assignSampleResultsToWorkspace(analyzer, mode)
    switch mode
        case '3D'
            for i = 1:size(analyzer.Eg, 3)
                assignin('base', sprintf('Er_%d', i), analyzer.ErList(:,:,i));
                assignin('base', sprintf('V_%d', i),  analyzer.VList(:,:,i));
                assignin('base', sprintf('Eg_%d', i), analyzer.Eg(:,:,i));
            end
        case '2D'
            for i = 1:size(analyzer.EgD, 3)
                assignin('base', sprintf('Er_2D_%d', i), analyzer.ErList2D(:,:,i));
                assignin('base', sprintf('V_2D_%d', i),  analyzer.VList2D(:,:,i));
                assignin('base', sprintf('Eg_2D_%d', i), analyzer.EgD(:,:,i));
            end
    end
end

%%
function combineAndAnalyzeEgVars(varA, varB, config)
    % 預設參數
    if nargin < 3
        config.CleanSmallValues = true;
        config.Threshold = 1e-10;
    end

    % 嘗試取得變數
    try
        Eg_A = evalin('base', varA);
        Eg_B = evalin('base', varB);
        
    catch
        fprintf('❌ 找不到其中一個變數：%s 或 %s\n', varA, varB);
        return;
    end
    V_A = tryGetVFromEgName(varA);
    V_B = tryGetVFromEgName(varB);
    % 執行 Eg 相乘
    Eg_combined = Eg_B * Eg_A;
    Eg_sym = (Eg_combined + Eg_combined') / 2;

    % 可選：清除極小值
    if config.CleanSmallValues
        Eg_combined = cleanMatrix(Eg_combined, config.Threshold);
        Eg_sym = cleanMatrix(Eg_sym, config.Threshold);
    end

    % 特徵值分析（並排序）
    [V, D] = eig(Eg_sym);
    [d, idx] = sort(diag(D), 'descend');
    D_sorted = diag(d);
    V_sorted = V(:, idx);

    % 找下一個命名編號
    existingVars = evalin('base', 'who');
    matched = regexp(existingVars, '^Eg_combined_\d+$', 'match');
    nextID = sum(~cellfun('isempty', matched)) + 1;
    varName = @(base) sprintf('%s_%d', base, nextID);

    % 顯示處理流程
    fprintf('\n🔗 Eg 組合分析（第 %d 次）: Eg_B * Eg_A\n', nextID);
    fprintf('Eg_A (%s):\n',varA), disp(Eg_A);
    fprintf('V_A(%s):\n',varA), disp(V_A);
    fprintf('V_A^T(%s):\n',varA),disp(V_A')
    fprintf('Eg_B (%s):\n',varB), disp(Eg_B);
    fprintf('V_B(%s):\n',varB), disp(V_B);
    fprintf('V_B^T(%s):\n',varB),disp(V_B')
    fprintf('Eg_combined = Eg_B * Eg_A:\n'), disp(Eg_combined);
    fprintf('對稱化 Eg_sym:\n'), disp(Eg_sym);
    fprintf('特徵值 (排序):\n'), disp(D_sorted);
    fprintf('特徵向量 (排序):\n'), disp(V_sorted);

    % 儲存結果到 base workspace
    assignin('base', varName('Eg_combined'), Eg_combined);
    assignin('base', varName('Eg_sym'), Eg_sym);
    assignin('base', varName('Eg_combined_V'), V_sorted);
    assignin('base', varName('Eg_combined_D'), D_sorted);

    % 儲存來源追蹤（可選）
    assignin('base', varName('Eg_combined_sources'), {varA, varB});

    % 完成通知
    fprintf('\n✅ 儲存變數：%s, %s, %s, %s\n來源記錄：%s\n', ...
        varName('Eg_combined'), varName('Eg_sym'), ...
        varName('Eg_combined_V'), varName('Eg_combined_D'), ...
        varName('Eg_combined_sources'));
end

%%

function combineAndAnalyzeEgAverage(varList, config)
    % 預設參數處理
    if nargin < 2
        config.CleanSmallValues = true;
        config.Threshold = 1e-10;
    end

    % 自動偵測變數
    if nargin < 1 || isempty(varList)
        allVars = evalin('base', 'who');
        matched = regexp(allVars, '^(Eg_|Eg_raw_)\d+$', 'match');
        varList = [matched{:}];
        if isempty(varList)
            disp('❌ 找不到任何 Eg_* 或 Eg_raw_* 變數可用於平均');
            return;
        end
        fprintf('📦 自動偵測到以下 Eg 變數將進行平均：\n');
        disp(varList');
    elseif ~iscell(varList)
        disp('❌ 請以 cell array 格式輸入，如 {''Eg_1'', ''Eg_raw_2''}');
        return;
    end

    % 讀取 Eg 矩陣
    Eg_matrices = [];
    for i = 1:length(varList)
        try
            Eg_i = evalin('base', varList{i});
            Eg_matrices(:, :, i) = Eg_i;
        catch
            fprintf('❌ 找不到變數：%s\n', varList{i});
            return;
        end
    end

    % 計算平均與對稱化
    Eg_avg = mean(Eg_matrices, 3);
    Eg_sym = (Eg_avg + Eg_avg') / 2;

    if config.CleanSmallValues
        Eg_avg = cleanMatrix(Eg_avg, config.Threshold);
        Eg_sym = cleanMatrix(Eg_sym, config.Threshold);
    end

    % 特徵分析
    [V, D] = eig(Eg_sym);
    [d, idx] = sort(diag(D), 'descend');
    D_sorted = diag(d);
    V_sorted = V(:, idx);

    % 命名與編號
    existingVars = evalin('base', 'who');
    matched = regexp(existingVars, '^Eg_avg_\d+$', 'match');
    nextID = sum(~cellfun('isempty', matched)) + 1;
    varName = @(base) sprintf('%s_%d', base, nextID);

    % 顯示處理過程
    fprintf('\n📊 平均 Eg (%d 個):\n', length(varList)); disp(Eg_avg);
    fprintf('\n🔄 對稱化:\n'); disp(Eg_sym);
    fprintf('\n📐 特徵值:\n'); disp(D_sorted);
    fprintf('\n🧭 特徵向量:\n'); disp(V_sorted);

    % 儲存至 base workspace
    assignin('base', varName('Eg_avg'), Eg_avg);
    assignin('base', varName('Eg_avg_sym'), Eg_sym);
    assignin('base', varName('Eg_avg_V'), V_sorted);
    assignin('base', varName('Eg_avg_D'), D_sorted);
    assignin('base', varName('Eg_avg_sources'), varList);

    % 回報儲存
    fprintf('\n✅ 平均後結果儲存為：\n%s\n%s\n%s\n%s\n來源清單：%s\n', ...
        varName('Eg_avg'), varName('Eg_avg_sym'), ...
        varName('Eg_avg_V'), varName('Eg_avg_D'), varName('Eg_avg_sources'));
end


function V = tryGetVFromEgName(egName)
    % 嘗試從 egName 推出對應的特徵向量名稱，支援 Eg_4 與 Eg_raw_4
    V = NaN;
    try
        % Eg_4 → V_4；Eg_raw_4 → V_4
        id = regexp(egName, '\d+$', 'match');
        if ~isempty(id)
            Vname = ['V_' id{1}];  % 產生 V_4
            V = evalin('base', Vname);
        end
    catch
        % 無對應值 → 傳回 NaN
    end
end
%%
function interactiveMenu(analyzer,config)
    while true
        disp('--- 分析功能選單 ---');
        disp('1. 顯示前 3 筆 3D 結果');
        %disp('2. 顯示前 3 筆 2D 結果');
        disp('3. 顯示任一筆樣本 3D 結果');
        %disp('4. 顯示任一筆樣本 2D 結果');
        disp('5. 合併兩筆 Eg 並分析（Eg_new = Eg_B * Eg_A）');
        disp('6. 對多個 Eg / Eg_raw 做平均後分析');
        disp('0. 離開');
        choice = input('請輸入功能編號：');
    
        switch choice
            case 1
                for i = 1:min(3, size(analyzer.Eg,3))
                    fprintf('\n樣本 %d（3D）:\n', i);
                    disp(['Er_' num2str(i)]); disp(analyzer.ErList(:,:,i));
                    disp(['V_' num2str(i)]);  disp(analyzer.VList(:,:,i));
                    disp(['Eg_' num2str(i)]); disp(analyzer.Eg(:,:,i));
                end
      %{     
            case 2
                for i = 1:min(3, size(analyzer.EgD,3))
                    fprintf('\n樣本 %d（2D）:\n', i);
                    disp(['Er_2D_' num2str(i)]); disp(analyzer.ErList2D(:,:,i));
                    disp(['V_2D_' num2str(i)]);  disp(analyzer.VList2D(:,:,i));
                    disp(['Eg_2D_' num2str(i)]); disp(analyzer.EgD(:,:,i));
                end
      %}
            case 3
                idx = input('請輸入樣本編號：');
                fprintf('\n樣本 %d（3D）:\n', idx);
                disp(['Er_' num2str(idx)]); disp(analyzer.ErList(:,:,idx));
                disp(['V_' num2str(idx)]);  disp(analyzer.VList(:,:,idx));
                disp(['Eg_' num2str(idx)]); disp(analyzer.Eg(:,:,idx));
      %{
            case 4
                idx = input('請輸入樣本編號：');
                fprintf('\n樣本 %d（2D）:\n', idx);
                disp(['Er_2D_' num2str(idx)]); disp(analyzer.ErList2D(:,:,idx));
                disp(['V_2D_' num2str(idx)]);  disp(analyzer.VList2D(:,:,idx));
                disp(['Eg_2D_' num2str(idx)]); disp(analyzer.EgD(:,:,idx));
      %}
            case 5
                varA = input('請輸入第一個 Eg 變數名稱：', 's');
                varB = input('請輸入第二個 Eg 變數名稱：', 's');
                combineAndAnalyzeEgVars(varA, varB, config);

            case 6
                varList = input('請輸入 Eg 變數（如 {''Eg_1'', ''Eg_raw_3''}，留空則自動）：');
                combineAndAnalyzeEgAverage(varList, config);

            case 0
                disp('✅ 離開選單，分析結束！');
                break;
            otherwise
                disp('❌ 請輸入有效選項');
        end
        if nargin < 2
        config.CleanSmallValues = true;
        config.Threshold = 1e-10;
        end
    end
end

