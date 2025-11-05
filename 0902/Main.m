% ==========================
% 主程式：AMS 應變分析器
% 北橫磁感率異向性應變增量分析器
% ==========================

function main()
    % 主程式：AMS 磁感率異向性應變增量分析器
    clc; close all;

    fprintf('\n╔══════════════════════════════════════╗\n');
    fprintf('║     北橫磁感率異向性應變增量分析     ║\n');
    fprintf('║        AMS Strain Analyzer v2.0      ║\n');
    fprintf('╚══════════════════════════════════════╝\n\n');

    % 檔案選擇
    [filename, pathname] = uigetfile({'*.xlsx;*.csv;*.txt'}, '請選擇 AMS 數據檔');
    if isequal(filename, 0)
        disp('使用者取消選擇。');
        return;
    end

    % 建立分析物件
    try
        analyzer = AMSStrainAnalyzer(fullfile(pathname, filename), ...
            'SlateCoeffA', 6.897, ...
            'SlateCoeffB', 0.007, ...
            'Verbose', true, ...
            'CleanThreshold', 1e-12, ...
            'ValidateResults', true);
    catch ME
        error('建立分析器失敗：%s', ME.message);
    end

    fprintf('\n=== 分析流程開始 ===\n');
    
    try
        % 計算原始磁化率張量與 Eg_raw（驗證用）
        fprintf('1. 計算原始磁化率張量...\n');
        analyzer.computeErRawAll();
        analyzer.computeEgFromErRaw();

        % 計算 3D 有限應變張量
        fprintf('2. 計算有限應變張量...\n');
        analyzer.computeFiniteStrainTensors();
        
        % 儲存結果到工作區
        assignSampleResultsToWorkspace(analyzer, '3D');
        
        % 列印分析摘要
        analyzer.printSummary();
        
        % 儲存分析物件供後續使用
        assignin('base', 'ams_analyzer', analyzer);

        fprintf('\n✅ 所有樣本應變結果已儲存至工作區\n');
        fprintf('🔍 範例：輸入 `Er_1`、`V_3` 或 `Eg_5` 存取對應結果\n');
        fprintf('📊 使用 `ams_analyzer.printSummary()` 查看分析摘要\n');

    catch ME
        error('分析過程發生錯誤：%s', ME.message);
    end

    % 啟動互動選單（選擇性）
    ask = input('\n是否要進入互動選單模式？(y/n)：', 's');
    if strcmpi(ask, 'y')
        interactiveMenu(analyzer);
    else
        disp('✅ 分析結束。結果已可直接使用。');
        
        % 詢問是否匯出結果
        exportAsk = input('是否要匯出分析結果？(y/n)：', 's');
        if strcmpi(exportAsk, 'y')
            exportResults(analyzer);
        end
    end
end

function assignSampleResultsToWorkspace(analyzer, mode)
    % 將樣本結果指派到工作區
    switch mode
        case '3D'
            n = size(analyzer.Eg, 3);
            for i = 1:n
                assignin('base', sprintf('Er_%d', i), analyzer.ErList(:,:,i));
                assignin('base', sprintf('V_%d', i),  analyzer.VList(:,:,i));
                assignin('base', sprintf('Eg_%d', i), analyzer.Eg(:,:,i));
            end
            if analyzer.config.Verbose
                fprintf('已建立 %d 組變數：Er_1~%d, V_1~%d, Eg_1~%d\n', n, n, n, n);
            end
    end
end

function exportResults(analyzer)
    % 匯出結果的子程式
    formats = {'Excel (.xlsx)', 'CSV (.csv)', 'MATLAB (.mat)'};
    [idx, tf] = listdlg('PromptString', '選擇匯出格式：', ...
                        'SelectionMode', 'single', ...
                        'ListString', formats);
    
    if ~tf
        disp('取消匯出。');
        return;
    end
    
    % 取得檔案名稱
    [~, baseName, ~] = fileparts(analyzer.filename);
    defaultName = sprintf('%s_應變分析結果', baseName);
    
    switch idx
        case 1  % Excel
            [filename, pathname] = uiputfile('*.xlsx', '儲存 Excel 檔案', [defaultName '.xlsx']);
            if ~isequal(filename, 0)
                analyzer.exportToExcel(fullfile(pathname, filename));
            end
        case 2  % CSV
            [filename, pathname] = uiputfile('*.csv', '儲存 CSV 檔案', [defaultName '.csv']);
            if ~isequal(filename, 0)
                analyzer.exportToCSV(fullfile(pathname, filename));
            end
        case 3  % MAT
            [filename, pathname] = uiputfile('*.mat', '儲存 MATLAB 檔案', [defaultName '.mat']);
            if ~isequal(filename, 0)
                analyzer.exportToMat(fullfile(pathname, filename));
            end
    end
end

function interactiveMenu(analyzer)
    % 互動選單介面
    while true
        fprintf('\n╔════════════════════════════════════════╗\n');
        fprintf('║              分析功能選單              ║\n');
        fprintf('╠════════════════════════════════════════╣\n');
        fprintf('║ 1. 顯示前 3 筆 3D 結果                 ║\n');
        fprintf('║ 2. 顯示指定樣本 3D 結果                ║\n');
        fprintf('║ 3. 合併兩筆 Eg 並分析                  ║\n');
        fprintf('║ 4. 對多個 Eg 做平均後分析              ║\n');
        fprintf('║ 5. 計算增量應變（兩樣本間）            ║\n');
        fprintf('║ 6. 繪製應變演化圖                      ║\n');
        fprintf('║ 7. 匯出分析結果                        ║\n');
        fprintf('║ 8. 顯示分析摘要                        ║\n');
        fprintf('║ 0. 離開                                ║\n');
        fprintf('╚════════════════════════════════════════╝\n');
        
        choice = input('請輸入功能編號：');
    
        try
            switch choice
                case 1
                    % 顯示前 3 筆結果
                    displaySampleResults(analyzer, 1:min(3, size(analyzer.Eg,3)));
                    
                case 2
                    % 顯示指定樣本結果
                    maxSample = size(analyzer.Eg, 3);
                    idx = input(sprintf('請輸入樣本編號 (1-%d)：', maxSample));
                    if idx >= 1 && idx <= maxSample
                        displaySampleResults(analyzer, idx);
                    else
                        fprintf('無效的樣本編號\n');
                    end
                    
                case 3
                    % 合併兩筆 Eg
                    varA = input('請輸入第一個 Eg 變數名稱：', 's');
                    varB = input('請輸入第二個 Eg 變數名稱：', 's');
                    config = struct('CleanSmallValues', true, 'Threshold', analyzer.config.CleanThreshold);
                    analyzer.combineAndAnalyzeEg(varA, varB, config);

                case 4
                    % 平均多個 Eg
                    fprintf('範例輸入：{''Eg_1'', ''Eg_raw_3'', ''Eg_5''}\n');
                    varList = input('請輸入 Eg 變數（留空則自動偵測）：');
                    config = struct('CleanSmallValues', true, 'Threshold', analyzer.config.CleanThreshold);
                    analyzer.averageAndAnalyzeEg(varList, config);

                case 5
                    % 增量應變分析
                    maxSample = size(analyzer.Eg, 3);
                    sampleA = input(sprintf('請輸入初始樣本編號 (1-%d)：', maxSample));
                    sampleB = input(sprintf('請輸入最終樣本編號 (1-%d)：', maxSample));
                    resultName = input('請輸入結果名稱（如 increment_1_to_5）：', 's');
                    
                    if sampleA >= 1 && sampleA <= maxSample && sampleB >= 1 && sampleB <= maxSample
                        analyzer.computeIncrementalStrain(sampleA, sampleB, resultName);
                    else
                        fprintf('無效的樣本編號\n');
                    end

                case 6
                    % 繪製應變演化圖
                    maxSample = size(analyzer.Eg, 3);
                    rangeStr = input(sprintf('請輸入樣本範圍（如 1:10，預設 1:%d）：', min(10, maxSample)), 's');
                    if isempty(rangeStr)
                        analyzer.plotStrainEvolution();
                    else
                        try
                            sampleRange = eval(rangeStr);
                            analyzer.plotStrainEvolution(sampleRange);
                        catch
                            fprintf('無效的範圍格式\n');
                        end
                    end

                case 7
                    % 匯出結果
                    exportResults(analyzer);

                case 8
                    % 顯示摘要
                    analyzer.printSummary();

                case 0
                    fprintf('✅ 離開選單，分析結束！\n');
                    break;
                    
                otherwise
                    fprintf('❌ 請輸入有效選項 (0-8)\n');
            end
        catch ME
            fprintf('❌ 操作失敗：%s\n', ME.message);
        end
        
        % 暫停讓使用者查看結果
        if choice ~= 0
            input('按 Enter 繼續...');
        end
    end
end

function displaySampleResults(analyzer, indices)
    % 顯示指定樣本的詳細結果
    if isscalar(indices)
        indices = indices:indices;
    end
    
    for i = indices
        if i > size(analyzer.Eg, 3)
            fprintf('樣本 %d 超出範圍\n', i);
            continue;
        end
        
        fprintf('\n━━━━━━━━ 樣本 %d 詳細結果 ━━━━━━━━\n', i);
        
        % 原始磁化率值
        fprintf('原始磁化率：K1=%.6f, K2=%.6f, K3=%.6f\n', ...
            analyzer.data.K1(i), analyzer.data.K2(i), analyzer.data.K3(i));
        
        % Er 張量
        fprintf('\nEr_%d (主軸座標系應變張量):\n', i);
        disp(analyzer.ErList(:,:,i));
        
        % V 矩陣
        fprintf('V_%d (特徵向量矩陣):\n', i);
        disp(analyzer.VList(:,:,i));
        
        % Eg 張量
        fprintf('Eg_%d (地理座標系應變張量):\n', i);
        disp(analyzer.Eg(:,:,i));
        
        % 特徵值分析
        [V_eigen, D_eigen] = analyzer.eigSorted(analyzer.Eg(:,:,i));
        eigenvals = diag(D_eigen);
        
        fprintf('特徵值（由大到小）：[%.6f, %.6f, %.6f]\n', eigenvals);
        fprintf('應變比：%.4f\n', eigenvals(1)/eigenvals(3));
        
        % 主應變軸方向
        fprintf('最大主應變軸方向：\n');
        trend = atan2d(V_eigen(2,1), V_eigen(1,1));
        plunge = asind(V_eigen(3,1));
        fprintf('  方位角：%.2f°, 傾角：%.2f°\n', trend, plunge);
        
        fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
    end
end

function quickAnalysis()
    % 快速分析模式（不需要互動）
    fprintf('\n=== 快速分析模式 ===\n');
    
    % 檔案選擇
    [filename, pathname] = uigetfile({'*.xlsx;*.csv;*.txt'}, '請選擇 AMS 數據檔');
    if isequal(filename, 0)
        return;
    end

    % 建立並執行分析
    analyzer = AMSStrainAnalyzer(fullfile(pathname, filename), 'Verbose', false);
    analyzer.computeFiniteStrainTensors();
    assignSampleResultsToWorkspace(analyzer, '3D');
    
    % 自動匯出結果
    [~, baseName, ~] = fileparts(filename);
    outputFile = sprintf('%s_應變分析結果.xlsx', baseName);
    analyzer.exportToExcel(fullfile(pathname, outputFile));
    
    % 繪製演化圖
    analyzer.plotStrainEvolution();
    
    % 儲存分析物件
    assignin('base', 'ams_analyzer', analyzer);
    
    fprintf('✅ 快速分析完成，結果已匯出至 %s\n', outputFile);
end

function batchAnalysis()
    % 批次分析多個檔案
    fprintf('\n=== 批次分析模式 ===\n');
    
    % 選擇多個檔案
    [filenames, pathname] = uigetfile({'*.xlsx;*.csv;*.txt'}, '請選擇 AMS 數據檔', 'MultiSelect', 'on');
    if isequal(filenames, 0)
        return;
    end
    
    if ~iscell(filenames)
        filenames = {filenames};
    end
    
    fprintf('將分析 %d 個檔案...\n', length(filenames));
    
    % 建立結果儲存結構
    batchResults = struct();
    
    for i = 1:length(filenames)
        try
            fprintf('\n處理檔案 %d/%d: %s\n', i, length(filenames), filenames{i});
            
            % 建立分析器
            analyzer = AMSStrainAnalyzer(fullfile(pathname, filenames{i}), 'Verbose', false);
            analyzer.computeFiniteStrainTensors();
            
            % 儲存結果
            [~, baseName, ~] = fileparts(filenames{i});
            batchResults.(sprintf('file_%d_%s', i, baseName)) = analyzer;
            
            % 匯出個別結果
            outputFile = sprintf('%s_應變分析結果.xlsx', baseName);
            analyzer.exportToExcel(fullfile(pathname, outputFile));
            
        catch ME
            warning('檔案 %s 分析失敗：%s', filenames{i}, ME.message);
        end
    end
    
    % 儲存批次結果
    assignin('base', 'batch_results', batchResults);
    fprintf('\n✅ 批次分析完成，結果儲存在 batch_results 變數中\n');
end

function demonstrateIncrementalAnalysis(analyzer)
    % 示範增量應變分析
    fprintf('\n=== 增量應變分析示範 ===\n');
    
    n = size(analyzer.Eg, 3);
    if n < 2
        fprintf('需要至少 2 個樣本才能進行增量應變分析\n');
        return;
    end
    
    % 示範：計算第1個樣本到最後一個樣本的增量應變
    analyzer.computeIncrementalStrain(1, n, 'demo_increment');
    
    % 顯示結果
    result = analyzer.incrementalResults.demo_increment;
    fprintf('增量變形梯度張量：\n');
    disp(result.F_increment);
    fprintf('對稱拉伸張量：\n');
    disp(result.U);
    fprintf('旋轉張量：\n');
    disp(result.R);
    fprintf('主拉伸值：[%.4f, %.4f, %.4f]\n', result.eigenvalues);
end

function runMode = selectRunMode()
    % 選擇執行模式
    modes = {'標準分析模式', '快速分析模式', '批次分析模式', '取消'};
    [idx, tf] = listdlg('PromptString', '請選擇執行模式：', ...
                        'SelectionMode', 'single', ...
                        'ListString', modes, ...
                        'ListSize', [200, 100]);
    
    if ~tf || idx == 4
        runMode = 'cancel';
    else
        runModes = {'standard', 'quick', 'batch'};
        runMode = runModes{idx};
    end
end

% 主執行邏輯
if ~exist('runMode', 'var')
    runMode = selectRunMode();
end

switch runMode
    case 'standard'
        main();
    case 'quick'
        quickAnalysis();
    case 'batch'
        batchAnalysis();
    case 'cancel'
        disp('程式結束。');
    otherwise
        main();  % 預設執行標準模式
end