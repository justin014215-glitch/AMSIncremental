%{
function main()
    % 主要執行函數
    clear; clc; close all;

    % 設定參數
    %filename = '/Users/lienrueiyu/Desktop/LAB/03collegeproject/01Matlab/file/NXHAMS.xlsx';
    [filename, pathname] = uigetfile('*.xlsx', '請選擇 AMS 數據檔');
    if isequal(filename, 0)
    disp('使用者取消選擇。');
    return;
    end

    fullpath = fullfile(pathname, filename);
    analyzer = AMSStrainAnalyzer(fullpath, ...
    'IntervalCount', 4, ...
    'SamplesPerInterval', 10, ...
    'Verbose', true);

    % 創建分析器物件
    %{
    analyzer = AMSStrainAnalyzer(filename, ...
        'IntervalCount', 4, ...
        'SamplesPerInterval', 10, ...
        'Verbose', true); %}
    %}
    % 執行完整分析
    fprintf('=== 開始北橫磁感率異向性應變增量分析 ===\n');

    % 可選擇分析模式
    choice = input('選擇分析模式 (1: 自動等分, 2: 自定義區間): ');

    switch choice
        case 1
            analyzer.runCompleteAnalysis('auto');
        case 2
            analyzer.runCompleteAnalysis('custom');
        otherwise
            fprintf('無效選擇，使用自動模式\n');
            analyzer.runCompleteAnalysis('auto');
    end

    % 儲存分析器物件供後續使用
    assignin('base', 'ams_analyzer', analyzer);
    fprintf('\n分析器物件已儲存至工作區變數 "ams_analyzer"\n');
end

% 如果尚未定義分析器，則執行主程式
if ~exist('ams_analyzer', 'var')
    main();
end
%}

%{
function main()
    % 使用者輸入檔案名稱
    [filename, pathname] = uigetfile('*.xlsx', '請選擇 AMS 數據檔');
    if isequal(filename, 0)
        disp('使用者取消選擇。');
        return;
    end

    % 建立分析物件
    analyzer = AMSStrainAnalyzer(fullfile(pathname, filename), ...
        'SlateCoeffA', 6.897, ...
        'SlateCoeffB', 0.007, ...
        'Verbose', true);
       %'IntervalCount', 4, ...
       %'SamplesPerInterval', 10, ...
        

    % 選擇分析功能
    while true
        disp('--- 分析功能選單 ---');
        disp('1. 計算 Eg（Er+V+Eg）並顯示前幾筆中間值');
        disp('2. 顯示任一樣本的中間值');
        disp('0. 結束');
        choice = input('請輸入功能編號：');

        switch choice
            case 1
                analyzer.computeFiniteStrainTensors();
                disp('✅ 完整應變張量 Eg 已計算完畢');
                for i = 1:min(3, size(analyzer.Eg,3))
                    fprintf('\n樣本 %d:\n', i);

                    disp('Er ='); disp(analyzer.ErList(:,:,i));
                    disp('V  ='); disp(analyzer.VList(:,:,i));
                    disp('Eg ='); disp(analyzer.Eg(:,:,i));
                end
            case 2
                idx = input('請輸入要檢查的樣本編號：');
                analyzer.displayIntermediateResults(idx);
            case 0
                disp('✅ 分析結束，感謝使用！');
                break;
            otherwise
                disp('請輸入正確選項');
        end
    end

end
%}

%{
function main()
    % 使用者選擇 AMS 數據檔案（xlsx）
    [filename, pathname] = uigetfile('*.xlsx', '請選擇 AMS 數據檔');
    if isequal(filename, 0)
        disp('使用者取消選擇。');
        return;
    end

    % 建立分析物件，並傳入參數
    analyzer = AMSStrainAnalyzer(fullfile(pathname, filename), ...
        'SlateCoeffA', 6.897, ...
        'SlateCoeffB', 0.007, ...
        'Verbose', true);

    % 功能選單迴圈
    while true
        disp('--- 分析功能選單 ---');
        disp('1. 計算 3D Eg (Er + V + Eg)');
        disp('2. 計算 2D Eg (Er + V + Eg)');
        disp('3. 顯示 3D 樣本中間結果 (Er, V, Eg)');
        disp('4. 顯示 2D 樣本中間結果 (Er, V, Eg)');
        disp('0. 結束');
        choice = input('請輸入功能編號：');

        switch choice
            case 1
                analyzer.computeFiniteStrainTensors(); % 3D 計算
                disp('✅ 3D 應變張量 Eg 已計算完畢');
                for i = 1:min(3, size(analyzer.Eg,3))
                    fprintf('\n3D 樣本 %d 中間結果:\n', i);
                    disp('Er ='); disp(analyzer.ErList(:,:,i));
                    disp('V  ='); disp(analyzer.VList(:,:,i));
                    disp('Eg ='); disp(analyzer.Eg(:,:,i));
                end
            case 2
                analyzer.computeFiniteStrainTensors2D(); % 2D 計算
                disp('✅ 2D 應變張量 Eg 已計算完畢');
                for i = 1:min(3, size(analyzer.EgD,3))
                    fprintf('\n2D 樣本 %d 中間結果:\n', i);
                    disp('Er ='); disp(analyzer.ErList2D(:,:,i));
                    disp('V  ='); disp(analyzer.VList2D(:,:,i));
                    disp('Eg ='); disp(analyzer.EgD(:,:,i));
                end
            case 3
                idx = input('請輸入要檢查的 3D 樣本編號：');
                if isempty(analyzer.Eg)
                    disp('請先執行 3D Eg 計算');
                else
                    fprintf('\n3D 樣本 %d 中間結果:\n', idx);
                    disp('Er ='); disp(analyzer.ErList(:,:,idx));
                    disp('V  ='); disp(analyzer.VList(:,:,idx));
                    disp('Eg ='); disp(analyzer.Eg(:,:,idx));
                end
            case 4
                idx = input('請輸入要檢查的 2D 樣本編號：');
                if isempty(analyzer.EgD)
                    disp('請先執行 2D Eg 計算');
                else
                    fprintf('\n2D 樣本 %d 中間結果:\n', idx);
                    disp('Er ='); disp(analyzer.ErList2D(:,:,idx));
                    disp('V  ='); disp(analyzer.VList2D(:,:,idx));
                    disp('Eg ='); disp(analyzer.EgD(:,:,idx));
                end
            case 0
                disp('✅ 分析結束，感謝使用！');
                break;
            otherwise
                disp('請輸入正確選項');
        end
    end
end
%}
%{
% ==========================
% 主程式：AMS 應變分析器
% ==========================
function main()
    % 使用者選擇 AMS 數據檔案（xlsx）
    [filename, pathname] = uigetfile({'*.xlsx;*.csv'}, '請選擇 AMS 數據檔');
    if isequal(filename, 0)
        disp('使用者取消選擇。');
        return;
    end

    % 建立分析物件，並傳入參數
    analyzer = AMSStrainAnalyzer(fullfile(pathname, filename), ...
        'SlateCoeffA', 6.897, ...
        'SlateCoeffB', 0.007, ...
        'Verbose', true);

    fprintf('\n=== 北橫磁感率異向性應變增量分析啟動 ===\n');
    %計算 Er_raw
    analyzer.computeErRawAll();
    % 計算 Eg_raw（使用未轉換的磁感率 Er_raw）
    analyzer.computeEgFromErRaw();

    

    % 計算 3D Eg 並儲存變數
    analyzer.computeFiniteStrainTensors();
    assignSampleResultsToWorkspace(analyzer, '3D');

    % 計算 2D Eg 並儲存變數
    analyzer.computeFiniteStrainTensors2D();
    assignSampleResultsToWorkspace(analyzer, '2D');

    fprintf('\n✅ 所有樣本應變結果已儲存至工作區，可直接使用。\n');
    fprintf('🔍 範例：輸入 `Er_1`、`V_2D_3` 或 `Eg_5` 存取對應結果。\n');

    % 儲存整個物件以供進一步使用
    assignin('base', 'ams_analyzer', analyzer);

    % 啟動選單模式（可選）
    ask = input('\n是否要進入互動選單模式？(y/n)：', 's');
    if strcmpi(ask, 'y')
        interactiveMenu(analyzer);
    else
        disp('✅ 分析已完成。可在命令列輸入參數名稱自由使用結果。');
    end
end


% =========================================
% 將分析結果逐樣本命名儲存至 Workspace
% =========================================
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


% =============================
% 選單功能（可選啟動）
% =============================
function interactiveMenu(analyzer)
    while true
        disp('\n--- 分析功能選單 ---');
        disp('1. 顯示前 3 筆 3D 結果');
        disp('2. 顯示前 3 筆 2D 結果');
        disp('3. 顯示任一筆樣本 3D 結果');
        disp('4. 顯示任一筆樣本 2D 結果');
        disp('0. 離開');
        choice = input('請輸入功能編號：');

        switch choice
            case 1
                for i = 1:min(3, size(analyzer.Eg,3))
                    fprintf('\n樣本 %d（3D）:\n', i);
                    disp(['Er_' num2str(i) ' =']); disp(analyzer.ErList(:,:,i));
                    disp(['V_' num2str(i)  ' =']); disp(analyzer.VList(:,:,i));
                    disp(['Eg_' num2str(i) ' =']); disp(analyzer.Eg(:,:,i));
                end
            case 2
                for i = 1:min(3, size(analyzer.EgD,3))
                    fprintf('\n樣本 %d（2D）:\n', i);
                    disp(['Er_2D_' num2str(i) ' =']); disp(analyzer.ErList2D(:,:,i));
                    disp(['V_2D_' num2str(i)  ' =']); disp(analyzer.VList2D(:,:,i));
                    disp(['Eg_2D_' num2str(i) ' =']); disp(analyzer.EgD(:,:,i));
                end
            case 3
                idx = input('請輸入樣本編號：');
                fprintf('\n樣本 %d（3D）:\n', idx);
                disp(['Er_' num2str(idx) ' =']); disp(analyzer.ErList(:,:,idx));
                disp(['V_' num2str(idx)  ' =']); disp(analyzer.VList(:,:,idx));
                disp(['Eg_' num2str(idx) ' =']); disp(analyzer.Eg(:,:,idx));
            case 4
                idx = input('請輸入樣本編號：');
                fprintf('\n樣本 %d（2D）:\n', idx);
                disp(['Er_2D_' num2str(idx) ' =']); disp(analyzer.ErList2D(:,:,idx));
                disp(['V_2D_' num2str(idx)  ' =']); disp(analyzer.VList2D(:,:,idx));
                disp(['Eg_2D_' num2str(idx) ' =']); disp(analyzer.EgD(:,:,idx));
            case 0
                disp('✅ 離開選單，分析結束！');
                break;
            otherwise
                disp('請輸入有效選項。');
        end
    end
end
%}
%{
classdef AMSStrainAnalyzer < handle
    properties
        filename
        data
        Er
        V
        Eg
        config
        results
    end

    methods
        function obj = AMSStrainAnalyzer(filename, varargin)
            obj.filename = filename;
            obj.config = obj.parseConfig(varargin{:});
            obj.loadData();
        end

        function config = parseConfig(obj, varargin)
            p = inputParser;
            addParameter(p, 'SlateCoeffA', 6.897, @isnumeric);
            addParameter(p, 'SlateCoeffB', 0.007, @isnumeric);
            addParameter(p, 'IntervalCount', 4, @isnumeric);
            addParameter(p, 'SamplesPerInterval', 10, @isnumeric);
            addParameter(p, 'OutputPrefix', '', @ischar);
            addParameter(p, 'Verbose', true, @islogical);
            parse(p, varargin{:});
            config = p.Results;

            if isempty(config.OutputPrefix)
                [~, name, ~] = fileparts(obj.filename);
                config.OutputPrefix = name;
            end
        end

        function loadData(obj)
            try
                obj.data = readtable(obj.filename, 'VariableNamingRule', 'preserve');
                if obj.config.Verbose
                    fprintf('成功載入 %d 筆AMS數據\n', height(obj.data));
                end
            catch ME
                error('載入數據失敗：%s', ME.message);
            end

            required_fields = {'K1', 'K2', 'K3', 'dK1geo', 'iK1geo', 'dK2geo', 'iK2geo', 'dK3geo', 'iK3geo'};
            missing_fields = setdiff(required_fields, obj.data.Properties.VariableNames);
            if ~isempty(missing_fields)
                error('缺少必要欄位：%s', strjoin(missing_fields, ', '));
            end
        end
  function computeFiniteStrainTensors(obj)
            K1 = obj.data.K1; K2 = obj.data.K2; K3 = obj.data.K3;
            dK1 = deg2rad(obj.data.dK1geo); iK1 = deg2rad(obj.data.iK1geo);
            dK2 = deg2rad(obj.data.dK2geo); iK2 = deg2rad(obj.data.iK2geo);
            dK3 = deg2rad(obj.data.dK3geo); iK3 = deg2rad(obj.data.iK3geo);

            num_samples = length(K1);
            obj.Eg = zeros(3, 3, num_samples);

            for i = 1:num_samples
                Er = obj.computeEr(K1(i), K2(i), K3(i));
                V = obj.computeV(dK1(i), iK1(i), dK2(i), iK2(i), dK3(i), iK3(i));
                obj.Eg(:,:,i) = obj.computeEg(Er, V);
            end
        end

        function Er = computeEr(obj, K1, K2, K3)
            K0 = (K1 * K2 * K3)^(1/3);
            a = obj.config.SlateCoeffA;
            b = obj.config.SlateCoeffB;
            ln1pe1 = a * ((K1 / K0) - 1) - b;
            ln1pe2 = a * ((K2 / K0) - 1) - b;
            ln1pe3 = a * ((K3 / K0) - 1) - b;
            e1 = exp(ln1pe1) - 1;
            e2 = exp(ln1pe2) - 1;
            e3 = exp(ln1pe3) - 1;
            omega = (1 + e1) * (1 + e2) * (1 + e3);
            Er = omega^2 * diag([(1+e1)^(-2), (1+e2)^(-2), (1+e3)^(-2)]);
        end

        function V = computeV(~, dK1, iK1, dK2, iK2, dK3, iK3)
            V = [cos(iK1)*cos(dK1), cos(iK2)*cos(dK2), cos(iK3)*cos(dK3);
                 cos(iK1)*sin(dK1), cos(iK2)*sin(dK2), cos(iK3)*sin(dK3);
                 sin(iK1),          sin(iK2),          sin(iK3)];
        end

        function Eg_tensor = computeEg(~, Er, V)
            Eg_tensor = V' * Er * V;
        end










%{
        function computeFiniteStrainTensors(obj)
            K1 = obj.data.K1; K2 = obj.data.K2; K3 = obj.data.K3;
            dK1 = deg2rad(obj.data.dK1geo); iK1 = deg2rad(obj.data.iK1geo);
            dK2 = deg2rad(obj.data.dK2geo); iK2 = deg2rad(obj.data.iK2geo);
            dK3 = deg2rad(obj.data.dK3geo); iK3 = deg2rad(obj.data.iK3geo);

            num_samples = length(K1);
            K0 = (K1 .* K2 .* K3).^(1/3);
            a = obj.config.SlateCoeffA; b = obj.config.SlateCoeffB;

            ln1pe1 = a .* ((K1 ./ K0) - 1) - b;
            ln1pe2 = a .* ((K2 ./ K0) - 1) - b;
            ln1pe3 = a .* ((K3 ./ K0) - 1) - b;
            e1 = exp(ln1pe1) - 1;
            e2 = exp(ln1pe2) - 1;
            e3 = exp(ln1pe3) - 1;
            omega = (1 + e1) .* (1 + e2) .* (1 + e3);

            obj.Eg = zeros(3, 3, num_samples);
            for i = 1:num_samples
                Er = omega(i)^2 * diag([(1+e1(i))^(-2), (1+e2(i))^(-2), (1+e3(i))^(-2)]);
                V = [cos(iK1(i))*cos(dK1(i)), cos(iK2(i))*cos(dK2(i)), cos(iK3(i))*cos(dK3(i));
                     cos(iK1(i))*sin(dK1(i)), cos(iK2(i))*sin(dK2(i)), cos(iK3(i))*sin(dK3(i));
                     sin(iK1(i)), sin(iK2(i)), sin(iK3(i))];
                obj.Eg(:,:,i) = V' * Er * V;
            end
        end
%}
        function performIntervalAnalysis(obj, mode)
            if nargin < 2
                mode = 'auto';
            end
            switch lower(mode)
                case 'auto'
                    obj.autoIntervalAnalysis();
                case 'custom'
                    obj.customIntervalAnalysis();
                otherwise
                    error('未知分析模式：%s', mode);
            end
        end

        function autoIntervalAnalysis(obj)
            total_samples = size(obj.Eg, 3);
            interval_count = obj.config.IntervalCount;
            samples_per_interval = obj.config.SamplesPerInterval;
            Eg_avg = zeros(3, 3, interval_count);
            sample_ranges = zeros(interval_count, 2);

            for i = 1:interval_count
                start_idx = (i-1)*samples_per_interval + 1;
                end_idx = min(i*samples_per_interval, total_samples);
                Eg_avg(:,:,i) = mean(obj.Eg(:,:,start_idx:end_idx), 3);
                sample_ranges(i,:) = [start_idx, end_idx];
            end

            obj.computeIncrementalStrain(Eg_avg, sample_ranges, 'auto');
        end

        function customIntervalAnalysis(obj)
            total_samples = size(obj.Eg, 3);
            fprintf('總樣本數：%d\n', total_samples);
            num_intervals = input('請輸入自訂區間數量：');
            Eg_avg = zeros(3, 3, num_intervals);
            sample_ranges = zeros(num_intervals, 2);

            for i = 1:num_intervals
                while true
                    start_idx = input(sprintf('區間 %d 起始樣本編號：', i));
                    end_idx = input('結束樣本編號：');
                    if start_idx >= 1 && end_idx <= total_samples && start_idx <= end_idx
                        break;
                    else
                        fprintf('樣本範圍錯誤，請重新輸入！\n');
                    end
                end
                sample_ranges(i,:) = [start_idx, end_idx];
                Eg_avg(:,:,i) = mean(obj.Eg(:,:,start_idx:end_idx), 3);
            end

            obj.computeIncrementalStrain(Eg_avg, sample_ranges, 'custom');
        end

        function computeIncrementalStrain(obj, Eg_avg, sample_ranges, mode)
            interval_count = size(Eg_avg, 3);
            obj.results = struct();
            obj.results.mode = mode;
            obj.results.sample_ranges = sample_ranges;
            obj.results.Eg_avg = Eg_avg;
            obj.results.Einc = zeros(3, 3, interval_count-1);
            obj.results.U = zeros(3, 3, interval_count-1);
            obj.results.R = zeros(3, 3, interval_count-1);
            obj.results.eigvals_U = zeros(3, interval_count-1);
            obj.results.eigvecs_U = zeros(3, 3, interval_count-1);
            obj.results.strain_ratios = zeros(2, interval_count-1);
            obj.results.total_strain = zeros(1, interval_count-1);

            for i = 1:interval_count-1
                E1 = Eg_avg(:,:,i);
                E2 = Eg_avg(:,:,i+1);
                F = E2 / E1;
                obj.results.Einc(:,:,i) = F;
                C = F' * F;
                U = sqrtm(C);
                obj.results.U(:,:,i) = U;
                obj.results.R(:,:,i) = F / U;
                [eigvec, eigval] = eig(U);
                [sorted_vals, idx] = sort(diag(eigval), 'descend');
                obj.results.eigvecs_U(:,:,i) = eigvec(:,idx);
                obj.results.eigvals_U(:,i) = sorted_vals;
                obj.results.strain_ratios(1,i) = sorted_vals(1) / sorted_vals(2);
                obj.results.strain_ratios(2,i) = sorted_vals(2) / sorted_vals(3);
                obj.results.total_strain(i) = sqrt(sum((sorted_vals - 1).^2));
            end
        end

        function displayResults(obj)
            if isempty(obj.results)
                fprintf('尚未進行分析，請先執行 performIntervalAnalysis\n');
                return;
            end
            for i = 1:size(obj.results.eigvals_U, 2)
                fprintf('\n--- 區間 %d → %d ---\n', i, i+1);
                fprintf('主應變值: [%.4f, %.4f, %.4f]\n', obj.results.eigvals_U(1,i), obj.results.eigvals_U(2,i), obj.results.eigvals_U(3,i));
                fprintf('應變比例 (L, F): [%.4f, %.4f]\n', obj.results.strain_ratios(1,i), obj.results.strain_ratios(2,i));
                fprintf('總應變強度: %.4f\n', obj.results.total_strain(i));
            end
        end

        function plotResults(obj)
            if isempty(obj.results)
                fprintf('尚未進行分析\n'); return;
            end
            figure;
            subplot(1,2,1);
            bar(obj.results.eigvals_U');
            title('主應變值'); xlabel('區間'); ylabel('U 值');
            legend('U1','U2','U3');
            subplot(1,2,2);
            plot(obj.results.strain_ratios(1,:), '-o'); hold on;
            plot(obj.results.strain_ratios(2,:), '-s');
            title('應變比例'); xlabel('區間'); ylabel('比例');
            legend('L', 'F');
        end

        function plotStrainEllipsoid(obj, interval_idx)
            if interval_idx > size(obj.results.eigvals_U,2)
                return;
            end
            vals = obj.results.eigvals_U(:, interval_idx);
            vecs = obj.results.eigvecs_U(:,:,interval_idx);
            [x, y, z] = sphere(20);
            points = [x(:), y(:), z(:)];
            ellipsoid_points = points * diag(vals) * vecs';
            xe = reshape(ellipsoid_points(:,1), size(x));
            ye = reshape(ellipsoid_points(:,2), size(y));
            ze = reshape(ellipsoid_points(:,3), size(z));
            figure; surf(xe, ye, ze, 'FaceAlpha', 0.7);
            axis equal; xlabel('X'); ylabel('Y'); zlabel('Z'); title('應變橢球體');
        end

        function exportResults(obj)
            if isempty(obj.results)
                return;
            end
            T = table();
            n_intervals = size(obj.results.eigvals_U, 2);
            T.From = (1:n_intervals)';
            T.To = (2:n_intervals+1)';
            T.U1 = obj.results.eigvals_U(1,:)';
            T.U2 = obj.results.eigvals_U(2,:)';
            T.U3 = obj.results.eigvals_U(3,:)';
            T.L_ratio = obj.results.strain_ratios(1,:)';
            T.F_ratio = obj.results.strain_ratios(2,:)';
            T.TotalStrain = obj.results.total_strain';
            filename = sprintf('%s_incremental_analysis.xlsx', obj.config.OutputPrefix);
            writetable(T, filename);
            fprintf('已匯出結果至 %s\n', filename);
        end

        function runCompleteAnalysis(obj, mode)
            if nargin < 2
                mode = 'auto';
            end
            try
                obj.computeFiniteStrainTensors();
                obj.performIntervalAnalysis(mode);
                obj.displayResults();
                obj.plotResults();
                obj.exportResults();
                fprintf('\n✅ 完整分析流程執行完畢！\n');
            catch ME
                fprintf('❌ 分析過程中出現錯誤：%s\n', ME.message);
                rethrow(ME);
            end
        end
    end
end
%}
%{
classdef AMSStrainAnalyzer < handle
    properties
        filename
        data
        Er
        V
        Eg
        config
        results
    end

    methods
        function obj = AMSStrainAnalyzer(filename, varargin)
            obj.filename = filename;
            obj.config = obj.parseConfig(varargin{:});
            obj.loadData();
        end

        function config = parseConfig(obj, varargin)
            p = inputParser;
            addParameter(p, 'SlateCoeffA', 6.897, @isnumeric);
            addParameter(p, 'SlateCoeffB', 0.007, @isnumeric);
            addParameter(p, 'IntervalCount', 4, @isnumeric);
            addParameter(p, 'SamplesPerInterval', 10, @isnumeric);
            addParameter(p, 'OutputPrefix', '', @ischar);
            addParameter(p, 'Verbose', true, @islogical);
            parse(p, varargin{:});
            config = p.Results;

            if isempty(config.OutputPrefix)
                [~, name, ~] = fileparts(obj.filename);
                config.OutputPrefix = name;
            end
        end

        function loadData(obj)
            try
                obj.data = readtable(obj.filename, 'VariableNamingRule', 'preserve');
                if obj.config.Verbose
                    fprintf('成功載入 %d 筆AMS數據\n', height(obj.data));
                end
            catch ME
                error('載入數據失敗：%s', ME.message);
            end

            required_fields = {'K1', 'K2', 'K3', 'dK1geo', 'iK1geo', 'dK2geo', 'iK2geo', 'dK3geo', 'iK3geo'};
            missing_fields = setdiff(required_fields, obj.data.Properties.VariableNames);
            if ~isempty(missing_fields)
                error('缺少必要欄位：%s', strjoin(missing_fields, ', '));
            end
        end
  function computeFiniteStrainTensors(obj)
            K1 = obj.data.K1; K2 = obj.data.K2; K3 = obj.data.K3;
            dK1 = deg2rad(obj.data.dK1geo); iK1 = deg2rad(obj.data.iK1geo);
            dK2 = deg2rad(obj.data.dK2geo); iK2 = deg2rad(obj.data.iK2geo);
            dK3 = deg2rad(obj.data.dK3geo); iK3 = deg2rad(obj.data.iK3geo);

            num_samples = length(K1);
            obj.Eg = zeros(3, 3, num_samples);

            for i = 1:num_samples
                Er = obj.computeEr(K1(i), K2(i), K3(i));
                V = obj.computeV(dK1(i), iK1(i), dK2(i), iK2(i), dK3(i), iK3(i));
                obj.Eg(:,:,i) = obj.computeEg(Er, V);
            end
        end

        function Er = computeEr(obj, K1, K2, K3)
            K0 = (K1 * K2 * K3)^(1/3);
            a = obj.config.SlateCoeffA;
            b = obj.config.SlateCoeffB;
            ln1pe1 = a * ((K1 / K0) - 1) - b;
            ln1pe2 = a * ((K2 / K0) - 1) - b;
            ln1pe3 = a * ((K3 / K0) - 1) - b;
            e1 = exp(ln1pe1) - 1;
            e2 = exp(ln1pe2) - 1;
            e3 = exp(ln1pe3) - 1;
            omega = (1 + e1) * (1 + e2) * (1 + e3);
            Er = omega^2 * diag([(1+e1)^(-2), (1+e2)^(-2), (1+e3)^(-2)]);
        end

        function V = computeV(~, dK1, iK1, dK2, iK2, dK3, iK3)
            V = [cos(iK1)*cos(dK1), cos(iK2)*cos(dK2), cos(iK3)*cos(dK3);
                 cos(iK1)*sin(dK1), cos(iK2)*sin(dK2), cos(iK3)*sin(dK3);
                 sin(iK1),          sin(iK2),          sin(iK3)];
        end

        function Eg_tensor = computeEg(~, Er, V)
            Eg_tensor = V' * Er * V;
        end










%{
        function computeFiniteStrainTensors(obj)
            K1 = obj.data.K1; K2 = obj.data.K2; K3 = obj.data.K3;
            dK1 = deg2rad(obj.data.dK1geo); iK1 = deg2rad(obj.data.iK1geo);
            dK2 = deg2rad(obj.data.dK2geo); iK2 = deg2rad(obj.data.iK2geo);
            dK3 = deg2rad(obj.data.dK3geo); iK3 = deg2rad(obj.data.iK3geo);

            num_samples = length(K1);
            K0 = (K1 .* K2 .* K3).^(1/3);
            a = obj.config.SlateCoeffA; b = obj.config.SlateCoeffB;

            ln1pe1 = a .* ((K1 ./ K0) - 1) - b;
            ln1pe2 = a .* ((K2 ./ K0) - 1) - b;
            ln1pe3 = a .* ((K3 ./ K0) - 1) - b;
            e1 = exp(ln1pe1) - 1;
            e2 = exp(ln1pe2) - 1;
            e3 = exp(ln1pe3) - 1;
            omega = (1 + e1) .* (1 + e2) .* (1 + e3);

            obj.Eg = zeros(3, 3, num_samples);
            for i = 1:num_samples
                Er = omega(i)^2 * diag([(1+e1(i))^(-2), (1+e2(i))^(-2), (1+e3(i))^(-2)]);
                V = [cos(iK1(i))*cos(dK1(i)), cos(iK2(i))*cos(dK2(i)), cos(iK3(i))*cos(dK3(i));
                     cos(iK1(i))*sin(dK1(i)), cos(iK2(i))*sin(dK2(i)), cos(iK3(i))*sin(dK3(i));
                     sin(iK1(i)), sin(iK2(i)), sin(iK3(i))];
                obj.Eg(:,:,i) = V' * Er * V;
            end
        end
%}
        function performIntervalAnalysis(obj, mode)
            if nargin < 2
                mode = 'auto';
            end
            switch lower(mode)
                case 'auto'
                    obj.autoIntervalAnalysis();
                case 'custom'
                    obj.customIntervalAnalysis();
                otherwise
                    error('未知分析模式：%s', mode);
            end
        end

        function autoIntervalAnalysis(obj)
            total_samples = size(obj.Eg, 3);
            interval_count = obj.config.IntervalCount;
            samples_per_interval = obj.config.SamplesPerInterval;
            Eg_avg = zeros(3, 3, interval_count);
            sample_ranges = zeros(interval_count, 2);

            for i = 1:interval_count
                start_idx = (i-1)*samples_per_interval + 1;
                end_idx = min(i*samples_per_interval, total_samples);
                Eg_avg(:,:,i) = mean(obj.Eg(:,:,start_idx:end_idx), 3);
                sample_ranges(i,:) = [start_idx, end_idx];
            end

            obj.computeIncrementalStrain(Eg_avg, sample_ranges, 'auto');
        end

        function customIntervalAnalysis(obj)
            total_samples = size(obj.Eg, 3);
            fprintf('總樣本數：%d\n', total_samples);
            num_intervals = input('請輸入自訂區間數量：');
            Eg_avg = zeros(3, 3, num_intervals);
            sample_ranges = zeros(num_intervals, 2);

            for i = 1:num_intervals
                while true
                    start_idx = input(sprintf('區間 %d 起始樣本編號：', i));
                    end_idx = input('結束樣本編號：');
                    if start_idx >= 1 && end_idx <= total_samples && start_idx <= end_idx
                        break;
                    else
                        fprintf('樣本範圍錯誤，請重新輸入！\n');
                    end
                end
                sample_ranges(i,:) = [start_idx, end_idx];
                Eg_avg(:,:,i) = mean(obj.Eg(:,:,start_idx:end_idx), 3);
            end

            obj.computeIncrementalStrain(Eg_avg, sample_ranges, 'custom');
        end

        function computeIncrementalStrain(obj, Eg_avg, sample_ranges, mode)
            interval_count = size(Eg_avg, 3);
            obj.results = struct();
            obj.results.mode = mode;
            obj.results.sample_ranges = sample_ranges;
            obj.results.Eg_avg = Eg_avg;
            obj.results.Einc = zeros(3, 3, interval_count-1);
            obj.results.U = zeros(3, 3, interval_count-1);
            obj.results.R = zeros(3, 3, interval_count-1);
            obj.results.eigvals_U = zeros(3, interval_count-1);
            obj.results.eigvecs_U = zeros(3, 3, interval_count-1);
            obj.results.strain_ratios = zeros(2, interval_count-1);
            obj.results.total_strain = zeros(1, interval_count-1);

            for i = 1:interval_count-1
                E1 = Eg_avg(:,:,i);
                E2 = Eg_avg(:,:,i+1);
                F = E2 / E1;
                obj.results.Einc(:,:,i) = F;
                C = F' * F;
                U = sqrtm(C);
                obj.results.U(:,:,i) = U;
                obj.results.R(:,:,i) = F / U;
                [eigvec, eigval] = eig(U);
                [sorted_vals, idx] = sort(diag(eigval), 'descend');
                obj.results.eigvecs_U(:,:,i) = eigvec(:,idx);
                obj.results.eigvals_U(:,i) = sorted_vals;
                obj.results.strain_ratios(1,i) = sorted_vals(1) / sorted_vals(2);
                obj.results.strain_ratios(2,i) = sorted_vals(2) / sorted_vals(3);
                obj.results.total_strain(i) = sqrt(sum((sorted_vals - 1).^2));
            end
        end

        function displayResults(obj)
            if isempty(obj.results)
                fprintf('尚未進行分析，請先執行 performIntervalAnalysis\n');
                return;
            end
            for i = 1:size(obj.results.eigvals_U, 2)
                fprintf('\n--- 區間 %d → %d ---\n', i, i+1);
                fprintf('主應變值: [%.4f, %.4f, %.4f]\n', obj.results.eigvals_U(1,i), obj.results.eigvals_U(2,i), obj.results.eigvals_U(3,i));
                fprintf('應變比例 (L, F): [%.4f, %.4f]\n', obj.results.strain_ratios(1,i), obj.results.strain_ratios(2,i));
                fprintf('總應變強度: %.4f\n', obj.results.total_strain(i));
            end
        end

        function plotResults(obj)
            if isempty(obj.results)
                fprintf('尚未進行分析\n'); return;
            end
            figure;
            subplot(1,2,1);
            bar(obj.results.eigvals_U');
            title('主應變值'); xlabel('區間'); ylabel('U 值');
            legend('U1','U2','U3');
            subplot(1,2,2);
            plot(obj.results.strain_ratios(1,:), '-o'); hold on;
            plot(obj.results.strain_ratios(2,:), '-s');
            title('應變比例'); xlabel('區間'); ylabel('比例');
            legend('L', 'F');
        end

        function plotStrainEllipsoid(obj, interval_idx)
            if interval_idx > size(obj.results.eigvals_U,2)
                return;
            end
            vals = obj.results.eigvals_U(:, interval_idx);
            vecs = obj.results.eigvecs_U(:,:,interval_idx);
            [x, y, z] = sphere(20);
            points = [x(:), y(:), z(:)];
            ellipsoid_points = points * diag(vals) * vecs';
            xe = reshape(ellipsoid_points(:,1), size(x));
            ye = reshape(ellipsoid_points(:,2), size(y));
            ze = reshape(ellipsoid_points(:,3), size(z));
            figure; surf(xe, ye, ze, 'FaceAlpha', 0.7);
            axis equal; xlabel('X'); ylabel('Y'); zlabel('Z'); title('應變橢球體');
        end

        function exportResults(obj)
            if isempty(obj.results)
                return;
            end
            T = table();
            n_intervals = size(obj.results.eigvals_U, 2);
            T.From = (1:n_intervals)';
            T.To = (2:n_intervals+1)';
            T.U1 = obj.results.eigvals_U(1,:)';
            T.U2 = obj.results.eigvals_U(2,:)';
            T.U3 = obj.results.eigvals_U(3,:)';
            T.L_ratio = obj.results.strain_ratios(1,:)';
            T.F_ratio = obj.results.strain_ratios(2,:)';
            T.TotalStrain = obj.results.total_strain';
            filename = sprintf('%s_incremental_analysis.xlsx', obj.config.OutputPrefix);
            writetable(T, filename);
            fprintf('已匯出結果至 %s\n', filename);
        end

        function runCompleteAnalysis(obj, mode)
            if nargin < 2
                mode = 'auto';
            end
            try
                obj.computeFiniteStrainTensors();
                obj.performIntervalAnalysis(mode);
                obj.displayResults();
                obj.plotResults();
                obj.exportResults();
                fprintf('\n✅ 完整分析流程執行完畢！\n');
            catch ME
                fprintf('❌ 分析過程中出現錯誤：%s\n', ME.message);
                rethrow(ME);
            end
        end
    end
end
%}
%{
classdef AMSStrainAnalyzer < handle
    properties
        filename
        data
        Eg
        ErList
        VList
        config
        Er_raw       % 儲存原始磁化率張量 Er_raw = diag(K1, K2, K3)
        EgRaw        % 儲存從 Er_raw 推出的 Eg_raw = V' * Er_raw * V

      % 新增2D結果
        EgD
        ErList2D
        VList2D
    end

    methods
        function obj = AMSStrainAnalyzer(filename, varargin)
            obj.filename = filename; %儲存檔案
            obj.config = obj.parseConfig(varargin{:});  %分析數據
            obj.loadData(); %載入資料

        end

        function config = parseConfig(obj, varargin)
            p = inputParser;
            addParameter(p, 'SlateCoeffA', 6.897, @isnumeric);
            addParameter(p, 'SlateCoeffB', 0.007, @isnumeric);
            addParameter(p, 'Verbose', true, @islogical);
            parse(p, varargin{:});
            config = p.Results;
        end

        function loadData(obj)
            try
                obj.data = readtable(obj.filename, 'VariableNamingRule', 'preserve');
                if obj.config.Verbose
                    fprintf('成功載入 %d 筆AMS數據\n', height(obj.data));
                end
            catch ME
                error('載入數據失敗：%s', ME.message);
            end

            required_fields = {'K1', 'K2', 'K3', 'dK1geo', 'iK1geo', 'dK2geo', 'iK2geo', 'dK3geo', 'iK3geo'};
            missing_fields = setdiff(required_fields, obj.data.Properties.VariableNames);
            if ~isempty(missing_fields)
                error('缺少必要欄位：%s', strjoin(missing_fields, ', '));
            end
        end

        function computeFiniteStrainTensors(obj)
            K1 = obj.data.K1; K2 = obj.data.K2; K3 = obj.data.K3;
%{            
            P = obj.data.P; T = obj.data.T; F = obj.data.F; L = obj.data.L; 
            Pj = obj.data.Pj; Int = (sqrt((L-1).^2+(F-1).^2));
%}
            %角度轉弧度
            dK1 = deg2rad(obj.data.dK1geo); iK1 = deg2rad(obj.data.iK1geo);
            dK2 = deg2rad(obj.data.dK2geo); iK2 = deg2rad(obj.data.iK2geo);
            dK3 = deg2rad(obj.data.dK3geo); iK3 = deg2rad(obj.data.iK3geo);
            
            %結果儲存矩陣
            num_samples = length(K1);
            obj.Eg = zeros(3, 3, num_samples);
            obj.ErList = zeros(3, 3, num_samples);
            obj.VList = zeros(3, 3, num_samples);
            

            %計算到第幾筆資料
            if obj.config.Verbose
               fprintf('開始計算 %d 筆有限應變張量...\n', num_samples);
            end  

            %
            for i = 1:num_samples
                Er = obj.computeEr(K1(i), K2(i), K3(i));
                V = obj.computeV(dK1(i), iK1(i), dK2(i), iK2(i), dK3(i), iK3(i));
                Eg = obj.computeEg(Er, V);
                obj.ErList(:,:,i) = Er;
                obj.VList(:,:,i) = V;
                obj.Eg(:,:,i) = Eg;
                if obj.config.Verbose && mod(i, 10) == 0
            fprintf('  → 已完成第 %d 筆樣本\n', i);
                end
            end
        end

        function Er = computeEr(obj, K1, K2, K3)
            K0 = (K1 * K2 * K3)^(1/3);
            a = obj.config.SlateCoeffA;
            b = obj.config.SlateCoeffB;
            ln1pe1 = a * ((K1 / K0) - 1) - b;
            ln1pe2 = a * ((K2 / K0) - 1) - b;
            ln1pe3 = a * ((K3 / K0) - 1) - b;
            e1 = exp(ln1pe1) - 1;
            e2 = exp(ln1pe2) - 1;
            e3 = exp(ln1pe3) - 1;
            omega = (1 + e1) * (1 + e2) * (1 + e3);
            Er = omega^2 * diag([(1+e1)^(-2), (1+e2)^(-2), (1+e3)^(-2)]);
            if obj.config.Verbose
               fprintf('Er = diag([%.4f %.4f %.4f]) × omega² = %.4f\n', ...
               (1+e1)^(-2), (1+e2)^(-2), (1+e3)^(-2), omega^2);
            end
        end

        function V = computeV(obj, dK1, iK1, dK2, iK2, dK3, iK3)
            V = [cos(iK1).*cos(dK1), cos(iK2).*cos(dK2), cos(iK3).*cos(dK3);
                 cos(iK1).*sin(dK1), cos(iK2).*sin(dK2), cos(iK3).*sin(dK3);
                 sin(iK1),          sin(iK2),          sin(iK3)];
            if isfield(obj, 'config') && obj.config.Verbose
               orthogonality = V' * V;
               fprintf('V^T * V =\n');
               disp(orthogonality);
            end
        end

        function Eg_tensor = computeEg( ~, Er, V)
            Eg_tensor = V' * Er * V;
            if norm(Eg_tensor - Eg_tensor') > 1e-6
               warning('Eg_tensor 非對稱，可能有計算誤差');
            end
        end

        function displayIntermediateResults(obj, sample_idx)
            if nargin < 2
                sample_idx = 1;
            end
            fprintf('\n--- 樣本 %d 的中間結果 ---\n', sample_idx);
            disp('Er:'); disp(obj.ErList(:,:,sample_idx));
            disp('V:'); disp(obj.VList(:,:,sample_idx));
            disp('Eg:'); disp(obj.Eg(:,:,sample_idx));
        end
        function computeFiniteStrainTensors2D(obj)
        % 只針對K1、K2計算2D版本 Er、V、Eg
            K1 = obj.data.K1;
            K2 = obj.data.K2;
            dK1 = deg2rad(obj.data.dK1geo);
            iK1 = deg2rad(obj.data.iK1geo);
            dK2 = deg2rad(obj.data.dK2geo);
            iK2 = deg2rad(obj.data.iK2geo);
        
            num_samples = length(K1);
        
            obj.EgD = zeros(2, 2, num_samples);
            obj.ErList2D = zeros(2, 2, num_samples);
            obj.VList2D = zeros(2, 2, num_samples);
        
            a = obj.config.SlateCoeffA;
            b = obj.config.SlateCoeffB;
    
            for i = 1:num_samples
                % 計算2D Er
                K0 = sqrt(K1(i)*K2(i));
                ln1pe1 = a*((K1(i)/K0)-1) - b;
                ln1pe2 = a*((K2(i)/K0)-1) - b;
                e1 = exp(ln1pe1) - 1;
                e2 = exp(ln1pe2) - 1;
                omega = (1 + e1) * (1 + e2);
        
                Er2D = omega^2 * diag([(1+e1)^(-2), (1+e2)^(-2)]);
        
                % 計算2D方向矩陣V，只取水平面投影 (X, Y)
                v1 = [cos(iK1(i))*cos(dK1(i)); cos(iK1(i))*sin(dK1(i))];
                v2 = [cos(iK2(i))*cos(dK2(i)); cos(iK2(i))*sin(dK2(i))];
                V2D = [v1, v2];
        
                % 計算2D Eg

                EgD = V2D' * Er2D * V2D;
        
                % 儲存結果
                obj.ErList2D(:,:,i) = Er2D;
                obj.VList2D(:,:,i) = V2D;
                obj.EgD(:,:,i) = EgD;
            end
        end
        function computeErRawAll(obj)
                K1 = obj.data.K1;
                K2 = obj.data.K2;
                K3 = obj.data.K3;
                n = length(K1);
                obj.Er_raw = zeros(3, 3, n);  % 預留儲存空間
            
                for i = 1:n
                    Er = obj.computeEr_raw(K1(i), K2(i), K3(i), i);
                    obj.Er_raw(:,:,i) = Er;  % 存入物件
                end
        end

         function computeEgFromErRaw(obj)
                % 從原始磁化率張量 (Er_raw = diag[K1, K2, K3]) 計算 Eg
                K1 = obj.data.K1;
                K2 = obj.data.K2;
                K3 = obj.data.K3;
            
                % 方位角與傾角
                dK1 = deg2rad(obj.data.dK1geo); iK1 = deg2rad(obj.data.iK1geo);
                dK2 = deg2rad(obj.data.dK2geo); iK2 = deg2rad(obj.data.iK2geo);
                dK3 = deg2rad(obj.data.dK3geo); iK3 = deg2rad(obj.data.iK3geo);
            
                n = length(K1);
                Eg_raw = zeros(3,3,n); % 儲存 Eg_raw 結果
            
                for i = 1:n
                    % 直接構建對角張量 Er_raw
                    Er_raw = diag([K1(i), K2(i), K3(i)]);
            
                    % 計算方向矩陣
                    V = obj.computeV(dK1(i), iK1(i), dK2(i), iK2(i), dK3(i), iK3(i));
            
                    % Eg = V' * Er_raw * V
                    Eg = V' * Er_raw * V;
                    Eg_raw(:,:,i) = Eg;
            
                    % 顯示中間結果
                    if obj.config.Verbose
                        fprintf('Eg_raw_%d 已計算完成。\n', i);
                    end
            
                    % 儲存至 workspace
                    assignin('base', sprintf('Eg_raw_%d', i), Eg);
                end
            
                % 若需存入物件中，可考慮擴充：
                obj.EgRaw = Eg_raw;
         end


    end
 end
%}
