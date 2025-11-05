classdef AMSStrainAnalyzer < handle
    properties
        filename
        data
        Eg              % 地理座標系中的有限應變張量
        ErList          % 主軸座標系中的應變張量
        VList           % 特徵向量矩陣
        config          % 配置參數
        Er_raw          % 原始磁化率構成的張量 Er_raw = diag(K1, K2, K3)
        EgRaw           % 由 Er_raw 推得的 Eg_raw = V' * Er_raw * V
        
        % 增量應變分析結果
        incrementalResults  % 儲存增量應變分析結果
    end

    
    methods
        function obj = AMSStrainAnalyzer(filename, varargin)
            % 建構函數
            % 輸入：
            %   filename - 資料檔案路徑
            %   varargin - 選用參數（SlateCoeffA, SlateCoeffB, Verbose等）
            
            obj.filename = filename;
            obj.config = obj.parseConfig(varargin{:});
            obj.loadData();
            obj.incrementalResults = struct();
        end

        function config = parseConfig(~, varargin)
            % 解析配置參數
            p = inputParser;
            addParameter(p, 'SlateCoeffA', 6.897, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'SlateCoeffB', 0.007, @isnumeric);
            addParameter(p, 'Verbose', true, @islogical);
            addParameter(p, 'CleanThreshold', 1e-12, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'ValidateResults', true, @islogical);
            addParameter(p, 'RockType', 'slate', @ischar);
            parse(p, varargin{:});
            config = p.Results;
        end

        function loadData(obj)
            % 載入並驗證資料檔案
            try
                [~, ~, ext] = fileparts(obj.filename);
                switch lower(ext)
                    case '.xlsx'
                        obj.data = readtable(obj.filename, 'VariableNamingRule', 'preserve');
                    case '.csv'
                        obj.data = readtable(obj.filename, 'VariableNamingRule', 'preserve');
                    case '.txt'
                        obj.data = readtable(obj.filename, 'VariableNamingRule', 'preserve', 'Delimiter', '\t');
                    otherwise
                        error('不支援的檔案格式：%s。請使用 .xlsx、.csv 或 .txt 格式', ext);
                end
                
                if obj.config.Verbose
                    fprintf('成功載入 %d 筆 AMS 數據（檔案：%s）\n', height(obj.data), obj.filename);
                end
            catch ME
                error('載入數據失敗：%s', ME.message);
            end

            % 檢查必要欄位
            required = {'K1','K2','K3','dK1geo','iK1geo','dK2geo','iK2geo','dK3geo','iK3geo'};
            missing = setdiff(required, obj.data.Properties.VariableNames);
            if ~isempty(missing)
                error('缺少必要欄位：%s', strjoin(missing, ', '));
            end
            
            % 基本資料驗證
            if any(obj.data.K1 <= 0) || any(obj.data.K2 <= 0) || any(obj.data.K3 <= 0)
                error('發現非正數的磁化率值，請檢查資料品質');
            end
            
            % 檢查磁化率順序
            invalidOrder = obj.data.K1 < obj.data.K2 | obj.data.K2 < obj.data.K3;
            if any(invalidOrder)
                warning('發現 %d 筆資料不符合 K1≥K2≥K3 順序，請檢查資料', sum(invalidOrder));
            end
        end

        function computeFiniteStrainTensors(obj)
            % 計算所有樣本的有限應變張量
            % 
            % 此方法對每個樣本：
            % 1. 計算主軸座標系中的應變張量 Er
            % 2. 計算特徵向量矩陣 V  
            % 3. 轉換至地理座標系得到 Eg
            
            K1 = obj.data.K1; K2 = obj.data.K2; K3 = obj.data.K3;
            dK1 = deg2rad(obj.data.dK1geo); iK1 = deg2rad(obj.data.iK1geo);
            dK2 = deg2rad(obj.data.dK2geo); iK2 = deg2rad(obj.data.iK2geo);
            dK3 = deg2rad(obj.data.dK3geo); iK3 = deg2rad(obj.data.iK3geo);

            n = length(K1);
            obj.Eg = zeros(3,3,n);
            obj.ErList = zeros(3,3,n);
            obj.VList = zeros(3,3,n);
            
            if n > 100 && obj.config.Verbose
                fprintf('處理 %d 個樣本，可能需要一些時間...\n', n);
            end

            for i = 1:n
                try
                    Er = obj.computeEr(K1(i), K2(i), K3(i));
                    V = obj.computeV(dK1(i), iK1(i), dK2(i), iK2(i), dK3(i), iK3(i));
                    Eg = obj.computeEg(Er, V);
                    
                    % 驗證結果
                    if obj.config.ValidateResults
                        obj.validateStrainTensor(Eg);
                    end
                    
                    obj.ErList(:,:,i) = Er;
                    obj.VList(:,:,i) = V;
                    obj.Eg(:,:,i) = Eg;
                    
                    if obj.config.Verbose && mod(i, 50) == 0
                        fprintf('已完成 %d/%d 樣本\n', i, n);
                    end
                catch ME
                    warning('樣本 %d 計算失敗：%s', i, ME.message);
                    continue;
                end
            end
            
            if obj.config.Verbose
                fprintf('有限應變張量計算完成\n');
            end
        end

        function Er = computeEr(obj, K1, K2, K3)
            % 計算磁化率應變張量
            % 
            % 輸入參數：
            %   K1, K2, K3 - 主磁化率值
            % 
            % 輸出參數：
            %   Er - 3x3 應變張量
            % 
            % 計算方法：
            %   基於 Slate 係數的指數變換
            %   Er = ω² * diag([(1+e₁)⁻², (1+e₂)⁻², (1+e₃)⁻²])
            
            % 驗證輸入
            if any([K1, K2, K3] <= 0)
                error('磁化率值必須為正數');
            end
            
            K0 = (K1 * K2 * K3)^(1/3);
            a = obj.config.SlateCoeffA;
            b = obj.config.SlateCoeffB;
            
            e1 = exp(a * ((K1/K0)-1) - b) - 1;
            e2 = exp(a * ((K2/K0)-1) - b) - 1;
            e3 = exp(a * ((K3/K0)-1) - b) - 1;
            
            omega = (1+e1)*(1+e2)*(1+e3);
            Er = omega^2 * diag([(1+e1)^-2, (1+e2)^-2, (1+e3)^-2]);
        end

        function V = computeV(obj, d1, i1, d2, i2, d3, i3)
            % 計算特徵向量矩陣
            %
            % 輸入參數：trend & plunge (dK1geo, dK2geo, dK3geo, iK1geo, iK2geo, iK3geo)
            %
            % 輸出參數:
            %   V - 特徵向量矩陣
            %
            % 計算方法：
            %   V = [cos(i1)*cos(d1), cos(i2)*cos(d2), cos(i3)*cos(d3); %N = X
            %        cos(i1)*sin(d1), cos(i2)*sin(d2), cos(i3)*sin(d3); %E = Y
            %        sin(i1),         sin(i2),         sin(i3)];        %D = Z
            
            V = [cos(i1)*cos(d1), cos(i2)*cos(d2), cos(i3)*cos(d3); %N = X
                 cos(i1)*sin(d1), cos(i2)*sin(d2), cos(i3)*sin(d3); %E = Y
                 sin(i1),         sin(i2),         sin(i3)];        %D = Z
            
            
%{
            % 檢查矩陣條件
            if abs(det(V)) < 1e-10
                warning('特徵向量矩陣接近奇異，可能影響計算精度');
            end
            %}            
            % 清理極小值
            threshold = obj.config.CleanThreshold;
            V(abs(V) < threshold) = 0;
        end

        function Eg = computeEg(~, Er, V)
            % 計算地理座標系中的應變張量
            %
            % 輸入參數：V, Er
            % 計算方法：Eg = V' * Er * V
            
            Eg = V' * Er * V;
        end

        function computeErRawAll(obj)
            % 計算原始磁化率張量（驗證用，不做磁感率轉應變）
            K1 = obj.data.K1; K2 = obj.data.K2; K3 = obj.data.K3;
            n = length(K1);
            obj.Er_raw = zeros(3,3,n);

            for i = 1:n
                Er = diag([K1(i), K2(i), K3(i)]);
                obj.Er_raw(:,:,i) = Er;
                assignin('base', sprintf('Er_raw_%d', i), Er);
                if obj.config.Verbose
                    fprintf('Er_raw_%d = diag([%.4e, %.4e, %.4e])\n', i, K1(i), K2(i), K3(i));
                end
            end
        end

        function computeEgFromErRaw(obj)
            % 從原始磁化率張量計算 Eg（驗證用）
            if isempty(obj.Er_raw)
                obj.computeErRawAll();
            end
            
            dK1 = deg2rad(obj.data.dK1geo); iK1 = deg2rad(obj.data.iK1geo);
            dK2 = deg2rad(obj.data.dK2geo); iK2 = deg2rad(obj.data.iK2geo);
            dK3 = deg2rad(obj.data.dK3geo); iK3 = deg2rad(obj.data.iK3geo);
            n = size(obj.Er_raw, 3);
            obj.EgRaw = zeros(3,3,n);

            for i = 1:n
                V = obj.computeV(dK1(i), iK1(i), dK2(i), iK2(i), dK3(i), iK3(i));
                Eg = V' * obj.Er_raw(:,:,i) * V;
                obj.EgRaw(:,:,i) = Eg;
                assignin('base', sprintf('Eg_raw_%d', i), Eg);
                if obj.config.Verbose
                    fprintf('Eg_raw_%d 已完成\n', i);
                end
            end
        end

        function [V_sorted, D_sorted] = eigSorted(~, A)
            % 對矩陣 A 做特徵分解並依特徵值由大到小排序
            %
            % 輸入：
            %   A - 方陣
            %
            % 輸出：
            %   V_sorted - 排序後的特徵向量矩陣（每一列是特徵向量）
            %   D_sorted - 排序後的特徵值對角矩陣（由大到小排序）
            
            [V, D] = eig(A);
            eigVals = diag(D);
            [sortedEigVals, idx] = sort(eigVals, 'descend');  % 由大到小排序
            D_sorted = diag(sortedEigVals);
            V_sorted = V(:, idx);
        end

        function A_clean = cleanMatrix(obj, A, threshold)
            % 清理矩陣中的極小值
            if nargin < 3
                threshold = obj.config.CleanThreshold;
            end
            A_clean = A;
            A_clean(abs(A) < threshold) = 0;
        end

        function validateStrainTensor(obj, Eg)
            % 驗證應變張量的物理合理性
%{
            eigenvals = eig(Eg);

            if any(eigenvals <= 0)
                warning('應變張量特徵值包含非正值，可能不符合物理意義');
            end
%}            
            % 檢查對稱性
            if max(max(abs(Eg - Eg'))) > obj.config.CleanThreshold * 100
                warning('應變張量不對稱程度超過容許範圍');
            end
%{            
            % 檢查條件數
            if cond(Eg) > 1e12
                warning('應變張量條件數過大，可能存在數值問題');
            end
%}
        end

        function [U, R, omega] = polarDecomposition(obj, F_increment)
            % 極分解：F = R * U
            % 
            % 輸入：
            %   F_increment - 增量變形梯度張量
            % 輸出：
            %   U - 對稱拉伸張量
            %   R - 旋轉張量  
            %   omega - 旋轉率張量
            
            % 計算右Cauchy-Green張量
            C = F_increment' * F_increment;
            
            % 特徵值分解
            [V, D] = obj.eigSorted(C);
            
            % 確保特徵值為正
            eigenvals = diag(D);
            if any(eigenvals <= 0)
                warning('Cauchy-Green張量特徵值非正，極分解可能不準確');
                eigenvals = max(eigenvals, 1e-10);
                D = diag(eigenvals);
            end
            
            % 計算對稱拉伸張量
            U = V * sqrt(D) * V';
            
            % 計算旋轉張量
            try
                R = F_increment / U;
            catch
                warning('旋轉張量計算遇到數值問題，使用偽逆');
                R = F_increment * pinv(U);
            end
            
            % 計算旋轉率張量
            omega = (R - R') / 2;
            
            if obj.config.Verbose
                fprintf('極分解完成，最大拉伸值：%.4f\n', sqrt(max(eigenvals)));
            end
        end

        function computeIncrementalStrain(obj, sampleA_idx, sampleB_idx, resultName)
            % 計算兩個樣本間的增量應變
            %
            % 輸入：
            %   sampleA_idx - 初始樣本索引
            %   sampleB_idx - 最終樣本索引
            %   resultName - 結果變數名稱
            
            if isempty(obj.Eg)
                error('請先執行 computeFiniteStrainTensors()');
            end
            
            if sampleA_idx > size(obj.Eg,3) || sampleB_idx > size(obj.Eg,3)
                error('樣本索引超出範圍');
            end
            
            Eg_initial = obj.Eg(:,:,sampleA_idx);
            Eg_final = obj.Eg(:,:,sampleB_idx);
            
            % 計算增量變形梯度
            try
                F_increment = Eg_final / Eg_initial;
            catch
                warning('使用偽逆計算增量變形梯度');
                F_increment = Eg_final * pinv(Eg_initial);
            end
            
            % 對稱化處理
            F_increment_sym = (F_increment + F_increment') / 2;
            
            % 極分解
            [U, R, omega] = obj.polarDecomposition(F_increment);
            
            % 特徵值分析
            [V_eigen, D_eigen] = obj.eigSorted(U);
            
            % 儲存結果
            result = struct();
            result.F_increment = F_increment;
            result.F_increment_sym = F_increment_sym;
            result.U = U;
            result.R = R;
            result.omega = omega;
            result.eigenvalues = diag(D_eigen);
            result.eigenvectors = V_eigen;
            result.sampleA_idx = sampleA_idx;
            result.sampleB_idx = sampleB_idx;
            
            obj.incrementalResults.(resultName) = result;
            
            % 輸出到工作區
            assignin('base', sprintf('%s_F', resultName), F_increment);
            assignin('base', sprintf('%s_U', resultName), U);
            assignin('base', sprintf('%s_R', resultName), R);
            assignin('base', sprintf('%s_omega', resultName), omega);
            assignin('base', sprintf('%s_eigenvals', resultName), diag(D_eigen));
            assignin('base', sprintf('%s_eigenvecs', resultName), V_eigen);
            
            if obj.config.Verbose
                fprintf('\n增量應變分析完成：%s\n', resultName);
                fprintf('樣本 %d → 樣本 %d\n', sampleA_idx, sampleB_idx);
                fprintf('主拉伸值：[%.4f, %.4f, %.4f]\n', result.eigenvalues);
            end
        end

        function combineAndAnalyzeEg(obj, varA, varB, config)
            % 合併兩個 Eg 變數並分析
            if nargin < 4
                config = struct('CleanSmallValues', true, 'Threshold', obj.config.CleanThreshold);
            end

            % 嘗試取得變數
            try
                Eg_A = evalin('base', varA);
                Eg_B = evalin('base', varB);
            catch
                fprintf('找不到其中一個變數：%s 或 %s\n', varA, varB);
                return;
            end
            
            V_A = obj.tryGetVFromEgName(varA);
            V_B = obj.tryGetVFromEgName(varB);
            
            % 執行 Eg 相乘
            Eg_combined = Eg_B * Eg_A;
            Eg_sym = (Eg_combined + Eg_combined') / 2;

            % 可選：清除極小值
            if config.CleanSmallValues
                Eg_combined = obj.cleanMatrix(Eg_combined, config.Threshold);
                Eg_sym = obj.cleanMatrix(Eg_sym, config.Threshold);
            end

            % 特徵值分析（並排序）
            [V_sorted, D_sorted] = obj.eigSorted(Eg_sym);

            % 找下一個命名編號
            existingVars = evalin('base', 'who');
            matched = regexp(existingVars, '^Eg_combined_\d+$', 'match');
            nextID = sum(~cellfun('isempty', matched)) + 1;
            varName = @(base) sprintf('%s_%d', base, nextID);

            % 顯示處理流程
            fprintf('\nEg 組合分析（第 %d 次）: Eg_B * Eg_A\n', nextID);
            fprintf('Eg_A (%s):\n', varA); disp(Eg_A);
            if ~isnan(V_A), fprintf('V_A (%s):\n', varA); disp(V_A); end
            fprintf('Eg_B (%s):\n', varB); disp(Eg_B);
            if ~isnan(V_B), fprintf('V_B (%s):\n', varB); disp(V_B); end
            fprintf('Eg_combined = Eg_B * Eg_A:\n'); disp(Eg_combined);
            fprintf('對稱化 Eg_sym:\n'); disp(Eg_sym);
            fprintf('特徵值 (排序):\n'); disp(D_sorted);
            fprintf('特徵向量 (排序):\n'); disp(V_sorted);

            % 儲存結果到 base workspace
            assignin('base', varName('Eg_combined'), Eg_combined);
            assignin('base', varName('Eg_sym'), Eg_sym);
            assignin('base', varName('Eg_combined_V'), V_sorted);
            assignin('base', varName('Eg_combined_D'), D_sorted);
            assignin('base', varName('Eg_combined_sources'), {varA, varB});

            fprintf('\n儲存變數：%s, %s, %s, %s\n來源記錄：%s\n', ...
                varName('Eg_combined'), varName('Eg_sym'), ...
                varName('Eg_combined_V'), varName('Eg_combined_D'), ...
                varName('Eg_combined_sources'));
        end

        function averageAndAnalyzeEg(obj, varList, config)
            % 對多個 Eg 變數進行平均並分析
            if nargin < 3
                config = struct('CleanSmallValues', true, 'Threshold', obj.config.CleanThreshold);
            end

            % 自動偵測變數
            if nargin < 2 || isempty(varList)
                allVars = evalin('base', 'who');
                matched = regexp(allVars, '^(Eg_|Eg_raw_)\d+$', 'match');
                varList = [matched{:}];
                if isempty(varList)
                    disp('找不到任何 Eg_* 或 Eg_raw_* 變數可用於平均');
                    return;
                end
                fprintf('自動偵測到以下 Eg 變數將進行平均：\n');
                disp(varList');
            elseif ~iscell(varList)
                disp('請以 cell array 格式輸入，如 {''Eg_1'', ''Eg_raw_2''}');
                return;
            end

            % 讀取 Eg 矩陣
            Eg_matrices = [];
            for i = 1:length(varList)
                try
                    Eg_i = evalin('base', varList{i});
                    Eg_matrices(:, :, i) = Eg_i;
                catch
                    fprintf('找不到變數：%s\n', varList{i});
                    return;
                end
            end

            % 計算平均與對稱化
            Eg_avg = mean(Eg_matrices, 3);
            Eg_sym = (Eg_avg + Eg_avg') / 2;

            if config.CleanSmallValues
                Eg_avg = obj.cleanMatrix(Eg_avg, config.Threshold);
                Eg_sym = obj.cleanMatrix(Eg_sym, config.Threshold);
            end

            % 特徵分析
            [V_sorted, D_sorted] = obj.eigSorted(Eg_sym);

            % 命名與編號
            existingVars = evalin('base', 'who');
            matched = regexp(existingVars, '^Eg_avg_\d+$', 'match');
            nextID = sum(~cellfun('isempty', matched)) + 1;
            varName = @(base) sprintf('%s_%d', base, nextID);

            % 顯示處理過程
            fprintf('\n平均 Eg (%d 個):\n', length(varList)); disp(Eg_avg);
            fprintf('對稱化:\n'); disp(Eg_sym);
            fprintf('特徵值:\n'); disp(D_sorted);
            fprintf('特徵向量:\n'); disp(V_sorted);

            % 儲存至 base workspace
            assignin('base', varName('Eg_avg'), Eg_avg);
            assignin('base', varName('Eg_avg_sym'), Eg_sym);
            assignin('base', varName('Eg_avg_V'), V_sorted);
            assignin('base', varName('Eg_avg_D'), D_sorted);
            assignin('base', varName('Eg_avg_sources'), varList);

            fprintf('\n平均後結果儲存為：\n%s\n%s\n%s\n%s\n來源清單：%s\n', ...
                varName('Eg_avg'), varName('Eg_avg_sym'), ...
                varName('Eg_avg_V'), varName('Eg_avg_D'), varName('Eg_avg_sources'));
        end

        function V = tryGetVFromEgName(~, egName)
            % 嘗試從 egName 推出對應的特徵向量名稱
            V = NaN;
            try
                id = regexp(egName, '\d+$', 'match');
                if ~isempty(id)
                    Vname = ['V_' id{1}];
                    V = evalin('base', Vname);
                end
            catch
                % 無對應值則傳回 NaN
            end
        end

        function exportResults(obj, outputPath, format)
            % 匯出分析結果
            if nargin < 3
                format = 'excel';
            end
            
            switch lower(format)
                case 'excel'
                    obj.exportToExcel(outputPath);
                case 'csv'
                    obj.exportToCSV(outputPath);
                case 'mat'
                    obj.exportToMat(outputPath);
                otherwise
                    error('不支援的輸出格式：%s', format);
            end
        end

        function exportToExcel(obj, filename)
            % 匯出到 Excel 檔案
            if isempty(obj.Eg)
                error('請先執行應變張量計算');
            end
            
            % 建立結果表格
            results = table();
            n = size(obj.Eg, 3);
            
            for i = 1:n
                [V, D] = obj.eigSorted(obj.Eg(:,:,i));
                eigenvals = diag(D);
                
                results.Sample(i) = i;
                results.MaxStrain(i) = eigenvals(1);
                results.IntStrain(i) = eigenvals(2);
                results.MinStrain(i) = eigenvals(3);
                results.StrainRatio(i) = eigenvals(1) / eigenvals(3);
                
                % 特徵向量方向（以度為單位）
                results.MaxStrainTrend(i) = atan2d(V(2,1), V(1,1));
                results.MaxStrainPlunge(i) = asind(V(3,1));
            end
            
            writetable(results, filename, 'Sheet', '應變分析結果');
            fprintf('結果已匯出至：%s\n', filename);
        end

        function exportToCSV(obj, filename)
            % 匯出到 CSV 檔案
            if isempty(obj.Eg)
                error('請先執行應變張量計算');
            end
            
            % 簡化的結果表格
            results = table();
            n = size(obj.Eg, 3);
            
            for i = 1:n
                [~, D] = obj.eigSorted(obj.Eg(:,:,i));
                eigenvals = diag(D);
                
                results.Sample(i) = i;
                results.MaxStrain(i) = eigenvals(1);
                results.IntStrain(i) = eigenvals(2);
                results.MinStrain(i) = eigenvals(3);
            end
            
            writetable(results, filename);
            fprintf('結果已匯出至：%s\n', filename);
        end

        function exportToMat(obj, filename)
            % 匯出到 MAT 檔案
            analyzer_data = struct();
            analyzer_data.Eg = obj.Eg;
            analyzer_data.ErList = obj.ErList;
            analyzer_data.VList = obj.VList;
            analyzer_data.config = obj.config;
            analyzer_data.filename = obj.filename;
            analyzer_data.incrementalResults = obj.incrementalResults;
            
            save(filename, 'analyzer_data');
            fprintf('完整分析結果已儲存至：%s\n', filename);
        end

        function plotStrainEvolution(obj, sampleRange)
            % 繪製應變演化圖
            if isempty(obj.Eg)
                error('請先執行應變張量計算');
            end
            
            if nargin < 2
                sampleRange = 1:min(10, size(obj.Eg,3));
            end
            
            figure;
            hold on;
            colors = lines(length(sampleRange));
            
            for idx = 1:length(sampleRange)
                i = sampleRange(idx);
                [~, D] = obj.eigSorted(obj.Eg(:,:,i));
                eigenvals = diag(D);
                plot3(eigenvals(1), eigenvals(2), eigenvals(3), 'o-', ...
                      'Color', colors(idx,:), 'MarkerSize', 8, 'LineWidth', 2);
                text(eigenvals(1), eigenvals(2), eigenvals(3), sprintf('  %d', i));
            end
            
            xlabel('最大主應變');
            ylabel('中間主應變');
            zlabel('最小主應變');
            title('應變橢球演化軌跡');
            grid on;
            view(3);
            legend(arrayfun(@(x) sprintf('樣本 %d', x), sampleRange, 'UniformOutput', false));
        end

        function printSummary(obj)
            % 列印分析摘要
            if isempty(obj.Eg)
                fprintf('尚未執行應變張量計算\n');
                return;
            end
            
            fprintf('\n=== AMS 應變分析摘要 ===\n');
            fprintf('檔案：%s\n', obj.filename);
            fprintf('樣本數量：%d\n', size(obj.Eg, 3));
            fprintf('Slate 係數 A：%.3f, B：%.3f\n', obj.config.SlateCoeffA, obj.config.SlateCoeffB);
            
            % 計算統計資訊
            all_eigenvals = [];
            for i = 1:size(obj.Eg, 3)
                [~, D] = obj.eigSorted(obj.Eg(:,:,i));
                all_eigenvals = [all_eigenvals; diag(D)'];
            end
            
            fprintf('應變值範圍：\n');
            fprintf('  最大主應變：%.4f - %.4f\n', min(all_eigenvals(:,1)), max(all_eigenvals(:,1)));
            fprintf('  中間主應變：%.4f - %.4f\n', min(all_eigenvals(:,2)), max(all_eigenvals(:,2)));
            fprintf('  最小主應變：%.4f - %.4f\n', min(all_eigenvals(:,3)), max(all_eigenvals(:,3)));
            
            if ~isempty(fieldnames(obj.incrementalResults))
                fprintf('增量應變分析：%d 個結果\n', length(fieldnames(obj.incrementalResults)));
            end
            
            fprintf('========================\n');
        end
    end
end