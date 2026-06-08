classdef AMSFormulas
    % AMSFormulas - AMS 應變分析公式庫
    % 
    % 此類別包含所有磁化率轉應變的核心公式計算
    % 不包含資料處理、UI 互動等功能
    
    methods (Static)
        
        %% ==================== 核心應變計算公式 ====================
        
        function Er = calculateEr(K1, K2, K3, slateA, slateB)
            % 計算磁化率應變張量 (Slate 模型)
            %
            % 輸入：
            %   K1, K2, K3  - 主磁化率值 (K1 ≥ K2 ≥ K3)
            %   slateA, slateB - Slate 係數 (預設: 6.897, 0.007)
            %
            % 輸出：
            %   Er - 3×3 主軸座標系應變張量
            %
            % 公式：
            %   K0 = (K1·K2·K3)^(1/3)
            %   ei = exp[a·(Ki/K0 - 1) - b] - 1
            %   ω = (1+e1)(1+e2)(1+e3)
            %   Er = ω² · diag[(1+e1)^(-2), (1+e2)^(-2), (1+e3)^(-2)]
            
            % 參數檢查
            if K1 <= 0 || K2 <= 0 || K3 <= 0
                error('磁化率值必須為正數');
            end
            
            % 計算平均磁化率
            K0 = (K1 * K2 * K3)^(1/3);
            
            % 計算應變分量
            e1 = exp(slateA * ((K1/K0) - 1) - slateB) - 1;
            e2 = exp(slateA * ((K2/K0) - 1) - slateB) - 1;
            e3 = exp(slateA * ((K3/K0) - 1) - slateB) - 1;
            
            % 計算體積因子
            omega = (1 + e1) * (1 + e2) * (1 + e3);
            
            % 組裝應變張量
            Er = omega^2 * diag([(1+e1)^-2, (1+e2)^-2, (1+e3)^-2]);
        end
        
        
        function V = calculateV(trend1, plunge1, trend2, plunge2, trend3, plunge3)
            % 計算特徵向量矩陣 (方向餘弦矩陣)
            % Eigen Vector
            %
            % [修改說明]：使用 cosd/sind 直接處理角度 (Degrees)
            %
            % 輸入：
            %   trend1, plunge1  - K1 軸的方位角和傾角 (度 Degree)
            %   trend2, plunge2  - K2 軸的方位角和傾角 (度 Degree)
            %   trend3, plunge3  - K3 軸的方位角和傾角 (度 Degree)
            %
            % 輸出：
            %   V - 3×3 特徵向量矩陣
            %
            % 公式 (地理座標系：X=北, Y=東, Z=下):
            %   V = [cosd(i1)·cosd(d1), cosd(i2)·cosd(d2), cosd(i3)·cosd(d3);
            %        cosd(i1)·sind(d1), cosd(i2)·sind(d2), cosd(i3)·sind(d3);
            %        sind(i1),          sind(i2),          sind(i3)]
            
            V = [cosd(plunge1)*cosd(trend1), cosd(plunge2)*cosd(trend2), cosd(plunge3)*cosd(trend3);
                 cosd(plunge1)*sind(trend1), cosd(plunge2)*sind(trend2), cosd(plunge3)*sind(trend3);
                 sind(plunge1),              sind(plunge2),              sind(plunge3)];
        end
        
        
        function Eg = calculateEg(Er, V)
            % 座標轉換：主軸座標系 → 地理座標系
            %
            % 輸入：
            %   Er - 主軸座標系應變張量
            %   V  - 特徵向量矩陣
            %
            % 輸出：
            %   Eg - 地理座標系應變張量
            %
            % 公式：
            %   Eg = V · Er · V'
            
            Eg = V * Er * V';
        end
        
        
        %% ==================== 增量應變計算公式 ====================
        
        function F = calculateF(Eg_initial, Eg_final)
            % 計算增量變形梯度張量
            %
            % 輸入：
            %   Eg_initial - 初始狀態應變張量
            %   Eg_final   - 最終狀態應變張量
            %
            % 輸出：
            %   F - 應變增量
            %
            % 公式：
            %   F = Eg_final · Eg_initial^(-1)
            
            F = Eg_final * Eg_initial^-1;
        end
        
        
        function [U, R, omega] = polarDecomposition(F)
            % 極分解：F = R · U
            %
            % 輸入：
            %   F - 變形梯度張量
            %
            % 輸出：
            %   U     - 右拉伸張量 (對稱正定)
            %   R     - 旋轉張量 (正交)
            %   omega - 旋轉率張量 (反對稱)
            %
            % 公式：
            %   C = F' · F  (右 Cauchy-Green 張量)
            %   U = sqrt(C) (特徵值分解)
            %   R = F · U^(-1)
            %   ω = (R - R') / 2
            
            % 計算右 Cauchy-Green 張量
            C = F' * F;
            
            % 特徵值分解
            [V, D] = eig(C);
            eigenvals = diag(D);
            
            % 確保特徵值為正
            if any(eigenvals <= 0)
                warning('Cauchy-Green 張量特徵值非正，進行修正');
                eigenvals = max(eigenvals, 1e-10);
                D = diag(eigenvals);
            end
            
            % 計算對稱拉伸張量
            U = V * sqrt(D) * V';
            
            % 計算旋轉張量
            try
                R = F / U;
            catch
                warning('使用偽逆計算旋轉張量');
                R = F * pinv(U);
            end
            
            % 計算旋轉率張量 (反對稱部分)
            omega = (R - R') / 2;
        end
        
        
        %% ==================== 張量運算工具 ====================
        
        function [V_sorted, D_sorted] = sortedEigenDecomposition(A)
            % 特徵值分解並依特徵值降序排列
            %
            % 輸入：
            %   A - 對稱矩陣
            %
            % 輸出：
            %   V_sorted - 排序後的特徵向量矩陣
            %   D_sorted - 排序後的特徵值對角矩陣 (降序)
            
            [V, D] = eig(A);
            eigenvals = diag(D);
            [sortedEigenvals, idx] = sort(eigenvals, 'descend');
            D_sorted = diag(sortedEigenvals);
            V_sorted = V(:, idx);
        end
        
        
        function A_sym = symmetrize(A)
            % 對稱化矩陣
            %
            % 輸入：
            %   A - 任意矩陣
            %
            % 輸出：
            %   A_sym - 對稱化後的矩陣
            %
            % 公式：
            %   A_sym = (A + A') / 2
            
            A_sym = (A + A') / 2;
        end
        
        
        %% ==================== 幾何計算工具 ====================
        
        function [trend, plunge] = vectorToTrendPlunge(v)
            % 1. 檢查是否指向「天上」(假設 Z 正向是向下，若 v(3) 為負則是指向天)
            %    或者簡單來說，只要 Plunge 算出來是負的，就反轉向量
            
            % 先計算傾角
            plunge = asind(v(3));
            
            % 計算方位角
            trend = atan2d(v(2), v(1));
            if trend < 0
                trend = trend + 360;
            end
            
            % [新增] 修正邏輯：強制 Plunge 為正 (Lower Hemisphere)
            if plunge < 0
                plunge = -plunge;       % 傾角變正
                trend = trend + 180;    % 方位角反轉 180 度
                
                % 確保 Trend 在 0-360 範圍內
                if trend >= 360
                    trend = trend - 360;
                end
            end
        end
        
        %{
        function v = trendPlungeToVector(trend, plunge)
            % 將方位角和傾角轉換為單位向量
            %
            % [修改說明]：移除 deg2rad，改用 cosd/sind 直接處理角度
            %
            % 輸入：
            %   trend  - 方位角 (度)
            %   plunge - 傾角 (度)
            %
            % 輸出：
            %   v - 3×1 單位向量 [x; y; z]
            
            % 直接使用 cosd 和 sind
            v = [cosd(plunge) * cosd(trend);
                 cosd(plunge) * sind(trend);
                 sind(plunge)];
        end
        %}
        
        %% ==================== 驗證與檢查 ====================
        
        function isValid = validateMagneticSusceptibility(K1, K2, K3)
            % 驗證磁化率資料的合理性
            %
            % 檢查項目：
            %   1. 數值為正
            %   2. 符合 K1 ≥ K2 ≥ K3 順序
            
            isValid = true;
            
            if K1 <= 0 || K2 <= 0 || K3 <= 0
                warning('磁化率值必須為正數');
                isValid = false;
            end
            
            if K1 < K2 || K2 < K3
                warning('磁化率不符合 K1 ≥ K2 ≥ K3 順序');
                isValid = false;
            end
        end
        
        
        function isSymmetric = checkSymmetry(A, tolerance)
            % 檢查矩陣對稱性
            %
            % 輸入：
            %   A         - 待檢查矩陣
            %   tolerance - 容差 (預設 1e-10)
            %
            % 輸出：
            %   isSymmetric - 是否對稱
            
            if nargin < 2
                tolerance = 1e-10;
            end
            
            isSymmetric = max(max(abs(A - A'))) < tolerance;
        end
        
        
        function isOrthogonal = checkOrthogonality(V, tolerance)
            % 檢查矩陣正交性
            %
            % 輸入：
            %   V         - 待檢查矩陣
            %   tolerance - 容差 (預設 1e-10)
            %
            % 輸出：
            %   isOrthogonal - 是否正交
            
            if nargin < 2
                tolerance = 1e-10;
            end
            
            I = eye(size(V));
            isOrthogonal = max(max(abs(V'*V - I))) < tolerance;
        end
        
        %% ==================== 模型驗證 ====================
        function Kg = calculateKg(K1, K2, K3, V)
            % 計算地理座標系中的磁化率張量 (Susceptibility Tensor)
            %
            % 這是在將磁化率轉換為應變(Er) *之前* 的原始張量
            %
            % 輸入：
            %   K1, K2, K3  - 主磁化率值 (K1 ≥ K2 ≥ K3)
            %   V           - 特徵向量矩陣 (來自 calculateV)
            %
            % 輸出：
            %   Kg - 3×3 地理座標系磁感率張量
            %
            % 公式：
            %   K_diag = diag([K1, K2, K3])  (主軸座標系張量)
            %   Kg = V · K_diag · V'          (轉換到地理座標系)
            
            % 1. 組裝主軸座標系的磁感率張量 (對角矩陣)
            K_diag = diag([K1, K2, K3]);
            
            % 2. 執行座標轉換 (主軸 -> 地理)
            % 注意：標準的張量座標轉換公式為 T_geo = V * T_diag * V'
            % 其中 V 是從主軸到地理的轉換矩陣 (如 calculateV 所計算)
            % V' 是 V 的轉置 (因為 V 是正交的, V' = inv(V))
            Kg = V' * K_diag * V;
        end
        
        function FK = calculateFK(Kg_initial, Kg_final)
            % 計算磁感率張量增量
            %
            % 輸入：
            %   Kg_initial - 初始狀態磁感率張量
            %   Kg_final   - 最終狀態磁感率張量
            %
            % 輸出：
            %   FK - 磁感率增量
            %
            % 公式：
            %   FK = Kg_final · Kg_initial^(-1)
            
            %try
                FK = Kg_final * Kg_initial^-1;
            %catch
                %warning('使用偽逆計算變形梯度');
                %FK = Kg_final * pinv(Kg_initial);
            %end
        end

        %% ==================== 統計與平均 (Statistics & Averaging) ====================
        
        function K_norm = standardizeTensor(K)
            % 標準化磁感率張量 (Standardize/Normalize Tensor)
            % 目的：消除標本間磁感率強度的差異，僅保留形狀與方向資訊
            % 依據：Jelinek (1978)
            %
            % 輸入：
            %   K - 3x3 磁感率張量 (Geographic or Specimen coordinates)
            %
            % 輸出：
            %   K_norm - 標準化後的張量 (每個元素除以平均磁感率)
            %
            % 公式：
            %   k_mean = (K1 + K2 + K3) / 3 = trace(K) / 3
            %   K_norm = K / k_mean
            
            k_mean = trace(K) / 3;
            
            % 避免除以零的錯誤
            if k_mean == 0
                warning('平均磁感率為 0，無法標準化，回傳原始張量。');
                K_norm = K;
            else
                K_norm = K / k_mean;
            end
        end
        
        function M = calculateMeanTensor(tensorList)
            % 計算平均張量 (Mean Tensor)
            % 目的：計算一組張量的分量平均值
            %
            % 輸入：
            %   tensorList - 3x3xN 矩陣
            %                (包含 N 個樣本的磁感率張量，需先堆疊成 3D 陣列)
            %
            % 輸出：
            %   M - 3x3 平均張量
            %
            % 公式：
            %   M_ij = (1/N) * sum(K_ij)_n  (對每個分量做算術平均)
            
            % 檢查輸入維度
            if ismatrix(tensorList)
                if all(size(tensorList) == [3 3])
                    % 只有傳入一個張量，平均就是自己
                    M = tensorList;
                    return;
                else
                    error('輸入必須是 3x3xN 的矩陣堆疊 (3D Array)');
                end
            end
            
            % 沿著第三維度 (樣本數) 取平均
            M = mean(tensorList, 3);
        end
        function results = performEigenAnalysis(Tensor)
            % 執行完整的特徵值分析 (Eigen-analysis)
            % 目的：從張量中提取主成分 (K1, K2, K3) 及其地理方向 (D, I)
            %
            % 輸入：
            %   Tensor - 3x3 對稱張量 (如 Kg, 平均張量, 或應變張量)
            %
            % 輸出：
            %   results - 包含分析結果的結構體 (Struct)，欄位如下：
            %       .vals   : [K1, K2, K3] 特徵值陣列 (由大到小)
            %       .V      : 3x3 特徵向量矩陣 (對應 K1, K2, K3)
            %       .K1_dir : [Trend, Plunge] K1 的方向
            %       .K2_dir : [Trend, Plunge] K2 的方向
            %       .K3_dir : [Trend, Plunge] K3 的方向
            
            % 1. 數學分解：計算特徵值與特徵向量 (並排序)
            [V_sorted, D_sorted] = AMSFormulas.sortedEigenDecomposition(Tensor);
            
            % 提取特徵值 (對角線元素)
            eigenvals = diag(D_sorted);
            
            % 2. 地質轉換：將特徵向量轉換為 Trend/Plunge
            % V_sorted 的每一行 (Column) 代表一個特徵向量
            [t1, p1] = AMSFormulas.vectorToTrendPlunge(V_sorted(:, 1)); %對應 K1
            [t2, p2] = AMSFormulas.vectorToTrendPlunge(V_sorted(:, 2)); %對應 K2
            [t3, p3] = AMSFormulas.vectorToTrendPlunge(V_sorted(:, 3)); %對應 K3
            
            % 3. 打包結果
            results.vals = eigenvals;      % [K1; K2; K3]
            results.V    = V_sorted;       % 排序後的矩陣
            
            results.K1_val = eigenvals(1);
            results.K1_dir = [t1, p1];     % [Trend, Plunge]
            
            results.K2_val = eigenvals(2);
            results.K2_dir = [t2, p2];
            
            results.K3_val = eigenvals(3);
            results.K3_dir = [t3, p3];
            
            % 額外計算一些形狀參數 (Jelinek參數) 方便參考
            % Pj (Corrected Anisotropy Degree) 等參數可視需求加入
            % 這裡先提供簡單的 L (Lineation), F (Foliation)
            if eigenvals(3) > 0
                results.L = eigenvals(1) / eigenvals(2); % 線理
                results.F = eigenvals(2) / eigenvals(3); % 葉理
                results.P = eigenvals(1) / eigenvals(3); % 異向性度
            else
                results.L = NaN; results.F = NaN; results.P = NaN;
            end
        end
        %% ==================== 1. Jelinek (1978) 統計分析 ====================
        function stats = calculateJelinekStats(tensorList)
            % 計算 Jelinek (1978) 信心橢圓角度 (95%)
            [~, ~, N] = size(tensorList);
            
            if N < 3
                stats = struct('E12',NaN, 'E23',NaN, 'E13',NaN, 'K1_err',NaN, 'K2_err',NaN, 'K3_err',NaN);
                return;
            end
            
            % 1. 平均張量與旋轉
            M = AMSFormulas.calculateMeanTensor(tensorList);
            [V_mean, D_mean] = AMSFormulas.sortedEigenDecomposition(M);
            lambda = diag(D_mean); 
            
            % 2. 旋轉至主軸座標系並提取分量
            k_rot_12 = zeros(N,1); k_rot_23 = zeros(N,1); k_rot_13 = zeros(N,1);
            k_rot_11 = zeros(N,1); k_rot_22 = zeros(N,1); k_rot_33 = zeros(N,1);
            
            for i = 1:N
                K_rot = V_mean' * tensorList(:,:,i) * V_mean;
                k_rot_12(i) = K_rot(1,2); k_rot_23(i) = K_rot(2,3); k_rot_13(i) = K_rot(1,3);
                k_rot_11(i) = K_rot(1,1); k_rot_22(i) = K_rot(2,2); k_rot_33(i) = K_rot(3,3);
            end
            
            % 3. 計算變異數與 F 因子 (95% 信心水準)
            var_12 = sum(k_rot_12.^2) / (N*(N-2));
            var_23 = sum(k_rot_23.^2) / (N*(N-2));
            var_13 = sum(k_rot_13.^2) / (N*(N-2));
            
            % 簡易 F 值查表近似
            if N <= 5, F=9.55; elseif N<=10, F=4.46; elseif N<=20, F=3.55; else, F=3.00; end
            factor = sqrt(2 * F);
            
            % 4. 計算信心角度
            stats.E12 = atand((factor * sqrt(var_12)) / abs(lambda(1)-lambda(2)));
            stats.E23 = atand((factor * sqrt(var_23)) / abs(lambda(2)-lambda(3)));
            stats.E13 = atand((factor * sqrt(var_13)) / abs(lambda(1)-lambda(3)));
            
            % 5. 數值標準誤
            stats.K1_err = std(k_rot_11)/sqrt(N);
            stats.K2_err = std(k_rot_22)/sqrt(N);
            stats.K3_err = std(k_rot_33)/sqrt(N);
        end

        %% ==================== 2. Bootstrap 統計分析 ====================
        function bootStats = calculateBootstrapStats(tensorList, numBootstraps)
            % 使用 Bootstrap 估算 K 值標準差與方向信心角
            if nargin < 2, numBootstraps = 1000; end
            [~, ~, N] = size(tensorList);
            
            if N < 3
                bootStats = struct('K1_std',NaN, 'V1_conf',NaN, 'K3_std',NaN, 'V3_conf',NaN);
                return;
            end
            
            % 基準方向 (用於校正向量反轉)
            Mean_Kg = AMSFormulas.calculateMeanTensor(tensorList);
            [V_ref, ~] = AMSFormulas.sortedEigenDecomposition(Mean_Kg);
            
            boot_K = zeros(numBootstraps, 3);
            boot_V1 = zeros(numBootstraps, 3);
            boot_V3 = zeros(numBootstraps, 3); % 通常最關心 K1 (線理) 和 K3 (極)
            
            for b = 1:numBootstraps
                indices = randi(N, N, 1);
                M_b = AMSFormulas.calculateMeanTensor(tensorList(:,:,indices));
                [V_b, D_b] = AMSFormulas.sortedEigenDecomposition(M_b);
                
                boot_K(b,:) = diag(D_b)';
                
                % 向量校正
                v1 = V_b(:,1); if dot(v1, V_ref(:,1)) < 0, v1 = -v1; end
                v3 = V_b(:,3); if dot(v3, V_ref(:,3)) < 0, v3 = -v3; end
                boot_V1(b,:) = v1';
                boot_V3(b,:) = v3';
            end
            
            % 統計結果
            bootStats.K1_std = std(boot_K(:,1));
            bootStats.K2_std = std(boot_K(:,2));
            bootStats.K3_std = std(boot_K(:,3));
            
            bootStats.V1_conf = AMSFormulas.calcConfidenceAngle(boot_V1, V_ref(:,1));
            bootStats.V3_conf = AMSFormulas.calcConfidenceAngle(boot_V3, V_ref(:,3));
        end
        
        function conf_angle = calcConfidenceAngle(boot_vectors, ref_vector)
            % 計算 95% 信心圓半徑
            num = size(boot_vectors, 1);
            angles = zeros(num, 1);
            for i = 1:num
                d = dot(boot_vectors(i,:)', ref_vector);
                angles(i) = acosd(max(min(d, 1), -1));
            end
            angles = sort(angles);
            conf_angle = angles(floor(0.95 * num));
        end
    end
end

