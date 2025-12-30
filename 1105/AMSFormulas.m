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
            %   Eg = V' · Er · V
            
            Eg = V' * Er * V;
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
            % 將單位向量轉換為方位角和傾角
            %
            % 輸入：
            %   v - 3×1 單位向量 [x; y; z] (北、東、下座標系)
            %
            % 輸出：
            %   trend  - 方位角 (度, 0-360°)
            %   plunge - 傾角 (度, -90° 到 +90°)
            
            trend = atan2d(v(2), v(1));
            if trend < 0
                trend = trend + 360;
            end
            plunge = asind(v(3));
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
            %   Kg - 3×3 地理座標系磁化率張量
            %
            % 公式：
            %   K_diag = diag([K1, K2, K3])  (主軸座標系張量)
            %   Kg = V · K_diag · V'          (轉換到地理座標系)
            
            % 1. 組裝主軸座標系的磁化率張量 (對角矩陣)
            K_diag = diag([K1, K2, K3]);
            
            % 2. 執行座標轉換 (主軸 -> 地理)
            % 注意：標準的張量座標轉換公式為 T_geo = V * T_diag * V'
            % 其中 V 是從主軸到地理的轉換矩陣 (如 calculateV 所計算)
            % V' 是 V 的轉置 (因為 V 是正交的, V' = inv(V))
            Kg = V * K_diag * V';
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
    end
end