classdef AMSFormulas
    % AMSFormulas - AMS 應變分析公式庫
    %
    % 座標系約定（全程統一）：
    %   地理座標系  X = 北, Y = 東, Z = 向下
    %   主軸座標系  軸1 = K1方向, 軸2 = K2方向, 軸3 = K3方向
    %
    % 張量轉換方向（V 的每一 column 是一個主軸在地理座標下的方向向量）：
    %   主軸 → 地理：T_geo = V  * T_princ * V'
    %   地理 → 主軸：T_princ = V' * T_geo   * V
    
    methods (Static)
        
        %% ==================== 核心應變計算公式 ====================
        
        function Er = calculateEr(K1, K2, K3, slateA, slateB)
            % 計算應變張量（主軸座標系，Slate 模型）
            %
            % 輸入：
            %   K1, K2, K3    - 主磁化率值 (K1 ≥ K2 ≥ K3)
            %   slateA, slateB - Slate 係數 (預設: 6.897, 0.007)
            %
            % 輸出：
            %   Er - 3×3 主軸座標系應變張量（對角矩陣）
            %
            % 公式：
            %   K0 = (K1·K2·K3)^(1/3)
            %   ei = exp[a·(Ki/K0 - 1) - b] - 1
            %   ω  = (1+e1)(1+e2)(1+e3)
            %   Er = ω² · diag[(1+e1)^-2, (1+e2)^-2, (1+e3)^-2]
            
            if K1 <= 0 || K2 <= 0 || K3 <= 0
                error('磁化率值必須為正數');
            end
            
            K0 = (K1 * K2 * K3)^(1/3);
            
            e1 = exp(slateA * ((K1/K0) - 1) - slateB) - 1;
            e2 = exp(slateA * ((K2/K0) - 1) - slateB) - 1;
            e3 = exp(slateA * ((K3/K0) - 1) - slateB) - 1;

            Er = diag([1+e1, (1+e2), (1+e3)]);
            %omega = (1 + e1) * (1 + e2) * (1 + e3);
             
            %Er = omega^2 * diag([(1+e1)^-2, (1+e2)^-2, (1+e3)^-2]);
        end
        
        
        function V = calculateV(trend1, plunge1, trend2, plunge2, trend3, plunge3)
            % 建構方向餘弦矩陣 V（主軸 → 地理座標轉換矩陣）
            %
            % 輸入：
            %   trend1, plunge1  - K1 軸的方位角和傾角（度）
            %   trend2, plunge2  - K2 軸的方位角和傾角（度）
            %   trend3, plunge3  - K3 軸的方位角和傾角（度）
            %
            % 輸出：
            %   V - 3×3 矩陣，每一 column 是對應主軸在地理座標下的單位向量
            %       V(:,1) = K1方向, V(:,2) = K2方向, V(:,3) = K3方向
            %
            % 地理座標系（X=北, Y=東, Z=向下）：
            %   Vx = cos(plunge)·cos(trend)
            %   Vy = cos(plunge)·sin(trend)
            %   Vz = sin(plunge)
            %
            % 張量轉換：
            %   Kg = V * diag([K1,K2,K3]) * V'   （主軸 → 地理）
            %   主軸張量 = V * Kg * V '            （地理 → 主軸）
            
            V = [cosd(plunge1)*cosd(trend1), cosd(plunge2)*cosd(trend2), cosd(plunge3)*cosd(trend3);
                 cosd(plunge1)*sind(trend1), cosd(plunge2)*sind(trend2), cosd(plunge3)*sind(trend3);
                 sind(plunge1),              sind(plunge2),              sind(plunge3)];
        end
        
        
        function Eg = calculateEg(Er, V_mean)
            % 座標轉換：主軸座標系應變張量 → 地理座標系應變張量
            %
            % 輸入：
            %   Er     - 主軸座標系應變張量（calculateEr 的輸出，對角矩陣）
            %   V_mean - 平均 Kg 特徵值分解得到的特徵向量矩陣
            %            每一 column 是對應主軸在地理座標下的單位向量
            %            （與 calculateV 的輸出格式相同）
            %
            % 輸出：
            %   Eg - 3×3 地理座標系應變張量
            %
            % 公式：
            %   Eg = V_mean * Er * V_mean'
            %
            % 注意：V_mean 來自 performEigenAnalysis 的 eigenRes.V，
            %       其每一 column 已對應排序後的主軸方向，
            %       與 calculateV 的 column 慣例一致。
            
            Eg = V_mean * Er * V_mean';
        end
        
        
        %% ==================== 增量應變計算公式 ====================
        
        function F = calculateF(Eg_initial, Eg_final)
            % 計算增量變形梯度張量
            %
            % 輸入：
            %   Eg_initial - 初始狀態地理座標系應變張量
            %   Eg_final   - 最終狀態地理座標系應變張量
            %
            % 輸出：
            %   F - 增量變形梯度 F = Eg_final · Eg_initial^(-1)
            
            F = Eg_final * Eg_initial^-1;
        end
        
        
        function [U, R, omega] = polarDecomposition(F)
            % 極分解：F = R · U
            %
            % 輸入：
            %   F - 變形梯度張量
            %
            % 輸出：
            %   U     - 右拉伸張量（對稱正定）
            %   R     - 旋轉張量（正交）
            %   omega - 旋轉率張量（反對稱部分）
            %
            % 公式：
            %   C = F' · F
            %   U = sqrt(C)
            %   R = F · U^(-1)
            %   ω = (R - R') / 2
            
            C = F' * F;
            [Vc, D] = eig(C);
            eigenvals = diag(D);
            
            if any(eigenvals <= 0)
                warning('Cauchy-Green 張量特徵值非正，進行修正');
                eigenvals = max(eigenvals, 1e-10);
                D = diag(eigenvals);
            end
            
            U = Vc * sqrt(D) * Vc';
            
            try
                R = F / U;
            catch
                warning('使用偽逆計算旋轉張量');
                R = F * pinv(U);
            end
            
            omega = (R - R') / 2;
        end
        
        
        %% ==================== 張量運算工具 ====================
        
        function [V_sorted, D_sorted] = sortedEigenDecomposition(A)
            % 特徵值分解並依特徵值降序排列
            %
            % 輸入：A - 對稱矩陣
            %
            % 輸出：
            %   V_sorted - 排序後的特徵向量矩陣（每 column 為一特徵向量）
            %   D_sorted - 排序後的特徵值對角矩陣（降序）
            %
            % 注意：輸出 V_sorted 的 column 慣例與 calculateV 相同，
            %       可直接用於 calculateEg 的 V_mean 參數。
            
            [V, D] = eig(A);
            eigenvals = diag(D);
            [sortedEigenvals, idx] = sort(eigenvals, 'descend');
            D_sorted = diag(sortedEigenvals);
            V_sorted = V(:, idx);
        end
        
        
        function A_sym = symmetrize(A)
            % 對稱化矩陣：A_sym = (A + A') / 2
            A_sym = (A + A') / 2;
        end
        
        
        %% ==================== 幾何計算工具 ====================
        
        function [trend, plunge] = vectorToTrendPlunge(v)
            % 將方向向量轉換為 Trend/Plunge（下半球投影）
            
            plunge = asind(v(3));
            trend  = atan2d(v(2), v(1));
            if trend < 0
                trend = trend + 360;
            end
            
            % 強制 Plunge 為正（下半球）
            if plunge < 0
                plunge = -plunge;
                trend  = trend + 180;
                if trend >= 360
                    trend = trend - 360;
                end
            end
        end
        
        
        %% ==================== 驗證與檢查 ====================
        
        function isValid = validateMagneticSusceptibility(K1, K2, K3)
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
            if nargin < 2, tolerance = 1e-10; end
            isSymmetric = max(max(abs(A - A'))) < tolerance;
        end
        
        function isOrthogonal = checkOrthogonality(V, tolerance)
            if nargin < 2, tolerance = 1e-10; end
            I = eye(size(V));
            isOrthogonal = max(max(abs(V'*V - I))) < tolerance;
        end
        
        
        %% ==================== 統計與平均 ====================
        
        function K_norm = standardizeTensor(K)
            % 標準化磁感率張量（Jelinek 1978）
            % K_norm = K / (trace(K)/3)
            k_mean = trace(K) / 3;
            if k_mean == 0
                warning('平均磁感率為 0，無法標準化，回傳原始張量。');
                K_norm = K;
            else
                K_norm = K / k_mean;
            end
        end
        
        
        function M = calculateMeanTensor(tensorList)
            % 計算平均張量（分量算術平均）
            %
            % 輸入：tensorList - 3×3×N 矩陣堆疊
            % 輸出：M         - 3×3 平均張量
            
            if ismatrix(tensorList)
                if all(size(tensorList) == [3 3])
                    M = tensorList;
                    return;
                else
                    error('輸入必須是 3×3×N 的矩陣堆疊');
                end
            end
            M = mean(tensorList, 3);
        end
        
        
        function results = performEigenAnalysis(Tensor)
            % 完整特徵值分析
            %
            % 輸入：Tensor - 3×3 對稱張量
            %
            % 輸出 results 欄位：
            %   .vals      [K1; K2; K3] 特徵值（降序）
            %   .V         3×3 特徵向量矩陣（每 column 為一主軸方向，降序）
            %              格式與 calculateV 一致，可直接傳入 calculateEg
            %   .K1_val / .K1_dir  [Trend, Plunge]
            %   .K2_val / .K2_dir
            %   .K3_val / .K3_dir
            %   .L  線理 = K1/K2
            %   .F  葉理 = K2/K3
            %   .P  異向性度 = K1/K3
            
            [V_sorted, D_sorted] = AMSFormulas.sortedEigenDecomposition(Tensor);
            eigenvals = diag(D_sorted);
            
            [t1, p1] = AMSFormulas.vectorToTrendPlunge(V_sorted(:,1));
            [t2, p2] = AMSFormulas.vectorToTrendPlunge(V_sorted(:,2));
            [t3, p3] = AMSFormulas.vectorToTrendPlunge(V_sorted(:,3));
            
            results.vals = eigenvals;
            results.V    = V_sorted;   % column = 主軸方向向量（與 calculateV 慣例相同）
            
            results.K1_val = eigenvals(1);
            results.K1_dir = [t1, p1];
            
            results.K2_val = eigenvals(2);
            results.K2_dir = [t2, p2];
            
            results.K3_val = eigenvals(3);
            results.K3_dir = [t3, p3];
            
            if eigenvals(3) > 0
                results.L = eigenvals(1) / eigenvals(2);
                results.F = eigenvals(2) / eigenvals(3);
                results.P = eigenvals(1) / eigenvals(3);
            else
                results.L = NaN; results.F = NaN; results.P = NaN;
            end
        end
        
        
        %% ==================== Jelinek (1978) 統計分析 ====================
        
        function stats = calculateJelinekStats(tensorList)
            % 計算 Jelinek (1978) 95% 信心橢圓角度
            
            [~, ~, N] = size(tensorList);
            
            if N < 3
                stats = struct('E12',NaN,'E23',NaN,'E13',NaN, ...
                               'K1_err',NaN,'K2_err',NaN,'K3_err',NaN);
                return;
            end
            
            M = AMSFormulas.calculateMeanTensor(tensorList);
            [V_mean, D_mean] = AMSFormulas.sortedEigenDecomposition(M);
            lambda = diag(D_mean);
            
            k_rot_12 = zeros(N,1); k_rot_23 = zeros(N,1); k_rot_13 = zeros(N,1);
            k_rot_11 = zeros(N,1); k_rot_22 = zeros(N,1); k_rot_33 = zeros(N,1);
            
            for i = 1:N
                % 地理 → 主軸：V' * K * V
                K_rot = V_mean' * tensorList(:,:,i) * V_mean;
                k_rot_12(i) = K_rot(1,2);
                k_rot_23(i) = K_rot(2,3);
                k_rot_13(i) = K_rot(1,3);
                k_rot_11(i) = K_rot(1,1);
                k_rot_22(i) = K_rot(2,2);
                k_rot_33(i) = K_rot(3,3);
            end
            
            var_12 = sum(k_rot_12.^2) / (N*(N-2));
            var_23 = sum(k_rot_23.^2) / (N*(N-2));
            var_13 = sum(k_rot_13.^2) / (N*(N-2));
            
            if N <= 5,    F = 9.55;
            elseif N<=10, F = 4.46;
            elseif N<=20, F = 3.55;
            else,         F = 3.00;
            end
            factor = sqrt(2 * F);
            
            stats.E12 = atand((factor * sqrt(var_12)) / abs(lambda(1)-lambda(2)));
            stats.E23 = atand((factor * sqrt(var_23)) / abs(lambda(2)-lambda(3)));
            stats.E13 = atand((factor * sqrt(var_13)) / abs(lambda(1)-lambda(3)));
            
            stats.K1_err = std(k_rot_11)/sqrt(N);
            stats.K2_err = std(k_rot_22)/sqrt(N);
            stats.K3_err = std(k_rot_33)/sqrt(N);
        end
        
        
        %% ==================== Bootstrap 統計分析 ====================
        
        function bootStats = calculateBootstrapStats(tensorList, numBootstraps)
            if nargin < 2, numBootstraps = 1000; end
            [~, ~, N] = size(tensorList);
            
            if N < 3
                bootStats = struct('K1_std',NaN,'K2_std',NaN,'K3_std',NaN, ...
                                   'V1_conf',NaN,'V3_conf',NaN);
                return;
            end
            
            Mean_Kg = AMSFormulas.calculateMeanTensor(tensorList);
            [V_ref, ~] = AMSFormulas.sortedEigenDecomposition(Mean_Kg);
            
            boot_K  = zeros(numBootstraps, 3);
            boot_V1 = zeros(numBootstraps, 3);
            boot_V3 = zeros(numBootstraps, 3);
            
            for b = 1:numBootstraps
                indices = randi(N, N, 1);
                M_b = AMSFormulas.calculateMeanTensor(tensorList(:,:,indices));
                [V_b, D_b] = AMSFormulas.sortedEigenDecomposition(M_b);
                
                boot_K(b,:) = diag(D_b)';
                
                v1 = V_b(:,1); if dot(v1, V_ref(:,1)) < 0, v1 = -v1; end
                v3 = V_b(:,3); if dot(v3, V_ref(:,3)) < 0, v3 = -v3; end
                boot_V1(b,:) = v1';
                boot_V3(b,:) = v3';
            end
            
            bootStats.K1_std  = std(boot_K(:,1));
            bootStats.K2_std  = std(boot_K(:,2));
            bootStats.K3_std  = std(boot_K(:,3));
            bootStats.V1_conf = AMSFormulas.calcConfidenceAngle(boot_V1, V_ref(:,1));
            bootStats.V3_conf = AMSFormulas.calcConfidenceAngle(boot_V3, V_ref(:,3));
        end
        
        
        function conf_angle = calcConfidenceAngle(boot_vectors, ref_vector)
            % 計算 95% 信心圓半徑（角度）
            num = size(boot_vectors, 1);
            angles = zeros(num, 1);
            for i = 1:num
                d = dot(boot_vectors(i,:)', ref_vector);
                angles(i) = acosd(max(min(d,1),-1));
            end
            angles = sort(angles);
            conf_angle = angles(floor(0.95 * num));
        end
        
    end
end