clc
function [Einc_intervals, U_intervals, R_intervals, U_eigenvalues_intervals, U_eigenvectors_intervals] = ...
    incremental_strain_interval_analysis(Eg, interval_count, samples_per_interval)


% INCREMENTAL_STRAIN_INTERVAL_ANALYSIS 計算指定區間的平均應變增量
%
% 輸入參數:
%   filename - Excel 檔案名稱，包含磁感率資料
%   interval_count - 要劃分的區間數量
%   samples_per_interval - 每個區間內的樣本數量
%
% 輸出參數:
%   Einc_intervals - 區間之間的應變增量
%   U_intervals - 區間之間的形變率張量
%   R_intervals - 區間之間的旋轉率張量
%   U_eigenvalues_intervals - 形變率張量的特徵值
%   U_eigenvectors_intervals - 形變率張量的特徵向量

    %% 讀取 Excel 檔案
    filename = 'NXHAMS.csv';
    data = readtable(filename, 'VariableNamingRule', 'preserve');
    
    %% 擷取並計算基本磁感率參數
    K1 = data.K1;   % 最大主軸磁感率
    K2 = data.K2;   % 中間主軸磁感率
    K3 = data.K3;   % 最小主軸磁感率
    Km = data.Km;   % 平均磁感率
    
    % 計算幾何平均磁感率 K0
    K0 = (K1 .* K2 .* K3).^(1/3);
    
    % 設定Wood et al. (1976)板岩參數
    a = 6.897;
    b = 0.007;
    
    %% 處理磁感率主軸方向數據
    % 獲取主軸方向
    dK1 = data.dK1geo;  % K1 方位角
    iK1 = data.iK1geo;  % K1 傾角
    dK2 = data.dK2geo;  % K2 方位角
    iK2 = data.iK2geo;  % K2 傾角
    dK3 = data.dK3geo;  % K3 方位角
    iK3 = data.iK3geo;  % K3 傾角
    
    % 轉換角度為弧度
    dK1 = deg2rad(dK1);
    iK1 = deg2rad(iK1);
    dK2 = deg2rad(dK2);
    iK2 = deg2rad(iK2);
    dK3 = deg2rad(dK3);
    iK3 = deg2rad(iK3);
    
    num_samples = length(K1);
    
    % 初始化特徵向量矩陣
    V = zeros(3, 3, num_samples);
    
    % 計算方向餘弦(特徵向量)矩陣
    for i = 1:num_samples
        % 計算方向餘弦 (l, m, n)
        l1 = cos(iK1(i)) * cos(dK1(i));
        m1 = cos(iK1(i)) * sin(dK1(i));
        n1 = sin(iK1(i));
        
        l2 = cos(iK2(i)) * cos(dK2(i));
        m2 = cos(iK2(i)) * sin(dK2(i));
        n2 = sin(iK2(i));
        
        l3 = cos(iK3(i)) * cos(dK3(i));
        m3 = cos(iK3(i)) * sin(dK3(i));
        n3 = sin(iK3(i));
        
        % 存儲特徵向量矩陣
        V(:,:,i) = [l1, l2, l3;
                    m1, m2, m3; 
                    n1, n2, n3];
    end
    
    %% 計算應變參數
    % 根據公式4.9計算主偽應變
    ln1pe1 = a .* ((K1 ./ K0) - 1) - b;
    ln1pe2 = a .* ((K2 ./ K0) - 1) - b;
    ln1pe3 = a .* ((K3 ./ K0) - 1) - b;
    
    % 計算有限應變
    e1 = exp(ln1pe1) - 1;
    e2 = exp(ln1pe2) - 1;
    e3 = exp(ln1pe3) - 1;
    
    % 計算體積變化因子
    omega = (1 + e1) .* (1 + e2) .* (1 + e3);
    
    % 初始化主軸座標系下偽應變橢球
    Er = zeros(3, 3, num_samples);
    
    % 計算主軸座標系下偽應變橢球(公式4.10)
    for i = 1:num_samples
        % 計算公式4.10中的幂次
        omega_i = omega(i);
        
        % 構建主軸座標系下偽應變橢球
        Er(:,:,i) = omega_i^(2) * [
            (1 + e1(i))^(-2), 0, 0;
            0, (1 + e2(i))^(-2), 0;
            0, 0, (1 + e3(i))^(-2)
        ];
    end
    
    %% 計算地理座標系下的有限應變橢球(公式4.11)
    % 初始化地理座標系下偽應變橢球
    Eg = zeros(3, 3, num_samples);
    
    for i = 1:num_samples
        % 取得特徵向量矩陣及其轉置
        Vi = V(:,:,i);
        ViT = Vi';
        
        % 計算地理座標系下偽應變橢球: Eg = V^T * Er * V
        Eg(:,:,i) = ViT * Er(:,:,i) * Vi;
    end
    
    %% 計算區間平均的Eg
    [Eg_avg, sample_ranges] = calculate_average_Eg(Eg, interval_count, samples_per_interval);
    
    %% 計算區間之間的應變增量
    % 初始化區間應變增量矩陣
    Einc_intervals = zeros(3, 3, interval_count-1);
    U_intervals = zeros(3, 3, interval_count-1);  % 形變率張量
    R_intervals = zeros(3, 3, interval_count-1);  % 旋轉率張量
    
    % 初始化形變率張量主軸和主值的矩陣
    U_eigenvalues_intervals = zeros(3, interval_count-1);
    U_eigenvectors_intervals = zeros(3, 3, interval_count-1);
    
    for i = 1:interval_count-1
        % 初始有限應變矩陣
        E1 = Eg_avg(:,:,i);
        % 變形後有限應變矩陣
        E2 = Eg_avg(:,:,i+1);
        
        % 計算應變增量: Einc = E2 * E1^(-1)
        Einc_intervals(:,:,i) = E2 * inv(E1);
        
        % 極分解(Polar Decomposition)
        % 獲取應變增量
        F = Einc_intervals(:,:,i);
        
        % 計算 C = F^T * F = U^2
        C = F' * F;
        
        % 計算 U = sqrt(C) = sqrt(F^T * F)
        U_intervals(:,:,i) = sqrtm(C);
        
        % 計算旋轉張量: R = F * U^(-1)
        R_intervals(:,:,i) = F * inv(U_intervals(:,:,i));
        
        % 計算形變率張量的特徵值和特徵向量
        [eigvec, eigval] = eig(U_intervals(:,:,i));
        
        % 按特徵值降序排序
        [eigval_diag, idx] = sort(diag(eigval), 'descend');
        eigvec = eigvec(:, idx);
        
        % 存儲排序後的特徵向量和特徵值
        U_eigenvectors_intervals(:,:,i) = eigvec;
        U_eigenvalues_intervals(:,i) = eigval_diag;
    end
    
    %% 輸出區間應變分析結果
    fprintf('=====================================================\n');
    fprintf('區間增量應變分析結果\n');
    fprintf('=====================================================\n\n');
    
    % 顯示區間劃分
    fprintf('區間劃分：\n');
    for i = 1:interval_count
        fprintf('區間 %d: 樣本 %d 到 %d\n', i, sample_ranges(i, 1), sample_ranges(i, 2));
    end
    fprintf('\n');
    
    % 顯示區間平均Eg矩陣
    fprintf('區間平均有限應變橢球(Eg)：\n');
    for i = 1:interval_count
        fprintf('區間 %d:\n', i);
        disp(Eg_avg(:,:,i));
    end
    fprintf('\n');
    
    % 顯示區間間的應變增量
    fprintf('相鄰區間之間的應變增量(Einc)：\n');
    for i = 1:interval_count-1
        fprintf('區間 %d 到 區間 %d:\n', i, i+1);
        disp(Einc_intervals(:,:,i));
    end
    fprintf('\n');
    
    % 顯示區間間的形變率張量和旋轉率張量
    fprintf('相鄰區間之間的形變率張量(U)和旋轉率張量(R)：\n');
    for i = 1:interval_count-1
        fprintf('區間 %d 到 區間 %d:\n', i, i+1);
        fprintf('形變率張量(U):\n');
        disp(U_intervals(:,:,i));
        fprintf('旋轉率張量(R):\n');
        disp(R_intervals(:,:,i));
    end
    fprintf('\n');
    
    % 顯示形變率張量的主軸和主值
    fprintf('相鄰區間之間形變率張量的主軸和主值：\n');
    for i = 1:interval_count-1
        fprintf('區間 %d 到 區間 %d:\n', i, i+1);
        fprintf('主值: [%.4f, %.4f, %.4f]\n', U_eigenvalues_intervals(1,i), U_eigenvalues_intervals(2,i), U_eigenvalues_intervals(3,i));
        fprintf('主軸方向:\n');
        disp(U_eigenvectors_intervals(:,:,i));
    end
    fprintf('\n');
    
    % 繪製結果圖
    figure;
    
    % 形變率張量主值比較
    subplot(2,2,1);
    bar(1:(interval_count-1), U_eigenvalues_intervals');
    title('各區間間形變率張量主值比較');
    xlabel('區間間隔');
    ylabel('形變率主值');
    legend('U1', 'U2', 'U3', 'Location', 'best');
    grid on;
    
    % 輸出分析結果到Excel
    % 創建輸出表格
    results = table();
    results.Interval_Start = sample_ranges(1:end-1, 1);
    results.Interval_End = sample_ranges(1:end-1, 2);
    results.NextInterval_Start = sample_ranges(2:end, 1);
    results.NextInterval_End = sample_ranges(2:end, 2);
    results.U1 = U_eigenvalues_intervals(1, :)';
    results.U2 = U_eigenvalues_intervals(2, :)';
    results.U3 = U_eigenvalues_intervals(3, :)';
    
    % 生成輸出文件名
    [~, fname, ~] = fileparts(filename);
    output_filename = [fname '_interval_analysis.xlsx'];
    
    % 寫入Excel文件
    writetable(results, output_filename);
    fprintf('分析結果已保存至文件: %s\n', output_filename);
end