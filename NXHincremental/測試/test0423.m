clear all
close all
clc

%% 讀取 Excel 檔案
filename = 'NXHAMStest.xlsx';
data = readtable(filename, 'VariableNamingRule', 'preserve');

%% 擷取並計算基本磁感率參數
K1 = data.K1;   % 最大主軸磁感率
K2 = data.K2;   % 中間主軸磁感率
K3 = data.K3;   % 最小主軸磁感率
Km = data.Km;   % 平均磁感率
F = data.F;     % 磁性葉理強度
L = data.L;     % 磁性線理強度
T = data.T;     % 形狀參數
P = data.P;     % 異向性程度
Pj = data.Pj;   % 校正異向性

% 計算幾何平均磁感率 K0
K0 = (K1 .* K2 .* K3).^(1/3);
% 計算變形強度 Int
Int = sqrt((L-1).^2+(F-1).^2);

% 將新計算的參數添加到資料表
data.K0 = K0;
data.Int = Int;

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

%% 計算應變增量(公式4.14)
% 初始化應變增量矩陣
Einc = zeros(3, 3, num_samples-1);

for i = 1:num_samples-1
    % 初始有限應變矩陣
    E1 = Eg(:,:,i);
    % 變形後有限應變矩陣
    E3 = Eg(:,:,i+1);
    
    % 計算應變增量: E2 = E3 * E1^(-1)
    Einc(:,:,i) = E3 * inv(E1);
end

%% 極分解(Polar Decomposition) - 符合圖示方法
% 初始化形變率張量和旋轉率張量
U = zeros(3, 3, num_samples-1);  % 形變率張量
R = zeros(3, 3, num_samples-1);  % 旋轉率張量

for i = 1:num_samples-1
    % 獲取第i個應變增量
    F = Einc(:,:,i);
    
    % 計算 C = F^T * F = U^2
    C = F' * F;
    
    % 計算 U = sqrt(C) = sqrt(F^T * F)
    U(:,:,i) = sqrtm(C);
    
    % 計算旋轉張量: R = F * U^(-1)
    R(:,:,i) = F * inv(U(:,:,i));
    
    % 驗證 F = R * U
    verify_F = R(:,:,i) * U(:,:,i);
    if norm(verify_F - F, 'fro') > 1e-10
        warning('樣本 %d 的極分解驗證失敗', i);
    end
end

%% 計算形變率張量主軸
% 初始化存儲主軸和主值的矩陣
U_eigenvectors = zeros(3, 3, num_samples-1);
U_eigenvalues = zeros(3, num_samples-1);

for i = 1:num_samples-1
    % 計算形變率張量的特徵值和特徵向量
    [eigvec, eigval] = eig(U(:,:,i));
    
    % 按特徵值降序排序
    [eigval_diag, idx] = sort(diag(eigval), 'descend');
    eigvec = eigvec(:, idx);
    
    % 存儲排序後的特徵向量和特徵值
    U_eigenvectors(:,:,i) = eigvec;
    U_eigenvalues(:,i) = eigval_diag;
end

%% 將形變率張量從地理座標系旋轉到主軸座標系
% 初始化主軸座標系下的形變率張量
U_principal = zeros(3, 3, num_samples-1);

for i = 1:num_samples-1
    % 獲取特徵向量矩陣
    V_U = U_eigenvectors(:,:,i);
    V_U_T = V_U';
    
    % 創建對角矩陣，包含特徵值 (Er)
    Er_U = diag(U_eigenvalues(:,i));
    
    % 將形變率張量旋轉到主軸座標系 U = V^T * Er * V
    U_principal(:,:,i) = V_U_T * U(:,:,i) * V_U;
    
    % 驗證主軸座標系下的U是否接近對角矩陣
    diff_diag = U_principal(:,:,i) - Er_U;
    off_diag_sum = sum(sum(abs(diff_diag - diag(diag(diff_diag)))));
    
    if off_diag_sum > 1e-10
        warning('樣本 %d 的主軸座標系轉換可能有誤 (非對角元素總和: %.6e)', i, off_diag_sum);
    end
end

%% 結果展示
fprintf('=====================================================\n');
fprintf('磁感率異向性分析與應變增量計算\n');
fprintf('=====================================================\n\n');

% 顯示樣本數量
fprintf('總樣本數: %d\n\n', num_samples);

% 顯示前幾個樣本的主要磁感率參數
fprintf('前3個樣本的磁感率參數:\n');
fprintf('%-10s %-10s %-10s %-10s %-10s %-10s\n', 'Sample', 'K1', 'K2', 'K3', 'K0', 'Pj');
for i = 1:min(3, num_samples)
    fprintf('%-10d %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f\n', i, K1(i), K2(i), K3(i), K0(i), Pj(i));
end
fprintf('\n');

% 顯示前幾個樣本的有限應變橢球
fprintf('前3個樣本的有限應變橢球(Eg):\n');
for i = 1:min(3, num_samples)
    fprintf('Sample %d:\n', i);
    disp(Eg(:,:,i));
end
fprintf('\n');

% 顯示前幾個樣本間的應變增量
fprintf('前2個相鄰樣本間的應變增量(Einc):\n');
for i = 1:min(2, num_samples-1)
    fprintf('Samples %d to %d:\n', i, i+1);
    disp(Einc(:,:,i));
end
fprintf('\n');

% 顯示前幾個樣本間的形變率張量和旋轉率張量
fprintf('前2個相鄰樣本間的形變率張量(U)和旋轉率張量(R):\n');
for i = 1:min(2, num_samples-1)
    fprintf('Samples %d to %d:\n', i, i+1);
    fprintf('形變率張量(U):\n');
    disp(U(:,:,i));
    fprintf('旋轉率張量(R):\n');
    disp(R(:,:,i));
end
fprintf('\n');

% 顯示前幾個樣本間形變率張量的主軸和主值
fprintf('前2個相鄰樣本間形變率張量的主軸和主值:\n');
for i = 1:min(2, num_samples-1)
    fprintf('Samples %d to %d:\n', i, i+1);
    fprintf('主值: [%.4f, %.4f, %.4f]\n', U_eigenvalues(1,i), U_eigenvalues(2,i), U_eigenvalues(3,i));
    fprintf('主軸方向:\n');
    disp(U_eigenvectors(:,:,i));
end
fprintf('\n');

% 顯示主軸座標系下的形變率張量
fprintf('前2個相鄰樣本間主軸座標系下的形變率張量:\n');
for i = 1:min(2, num_samples-1)
    fprintf('Samples %d to %d:\n', i, i+1);
    disp(U_principal(:,:,i));
end
fprintf('\n');

% 繪製結果圖
figure;
subplot(2,2,1);
plot(1:num_samples, e1, 'r-o', 1:num_samples, e2, 'g-s', 1:num_samples, e3, 'b-^');
title('主應變值隨樣本變化');
xlabel('樣本編號');
ylabel('主應變值');
legend('e1', 'e2', 'e3', 'Location', 'best');
grid on;

% 如果有多個樣本，顯示形變率張量主值隨樣本變化
if num_samples > 1
    subplot(2,2,2);
    plot(1.5:num_samples-0.5, U_eigenvalues(1,:), 'r-o', 1.5:num_samples-0.5, U_eigenvalues(2,:), 'g-s', 1.5:num_samples-0.5, U_eigenvalues(3,:), 'b-^');
    title('形變率張量主值隨樣本變化');
    xlabel('樣本間隔');
    ylabel('形變率主值');
    legend('U1', 'U2', 'U3', 'Location', 'best');
    grid on;
end

% 磁感率橢球形狀和應變橢球形狀的對比
subplot(2,2,3);
plot(Pj, T, 'ko', 'MarkerFaceColor', 'b');
title('磁感率橢球形狀參數');
xlabel('校正異向性 (Pj)');
ylabel('形狀參數 (T)');
grid on;

% 保存結果到 Excel 檔案
results = struct();
results.K0 = K0;
results.Int = Int;
results.e1 = e1;
results.e2 = e2;
results.e3 = e3;
results.omega = omega;

% 添加形變率張量主值到結果
if num_samples > 1
    for i = 1:min(num_samples-1, 10)  % 最多保存前10個間隔的結果
        field_name = ['U1_' num2str(i) '_' num2str(i+1)];
        results.(field_name) = U_eigenvalues(1,i);
        field_name = ['U2_' num2str(i) '_' num2str(i+1)];
        results.(field_name) = U_eigenvalues(2,i);
        field_name = ['U3_' num2str(i) '_' num2str(i+1)];
        results.(field_name) = U_eigenvalues(3,i);
    end
end

