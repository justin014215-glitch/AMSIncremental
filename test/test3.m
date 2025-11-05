%% 讀取 Excel 檔案
filename = 'NCIH.xlsx';
data = readtable(filename, 'VariableNamingRule', 'preserve');

% 讀取 AMS 參數 (磁感率主軸 & 主軸方向)
K1 = data.K1;  K2 = data.K2;  K3 = data.K3;  % 磁感率主軸
dK1 = data.dK1geo; iK1 = data.iK1geo;  % K1 方位角 & 傾角
dK2 = data.dK2geo; iK2 = data.iK2geo;  % K2 方位角 & 傾角
dK3 = data.dK3geo; iK3 = data.iK3geo;  % K3 方位角 & 傾角

num_samples = length(K1); % 樣本數 (應為 124)

% 參數設定 (根據 Wood et al. 1976)
a = 6.897;
b = 0.007;

% 初始化變數
S_matrices = zeros(num_samples, 3, 3); % AMS 矩陣
E_s_matrices = zeros(num_samples, 3, 3); % 有限應變張量

%% 逐筆計算 AMS 張量 -> 有限應變張量
for i = 1:num_samples
    % 計算 K1, K2, K3 軸的方向餘弦
    l1 = cosd(dK1(i)) * cosd(iK1(i));  m1 = sind(dK1(i)) * cosd(iK1(i));  n1 = sind(iK1(i));
    l2 = cosd(dK2(i)) * cosd(iK2(i));  m2 = sind(dK2(i)) * cosd(iK2(i));  n2 = sind(iK2(i));
    l3 = cosd(dK3(i)) * cosd(iK3(i));  m3 = sind(dK3(i)) * cosd(iK3(i));  n3 = sind(iK3(i));

    % 組成旋轉矩陣 V
    V = [l1, l2, l3;
         m1, m2, m3;
         n1, n2, n3];

    % AMS 張量 (對角矩陣)
    D = diag([K1(i), K2(i), K3(i)]);

    % 計算 AMS 矩陣 S = V * D * V'
    S = V * D * V';
    S_matrices(i, :, :) = S;

    % 計算 AMS 張量的特徵值與特徵向量
    [V_eigen, D_eigen] = eig(S);  
    lambda = diag(D_eigen); % 特徵值

    % 計算主有限應變值 e1, e2, e3
    e = log(lambda / (prod(lambda)^(1/3))) / a - b;
    
    % 組成有限應變張量 E_s
    E_s = diag(1 + e);
    E_s_matrices(i, :, :) = V_eigen * E_s * V_eigen'; % 旋轉至地理座標系
end

%% 儲存結果到 Excel
output_filename = 'Finite_Strain_Tensor.xlsx';
E_s_output = reshape(E_s_matrices, num_samples, 9);
writematrix(E_s_output, output_filename);
disp(['有限應變張量已儲存至 ' output_filename]);

disp(data);

%% 4.3.2 應變增量計算 (分期：每 31 顆樣本為一期)
% 假設 E_s_matrices 已從 4.3.1 計算完，尺寸為 [124 x 3 x 3]

num_samples = size(E_s_matrices, 1);
group_size = 31;
num_groups = num_samples / group_size;  % 這裡假設可整除

% 初始化存放應變增量及其對稱、反對稱部分
E_inc_matrices = zeros(num_samples, 3, 3);
E_inc_sym = zeros(num_samples, 3, 3); % 對稱部分 (形變率張量 D*)
E_inc_asym = zeros(num_samples, 3, 3); % 反對稱部分 (旋轉率張量 W*)

for g = 1:num_groups
    % 取得該期的樣本索引
    idx_start = (g - 1) * group_size + 1;
    idx_end = g * group_size;
    
    % 使用該期第一顆樣本作為初始有限應變張量 E_initial
    E_initial = squeeze(E_s_matrices(idx_start, :, :));
    
    % 逐顆計算該期內的應變增量
    for i = idx_start:idx_end
        E_final = squeeze(E_s_matrices(i, :, :)); % 該樣本有限應變張量
        % 根據公式：E_inc * E_initial = E_final  =>  E_inc = E_final * inv(E_initial)
        E_inc = E_final * inv(E_initial);
        E_inc_matrices(i, :, :) = E_inc;
        
        % 將 E_inc 分解為對稱部分 (D*) 與反對稱部分 (W*)
        E_inc_sym(i, :, :) = (E_inc + E_inc') / 2;
        E_inc_asym(i, :, :) = (E_inc - E_inc') / 2;
    end
end

% 顯示第一期（第一 31 顆樣本）第一顆樣本的結果
disp('第一期第一顆樣本的應變增量矩陣 E_inc：');
disp(squeeze(E_inc_matrices(1, :, :)));
disp('第一期第一顆樣本的對稱部分 D*：');
disp(squeeze(E_inc_sym(1, :, :)));
disp('第一期第一顆樣本的反對稱部分 W*：');
disp(squeeze(E_inc_asym(1, :, :)));

% 儲存所有期的應變增量矩陣到 Excel（每筆資料展平為 9 個數值）
output_filename_inc = 'Incremental_Strain_Tensor_Grouped.xlsx';
E_inc_output = reshape(E_inc_matrices, num_samples, 9);
writematrix(E_inc_output, output_filename_inc);
disp(['應變增量矩陣已儲存至 ' output_filename_inc]);
