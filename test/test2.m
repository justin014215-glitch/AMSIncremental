% 讀取 Excel 檔案
filename = 'NCIH.xlsx';
data = readtable(filename, 'VariableNamingRule', 'preserve');

% 讀取 AMS 參數
K1 = data.K1;   % 最大主軸磁感率
K2 = data.K2;   % 中間主軸磁感率
K3 = data.K3;   % 最小主軸磁感率

dK1 = data.dK1geo; iK1 = data.iK1geo; % K1 軸的 Declination & Inclination
dK2 = data.dK2geo; iK2 = data.iK2geo; % K2 軸的 Declination & Inclination
dK3 = data.dK3geo; iK3 = data.iK3geo; % K3 軸的 Declination & Inclination

num_samples = length(K1); % 樣本數
S_matrices = zeros(num_samples, 3, 3); % 存儲所有樣本的 S 矩陣

% 逐筆計算 AMS 矩陣 S
for i = 1:num_samples
    % 計算 K1 軸的方向餘弦
    l1 = cosd(dK1(i)) * cosd(iK1(i));
    m1 = sind(dK1(i)) * cosd(iK1(i));
    n1 = sind(iK1(i));

    % 計算 K2 軸的方向餘弦
    l2 = cosd(dK2(i)) * cosd(iK2(i));
    m2 = sind(dK2(i)) * cosd(iK2(i));
    n2 = sind(iK2(i));

    % 計算 K3 軸的方向餘弦
    l3 = cosd(dK3(i)) * cosd(iK3(i));
    m3 = sind(dK3(i)) * cosd(iK3(i));
    n3 = sind(iK3(i));

    % 構造旋轉矩陣 V
    V = [l1, l2, l3;
         m1, m2, m3;
         n1, n2, n3];

    % AMS 張量 (對角矩陣)
    K_tensor = diag([K1(i), K2(i), K3(i)]);

    % 計算 AMS 張量 S = V * K_tensor * V'
    S = V * K_tensor * V';

    % 存儲結果
    S_matrices(i, :, :) = S;
end

% 顯示第一個樣本的 AMS 矩陣
disp('第一個樣本的 AMS 矩陣 S：');
disp(squeeze(S_matrices(1, :, :)));

% 儲存結果到 Excel（每個樣本 9 個數值，展平為 1 列）
output_filename = 'AMS_Tensor.xlsx';
writematrix(reshape(S_matrices, num_samples, 9), output_filename);
disp(['AMS 矩陣已儲存至 ' output_filename]);
