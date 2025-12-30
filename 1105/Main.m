%% 腳本：從單一檔案驗證磁感率張量 (Kg) 與 應變張量 (Er/Eg)
%
% 修改說明：
% 1. [已更新] V (特徵向量矩陣) 直接在 Main 計算
% 2. [新增] 加入 calculateEr 計算應變張量 (Strain Tensor)
% 3. [新增] 顯示 Er (主軸應變) 與 Eg (地理應變)
% 4. 顯示順序：K_diag -> V -> Kg -> [Er -> Eg] -> FK
clear; clc;

%% ===== 參數設定 =====
% 1. 指定您的資料檔案名稱
DATA_FILENAME = 'test.xlsx';

% 2. 指定要比較的資料行 (Row Index)
INITIAL_ROW_INDEX = 1;  % <--- 指定初始狀態的「行號」
FINAL_ROW_INDEX   = 2;  % <--- 指定最終狀態的「行號」

% 3. Slate 模型參數 (用於 calculateEr)
%    通常 slateA=6.897, slateB=0.007 (依據實驗經驗值)
SLATE_A = 6.897; 
SLATE_B = 0.007;
% =========================

%% 1. 載入資料
try
    data_table = readtable(DATA_FILENAME);
    fprintf('成功讀取檔案：%s\n', DATA_FILENAME);
catch e
    fprintf('讀取檔案時發生錯誤：%s\n', e.message);
    fprintf('請確保 CSV/XLSX 檔案 (%s) 與此腳本在同一資料夾中。\n', DATA_FILENAME);
    return;
end

% 檢查行號是否有效
if INITIAL_ROW_INDEX > height(data_table) || FINAL_ROW_INDEX > height(data_table)
    error('指定的行號 (INITIAL_ROW_INDEX 或 FINAL_ROW_INDEX) 超出了檔案的總行數。');
end

%% 2. 選取要比較的資料
row_initial = data_table(INITIAL_ROW_INDEX, :);
row_final   = data_table(FINAL_ROW_INDEX, :);

fprintf('正在比較：\n');
fprintf('  - 初始狀態 (Row %d): %s\n', INITIAL_ROW_INDEX, row_initial.Spec_Name);
fprintf('  - 最終狀態 (Row %d): %s\n', FINAL_ROW_INDEX, row_final.Spec_Name);

%% 3. 計算初始狀態 (Initial)
% --- 3.1 磁感率部分 ---
K1_i = row_initial.K1;
K2_i = row_initial.K2;
K3_i = row_initial.K3;
K_diag_initial = diag([K1_i, K2_i, K3_i]);

% 提取方向 & 計算 V 矩陣
d1_i = row_initial.dK1geo; i1_i = row_initial.iK1geo;
d2_i = row_initial.dK2geo; i2_i = row_initial.iK2geo;
d3_i = row_initial.dK3geo; i3_i = row_initial.iK3geo;

v1_i = [cosd(i1_i)*cosd(d1_i); cosd(i1_i)*sind(d1_i); sind(i1_i)];
v2_i = [cosd(i2_i)*cosd(d2_i); cosd(i2_i)*sind(d2_i); sind(i2_i)];
v3_i = [cosd(i3_i)*cosd(d3_i); cosd(i3_i)*sind(d3_i); sind(i3_i)];
V_i = [v1_i, v2_i, v3_i];

% 計算地理磁感率張量 Kg
Kg_initial = AMSFormulas.calculateKg(K1_i, K2_i, K3_i, V_i);

% --- 3.2 應變部分 (新增) ---
% 計算主軸應變張量 (Er)
Er_initial = AMSFormulas.calculateEr(K1_i, K2_i, K3_i, SLATE_A, SLATE_B);

% 計算地理座標應變張量 (Eg)
% 公式：Eg = V' * Er * V (或是 V * Er * V'，視定義而定，此處使用 AMSFormulas.calculateEg)
Eg_initial = AMSFormulas.calculateEg(Er_initial, V_i);

%% 4. 計算最終狀態 (Final)
% --- 4.1 磁感率部分 ---
K1_f = row_final.K1;
K2_f = row_final.K2;
K3_f = row_final.K3;
K_diag_final = diag([K1_f, K2_f, K3_f]);

d1_f = row_final.dK1geo; i1_f = row_final.iK1geo;
d2_f = row_final.dK2geo; i2_f = row_final.iK2geo;
d3_f = row_final.dK3geo; i3_f = row_final.iK3geo;

v1_f = [cosd(i1_f)*cosd(d1_f); cosd(i1_f)*sind(d1_f); sind(i1_f)];
v2_f = [cosd(i2_f)*cosd(d2_f); cosd(i2_f)*sind(d2_f); sind(i2_f)];
v3_f = [cosd(i3_f)*cosd(d3_f); cosd(i3_f)*sind(d3_f); sind(i3_f)];
V_f = [v1_f, v2_f, v3_f];

% 計算地理磁感率張量 Kg
Kg_final = AMSFormulas.calculateKg(K1_f, K2_f, K3_f, V_f);

% --- 4.2 應變部分 (新增) ---
% 計算主軸應變張量 (Er)
Er_final = AMSFormulas.calculateEr(K1_f, K2_f, K3_f, SLATE_A, SLATE_B);

% 計算地理座標應變張量 (Eg)
Eg_final = AMSFormulas.calculateEg(Er_final, V_f);


%% 5. 計算增量磁感率張量 (FK)
FK = AMSFormulas.calculateFK(Kg_initial, Kg_final);

%% 6. 顯示結果
disp(' ');
disp('============================================================');
disp('                    計算結果報告');
disp('============================================================');

% --- Part A: 初始狀態 ---
disp('----------- [1. 初始狀態 Initial State] -----------');
disp('A) 原始主磁感率矩陣 (K_diag):');
disp(K_diag_initial);

disp('B) 特徵向量矩陣 (V):');
disp(V_i);

disp('C) 地理座標系磁感率張量 (Kg):');
disp(Kg_initial);

disp('>>> D) 主軸應變張量 (Er) [Slate Model]:');
disp(Er_initial);

disp('>>> E) 地理座標系應變張量 (Eg):');
disp(Eg_initial);

disp(' '); 

% --- Part B: 最終狀態 ---
disp('----------- [2. 最終狀態 Final State] -----------');
disp('A) 原始主磁感率矩陣 (K_diag):');
disp(K_diag_final);

disp('B) 特徵向量矩陣 (V):');
disp(V_f);

disp('C) 地理座標系磁感率張量 (Kg):');
disp(Kg_final);

disp('>>> D) 主軸應變張量 (Er) [Slate Model]:');
disp(Er_final);

disp('>>> E) 地理座標系應變張量 (Eg):');
disp(Eg_final);

disp(' ');

% --- Part C: 增量分析 ---
disp('----------- [3. 增量分析 (FK)] -----------');
disp('增量磁感率張量 (FK = Kg_final * inv(Kg_initial)):');
disp(FK);

% --- Part D: 分解 ---
[U_k, R_k, ~] = AMSFormulas.polarDecomposition(FK);
disp('FK 的右拉伸張量 (Right Stretch Tensor, U_k):');
disp(U_k);
disp('FK 的旋轉張量 (Rotation Tensor, R_k):');
disp(R_k);

disp('============================================================');
%% 7. FK 的特徵值分解
disp(' ');
disp('----------- [4. 增量張量 FK 的特徵值分析] -----------');

[V_FK, D_FK] = eig(FK);

disp('A) FK 的特徵向量矩陣 (Eigenvectors):');
disp(V_FK);

disp('B) FK 的特徵值矩陣 (Eigenvalues):');
disp(D_FK);

disp('============================================================');

%% 8. 驗證模型
disp(' ');
disp('----------- [5. 疊加驗證] -----------');
Kg_calculated_final = FK * Kg_initial;
difference = Kg_calculated_final - Kg_final;
diff_val = max(max(abs(difference)));

fprintf('驗證誤差 (最大絕對誤差): %.4e\n', diff_val);
if diff_val < 1e-10
    disp('   => 驗證成功！');
else
    disp('   => 驗證警告：存在數值誤差。');
end
disp('============================================================');

%% 9. 3D 繪圖展示
disp(' ');
disp('----------- [繪圖展示] -----------');
% ... (繪圖代碼保持不變，主要顯示磁感率橢球)
% 若要繪製應變橢球，可將傳入 plotEllipsoid3D 的參數改為 Er 的特徵值
plotEllipsoid3D(K1_i, K2_i, K3_i, V_i, ['Initial Susceptibility (Row ' num2str(INITIAL_ROW_INDEX) ')']);
plotEllipsoid3D(K1_f, K2_f, K3_f, V_f, ['Final Susceptibility (Row ' num2str(FINAL_ROW_INDEX) ')']);

lambda = diag(D_FK); 
plotEllipsoid3D(lambda(1), lambda(2), lambda(3), V_FK, 'Incremental Tensor (FK)');