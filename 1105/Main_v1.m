%% 腳本：從單一檔案驗證磁感率張量 (Kg)
%
% 修改說明：
% 1. [已更新] V (特徵向量矩陣) 直接在 Main 計算
% 2. [新增] 在結果報告中加入 V 矩陣的顯示
% 3. 顯示順序：K_diag (主軸) -> V (方向) -> Kg (地理) -> FK (增量)
% 4. 加上初始、最終、合併、應變增量橢球之3D圖
clear; clc;

%% ===== 參數設定 =====
% 1. 指定您的資料檔案名稱
DATA_FILENAME = 'test.xlsx';

% 2. 指定要比較的資料行 (Row Index)
INITIAL_ROW_INDEX = 7;  % <--- 指定初始狀態的「行號」(北東45度橢圓)
FINAL_ROW_INDEX   = 8;  % <--- 指定最終狀態的「行號」(北西45度橢圓)
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

%% 3. 計算初始狀態的 Kg (Kg_initial)
% 提取主磁化率
K1_i = row_initial.K1;
K2_i = row_initial.K2;
K3_i = row_initial.K3;

% --- [顯示用] 建立主軸對角矩陣 ---
K_diag_initial = diag([K1_i, K2_i, K3_i]);

% 提取主軸方向 (直接使用 Excel 中的角度 Degree)
d1_i = row_initial.dK1geo; i1_i = row_initial.iK1geo;
d2_i = row_initial.dK2geo; i2_i = row_initial.iK2geo;
d3_i = row_initial.dK3geo; i3_i = row_initial.iK3geo;

% [修改] 直接在 Main 中計算 V 矩陣 (使用 AMSFormulas 的邏輯)
% 將角度轉換為單位向量 (方向餘弦)
v1_i = [cosd(i1_i)*cosd(d1_i); cosd(i1_i)*sind(d1_i); sind(i1_i)]; % K1 軸向量
v2_i = [cosd(i2_i)*cosd(d2_i); cosd(i2_i)*sind(d2_i); sind(i2_i)]; % K2 軸向量
v3_i = [cosd(i3_i)*cosd(d3_i); cosd(i3_i)*sind(d3_i); sind(i3_i)]; % K3 軸向量

% 組合為 3x3 特徵向量矩陣 V = [v1, v2, v3]
V_i = [v1_i, v2_i, v3_i];

% 計算 Kg (仍然使用公式庫進行張量轉換)
Kg_initial = AMSFormulas.calculateKg(K1_i, K2_i, K3_i, V_i);

%% 4. 計算最終狀態的 Kg (Kg_final)
% 提取主磁化率
K1_f = row_final.K1;
K2_f = row_final.K2;
K3_f = row_final.K3;

% --- [顯示用] 建立主軸對角矩陣 ---
K_diag_final = diag([K1_f, K2_f, K3_f]);

% 提取主軸方向
d1_f = row_final.dK1geo; i1_f = row_final.iK1geo;
d2_f = row_final.dK2geo; i2_f = row_final.iK2geo;
d3_f = row_final.dK3geo; i3_f = row_final.iK3geo;

% [修改] 直接在 Main 中計算 V 矩陣
v1_f = [cosd(i1_f)*cosd(d1_f); cosd(i1_f)*sind(d1_f); sind(i1_f)];
v2_f = [cosd(i2_f)*cosd(d2_f); cosd(i2_f)*sind(d2_f); sind(i2_f)];
v3_f = [cosd(i3_f)*cosd(d3_f); cosd(i3_f)*sind(d3_f); sind(i3_f)];

V_f = [v1_f, v2_f, v3_f];

% 計算 Kg
Kg_final = AMSFormulas.calculateKg(K1_f, K2_f, K3_f, V_f);


%% 5. 計算增量磁感率張量 (FK)
FK = AMSFormulas.calculateFK(Kg_initial, Kg_final);

%% 6. 顯示結果 (依需求調整順序)
disp(' ');
disp('============================================================');
disp('                    計算結果報告');
disp('============================================================');

% --- Part A: 初始狀態 ---
disp('----------- [1. 初始狀態 Initial State] -----------');
disp('A) 原始主磁感率矩陣 (主軸座標系, K_diag):');
disp('   [說明: 這是 Excel 中原始的 K1, K2, K3]');
disp(K_diag_initial);

disp('B) 特徵向量矩陣 (Eigen Vector Matrix, V):');
disp('   [說明: 由 dK, iK 計算出的方向餘弦矩陣');%, Cols=[v1,v2,v3]]
disp(V_i);

disp('C) 地理座標系磁感率張量 (Geographic, Kg):');
disp('   [說明: 經由 V * K_diag * V'' 轉換後]');
disp(Kg_initial);

disp(' '); % 空行分隔

% --- Part B: 最終狀態 ---
disp('----------- [2. 最終狀態 Final State] -----------');
disp('A) 原始主磁感率矩陣 (主軸座標系, K_diag):');
disp('   [說明: 這是 Excel 中原始的 K1, K2, K3]');
disp(K_diag_final);

disp('B) 特徵向量矩陣 (Eigen Vector Matrix, V):');
disp('   [說明: 由 dK, iK 計算出的方向餘弦矩陣'); %, Cols=[v1,v2,v3]]
disp(V_f);

disp('C) 地理座標系磁感率張量 (Geographic, Kg):');
disp('   [說明: 經由 V * K_diag * V'' 轉換後]');
disp(Kg_final);

disp(' '); % 空行分隔

% --- Part C: 合併與增量 ---
disp('----------- [3. 增量與合併分析] -----------');
disp('A) 合併磁感率張量 (Kg_initial * Kg_final):');
disp(Kg_initial * Kg_final);

disp('B) 增量磁感率張量 (FK = Kg_final * inv(Kg_initial)):');
disp(FK);

% --- Part D: 分解 ---
%[U_k, R_k, ~] = AMSFormulas.polarDecomposition(FK);
%disp('C) FK 的右拉伸張量 (Right Stretch Tensor, U_k):');
%disp(U_k);
%disp('D) FK 的旋轉張量 (Rotation Tensor, R_k):');
%disp(R_k);

disp('============================================================');
%% 7. FK 的特徵值分解 (Eigen Analysis of FK)
disp(' ');
disp('----------- [4. 增量張量 FK 的特徵值分析] -----------');

% 計算 FK 的特徵向量 (V_FK) 與特徵值 (D_FK)
[V_FK, D_FK] = eig(FK);

disp('A) FK 的特徵向量矩陣 (Eigenvectors of FK):');
disp('   [說明: 每一行 (Column) 代表一個特徵向量]');
disp(V_FK);

disp('B) FK 的特徵值矩陣 (Eigenvalues of FK):');
disp('   [說明: 對角線上的數值即為特徵值]');
disp(D_FK);

disp('============================================================');

%% 8. 驗證模型：增量 * 初始 = 最終 (Validation: Increment * Initial = Final)
disp(' ');
disp('----------- [5. 疊加驗證：增量 * 初始 = 最終] -----------');

% 公式：Kg_final_calc = FK * Kg_initial
Kg_calculated_final = FK * Kg_initial;

disp('A) 計算出的最終磁感率張量 (Calculated Final = FK * Initial):');
disp('   [說明: 由增量張量與初始張量相乘得到的結果]');
disp(Kg_calculated_final);

disp('B) 原始資料的最終磁感率張量 (Original Final Kg):');
disp('   [說明: 直接從 Excel 數據計算出的結果，用於比對]');
disp(Kg_final);

% 計算兩者差異（驗證精確度）
difference = Kg_calculated_final - Kg_final;
diff_val = max(max(abs(difference)));

fprintf('C) 驗證誤差 (最大絕對誤差): %.4e\n', diff_val);
if diff_val < 1e-10
    disp('   => 驗證成功！兩者完全一致。');
else
    disp('   => 驗證警告：存在些微數值誤差。');
end
disp('============================================================');

%% 9. 3D 繪圖展示 (Visualization)
disp(' ');
disp('----------- [繪圖展示] -----------');
disp('正在繪製 3D 橢球圖...');

% 1. 繪製初始狀態 (Initial State)
% 使用變數：K1_i, K2_i, K3_i (主軸長度) 和 V_i (方向矩陣)
plotEllipsoid3D(K1_i, K2_i, K3_i, V_i, ['Initial State (Row ' num2str(INITIAL_ROW_INDEX) ')']);

% 2. 繪製最終狀態 (Final State)
% 使用變數：K1_f, K2_f, K3_f (主軸長度) 和 V_f (方向矩陣)
plotEllipsoid3D(K1_f, K2_f, K3_f, V_f, ['Final State (Row ' num2str(FINAL_ROW_INDEX) ')']);

% 3. 繪製合併磁感率張量
Fmerge = Kg_initial * Kg_final;
[V_Fmerge, D_Fmerge] = eig(Fmerge);
lambda_merge = diag(D_Fmerge);
[lambda_merge_sorted, sort_idx] = sort(lambda_merge, 'descend');
V_Fmerge_sorted = V_Fmerge(:, sort_idx);

% 4. 繪製合併磁感率張量 (Merged Tensor)
plotEllipsoid3D(lambda_merge_sorted(1), lambda_merge_sorted(2), lambda_merge_sorted(3), V_Fmerge_sorted, 'Merged Tensor (Kg_initial * Kg_final)');


% 5. 取得特徵值與向量
[V_FK, D_FK] = eig(FK);
lambda = diag(D_FK);

% 6. [新增] 排序：從大到小 (讓 K1 為最大值)
[lambda_sorted, sort_idx] = sort(lambda, 'descend'); % 2.0, 1.0, 0.5
V_FK_sorted = V_FK(:, sort_idx); % 對應調整向量順序

% 7. 繪圖 (使用排序後的數據)
% 這樣 lambda(1) 就是 2.0 (Max)，對應紅色
% lambda(3) 就是 0.5 (Min)，對應藍色
plotEllipsoid3D(lambda_sorted(1), lambda_sorted(2), lambda_sorted(3), V_FK_sorted, 'Incremental Tensor (FK)');

