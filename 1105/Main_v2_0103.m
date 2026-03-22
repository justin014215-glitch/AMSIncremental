%% 腳本：多樣本群組分析 + 指定樣本增量分析 (v4)
% 
% 功能：
% 1. [批次] 計算所有群組的平均張量 (Mean Tensor) 與 Jelinek 參數
% 2. [指定] 選擇兩組樣本，計算增量應變 (Incremental Strain)
% 3. [繪圖] 繪製初始、最終與增量橢球
%
% 輸出：
% - 完整的 AMS 平均數據表
% - 正推 (Forward) 與 逆推 (Inverse) 的增量張量分析
% - 3D 視覺化圖形

clear; clc; close all;

%% ===== 1. 參數設定 (User Settings) =====
FILENAME = 'Debby_rawdata.xlsx'; % 資料檔名

% Slate 模型參數 (用於將 K 轉為應變)
SLATE_A = 6.897; 
SLATE_B = 0.007;

% --- [設定] 指定要進行「增量分析」的兩個群組編號 (No.) ---
ID_INITIAL = 1;  % 初始狀態 (Initial State) 的 No.
ID_FINAL   = 4;  % 最終狀態 (Final State) 的 No.
% ==============================================

%% ===== 2. 讀取與前處理 =====
fprintf('正在讀取檔案...\n');
try
    raw_table = readtable(FILENAME, 'VariableNamingRule', 'preserve');
catch
    warning('無法使用 preserve 模式，嘗試使用預設模式讀取...');
    raw_table = readtable(FILENAME); 
end

% 自動偵測樣本編號欄位
possible_id_names = {'Number', 'Number.', 'No.', 'No', 'No_', 'Number_'};
col_ID = '';
var_names = raw_table.Properties.VariableNames;
for i = 1:length(possible_id_names)
    if ismember(possible_id_names{i}, var_names)
        col_ID = possible_id_names{i};
        break;
    end
end
if isempty(col_ID), error('找不到樣本編號欄位。'); end

unique_ids = unique(raw_table.(col_ID));

% 準備儲存容器
results_struct = struct();
tensor_db = containers.Map('KeyType', 'double', 'ValueType', 'any'); % 用來存 3x3 張量
result_count = 0;

fprintf('共 %d 個群組，開始計算平均張量...\n', length(unique_ids));

%% ===== 3. 群組平均計算 (Batch Processing) =====
for i = 1:length(unique_ids)
    curr_id = unique_ids(i);
    sub_table = raw_table(raw_table.(col_ID) == curr_id, :);
    num_subs = height(sub_table);
    
    tensor_stack_norm = zeros(3, 3, num_subs);
    
    % --- 內層迴圈：處理單一標本 ---
    for j = 1:num_subs
        row = sub_table(j, :);
        % 提取 K 值與方向
        K1 = row.K1; K2 = row.K2; K3 = row.K3;
        % 自動對應欄位 (若名稱不同請在此修改)
        d1 = row.dK1geo; i1 = row.iK1geo;
        d2 = row.dK2geo; i2 = row.iK2geo;
        d3 = row.dK3geo; i3 = row.iK3geo;
        
        % 計算 V
        v1 = [cosd(i1)*cosd(d1); cosd(i1)*sind(d1); sind(i1)];
        v2 = [cosd(i2)*cosd(d2); cosd(i2)*sind(d2); sind(i2)];
        v3 = [cosd(i3)*cosd(d3); cosd(i3)*sind(d3); sind(i3)];
        V = [v1, v2, v3];
        
        % 計算 Kg 並標準化
        Kg = AMSFormulas.calculateKg(K1, K2, K3, V);
        Kg_norm = AMSFormulas.standardizeTensor(Kg);
        tensor_stack_norm(:, :, j) = Kg_norm;
    end
    
    % --- 計算平均張量 (Mean Tensor) ---
    Mean_Kg = AMSFormulas.calculateMeanTensor(tensor_stack_norm);
    
    % [重要] 將平均張量存入 Map，方便稍後提取
    tensor_db(curr_id) = Mean_Kg;
    
    % --- 特徵分析 ---
    mean_analysis = AMSFormulas.performEigenAnalysis(Mean_Kg);
    
    % --- 應變計算 (Slate) ---
    K_mean_vals = mean_analysis.vals;
    Er_mean = AMSFormulas.calculateEr(K_mean_vals(1), K_mean_vals(2), K_mean_vals(3), SLATE_A, SLATE_B);
    Eg_mean = AMSFormulas.calculateEg(Er_mean, mean_analysis.V);
    strain_analysis = AMSFormulas.performEigenAnalysis(Eg_mean);
    
    % --- 存入結果表 ---
    result_count = result_count + 1;
    results_struct(result_count).ID = curr_id;
    results_struct(result_count).N  = num_subs;
    
    % K1, K2, K3 (平均張量特徵值)
    results_struct(result_count).K1_Val = mean_analysis.vals(1);
    results_struct(result_count).K2_Val = mean_analysis.vals(2);
    results_struct(result_count).K3_Val = mean_analysis.vals(3);
    
    % 方向 (Trend/Plunge)
    results_struct(result_count).K1_Tr = mean_analysis.K1_dir(1); results_struct(result_count).K1_Pl = mean_analysis.K1_dir(2);
    results_struct(result_count).K2_Tr = mean_analysis.K2_dir(1); results_struct(result_count).K2_Pl = mean_analysis.K2_dir(2);
    results_struct(result_count).K3_Tr = mean_analysis.K3_dir(1); results_struct(result_count).K3_Pl = mean_analysis.K3_dir(2);
    
    % 形狀因子
    results_struct(result_count).L = mean_analysis.L;
    results_struct(result_count).F = mean_analysis.F;
    results_struct(result_count).P = mean_analysis.P;
    
    % 應變軸
    results_struct(result_count).S_X = strain_analysis.vals(1);
    results_struct(result_count).S_Y = strain_analysis.vals(2);
    results_struct(result_count).S_Z = strain_analysis.vals(3);
end

% 顯示簡易列表
Final_Table = struct2table(results_struct);
disp(' ');
disp('=== 所有群組計算完成 (前5筆) ===');
if height(Final_Table)>5, disp(Final_Table(1:5,:)); else, disp(Final_Table); end


%% ===== 4. 指定增量分析 (Incremental Analysis) =====
disp(' ');
disp('============================================================');
fprintf('                增量分析: No.%d -> No.%d\n', ID_INITIAL, ID_FINAL);
disp('============================================================');

% 檢查 ID 是否存在
if ~isKey(tensor_db, ID_INITIAL) || ~isKey(tensor_db, ID_FINAL)
    error('錯誤：指定的 ID (%d 或 %d) 不在資料庫中，請檢查編號。', ID_INITIAL, ID_FINAL);
end

% 提取張量
Kg_init = tensor_db(ID_INITIAL);
Kg_final = tensor_db(ID_FINAL);

% 提取特徵 (為了顯示和繪圖)
Analysis_Init = AMSFormulas.performEigenAnalysis(Kg_init);
Analysis_Final = AMSFormulas.performEigenAnalysis(Kg_final);

% --- A. 正推增量 (Forward): FK = Final * inv(Initial) ---
FK_forward = AMSFormulas.calculateFK(Kg_init, Kg_final);
Analysis_FK_Fwd = AMSFormulas.performEigenAnalysis(FK_forward);

% --- B. 逆推增量 (Inverse): FK_inv = Initial * inv(Final) ---
% 代表「如何從最終變回初始」的變形
FK_inverse = AMSFormulas.calculateFK(Kg_final, Kg_init); 
Analysis_FK_Inv = AMSFormulas.performEigenAnalysis(FK_inverse);

% --- 顯示報告 ---
disp('--- [1] 正推增量 (Forward Increment: Initial -> Final) ---');
disp('特徵值 (Principal Values of Increment):');
disp(Analysis_FK_Fwd.vals');
disp('主軸方向 (Principal Directions):');
fprintf('  Max: %.1f / %.1f\n', Analysis_FK_Fwd.K1_dir(1), Analysis_FK_Fwd.K1_dir(2));
fprintf('  Min: %.1f / %.1f\n', Analysis_FK_Fwd.K3_dir(1), Analysis_FK_Fwd.K3_dir(2));
disp('張量矩陣 (Tensor):');
disp(FK_forward);

disp(' ');
disp('--- [2] 逆推增量 (Inverse Increment: Final -> Initial) ---');
disp('特徵值 (Principal Values):');
disp(Analysis_FK_Inv.vals');
disp('主軸方向 (Principal Directions):');
fprintf('  Max: %.1f / %.1f\n', Analysis_FK_Inv.K1_dir(1), Analysis_FK_Inv.K1_dir(2));
fprintf('  Min: %.1f / %.1f\n', Analysis_FK_Inv.K3_dir(1), Analysis_FK_Inv.K3_dir(2));


%% ===== 5. 3D 繪圖 (Visualization) =====
disp(' ');
disp('正在繪製 3D 橢球圖...');

% 設定圖表排列
figure('Name', 'Incremental Strain Analysis', 'Position', [100, 100, 1200, 400]);

% 1. 初始狀態
subplot(1, 3, 1);
plotEllipsoid3D_subplot(Analysis_Init.vals, Analysis_Init.V, ...
    sprintf('Initial State (No.%d)', ID_INITIAL));

% 2. 增量 (正推)
subplot(1, 3, 2);
plotEllipsoid3D_subplot(Analysis_FK_Fwd.vals, Analysis_FK_Fwd.V, ...
    'Increment (Initial -> Final)');

% 3. 最終狀態
subplot(1, 3, 3);
plotEllipsoid3D_subplot(Analysis_Final.vals, Analysis_Final.V, ...
    sprintf('Final State (No.%d)', ID_FINAL));

disp('完成！');


%% ===== 輔助函式: 支援 subplot 的繪圖 =====
function plotEllipsoid3D_subplot(vals, V, titleText)
    % 這是 plotEllipsoid3D 的簡化版，專門用於 subplot
    a = vals(1); b = vals(2); c = vals(3);
    
    % 產生網格
    [u, v] = meshgrid(linspace(0, 2*pi, 30), linspace(0, pi, 15));
    x = a * cos(u) .* sin(v);
    y = b * sin(u) .* sin(v);
    z = c * cos(v);
    
    % 旋轉
    pts = V * [x(:)'; y(:)'; z(:)'];
    x_rot = reshape(pts(1, :), size(x));
    y_rot = reshape(pts(2, :), size(y));
    z_rot = reshape(pts(3, :), size(z));
    
    % 繪圖
    surf(x_rot, y_rot, z_rot, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'FaceColor', 'cyan');
    hold on; axis equal; grid on;
    set(gca, 'XDir', 'reverse'); % 地質習慣 X 反向
    
    % 畫軸
    max_r = max(vals)*1.3;
    plot3([-max_r max_r], [0 0], [0 0], 'k--'); text(max_r,0,0,'X');
    plot3([0 0], [-max_r max_r], [0 0], 'k--'); text(0,max_r,0,'Y');
    plot3([0 0], [0 0], [-max_r max_r], 'k--'); text(0,0,max_r,'Z');
    
    % 畫主軸箭頭
    quiver3(0,0,0, V(1,1)*a, V(2,1)*a, V(3,1)*a, 'r', 'LineWidth', 2); % K1
    quiver3(0,0,0, V(1,3)*c, V(2,3)*c, V(3,3)*c, 'b', 'LineWidth', 2); % K3
    
    title(titleText);
    xlabel('N'); ylabel('E'); zlabel('D');
    view(135, 30);
    camlight; lighting gouraud;
    hold off;
end