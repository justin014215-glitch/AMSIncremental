%% 腳本：多樣本群組分析 (標準化流程 v6)
%
% 執行步驟：
% 1. [Geo] 轉座標系統 (Specimen -> Geographic)
% 2. [Avg] Jelinek 標準化與平均計算
% 3. [Strain] 轉換為應變張量 (Slate Model)
% 4. [Increment] 計算應變增量 (指定樣本)
%
% 額外輸出：AMS 主軸特徵 (方向與數值)

clear; clc; close all;

%% ===== 1. 參數設定 =====
FILENAME = 'Debby_rawdata.xlsx'; 

% Slate 模型參數 (用於將磁感率轉為應變)
SLATE_A = 6.897; 
SLATE_B = 0.007;

% [步驟 4 設定] 指定要做增量分析的群組編號
ID_INITIAL = 1;  % 初始狀態 (Initial)
ID_FINAL   = 4;  % 最終狀態 (Final)

%% ===== 2. 讀取與前處理 =====
fprintf('正在讀取檔案...\n');
try
    raw_table = readtable(FILENAME, 'VariableNamingRule', 'preserve');
catch
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
tensor_db = containers.Map('KeyType', 'double', 'ValueType', 'any'); % 存 AMS 平均張量
result_count = 0;

fprintf('共 %d 個群組，開始依序執行步驟 1~3...\n', length(unique_ids));

%% ===== 3. 迴圈處理 (步驟 1, 2, 3 + 特徵提取) =====
for i = 1:length(unique_ids)
    curr_id = unique_ids(i);
    sub_table = raw_table(raw_table.(col_ID) == curr_id, :);
    num_subs = height(sub_table);
    
    tensor_stack_norm = zeros(3, 3, num_subs);
    
    % --- [步驟 1] 先轉座標系統 (Specimen -> Geo) ---
    for j = 1:num_subs
        row = sub_table(j, :);
        % 提取 K 值
        K1 = row.K1; K2 = row.K2; K3 = row.K3;
        % 提取方向 (並建立旋轉矩陣 V)
        v1 = [cosd(row.iK1geo)*cosd(row.dK1geo); cosd(row.iK1geo)*sind(row.dK1geo); sind(row.iK1geo)];
        v2 = [cosd(row.iK2geo)*cosd(row.dK2geo); cosd(row.iK2geo)*sind(row.dK2geo); sind(row.iK2geo)];
        v3 = [cosd(row.iK3geo)*cosd(row.dK3geo); cosd(row.iK3geo)*sind(row.dK3geo); sind(row.iK3geo)];
        V_specimen = [v1, v2, v3];
        
        % 計算地理座標下的磁感率張量 Kg
        Kg = AMSFormulas.calculateKg(K1, K2, K3, V_specimen);
        
        % --- [步驟 2-1] 做標準化 (Jelinek Normalization) ---
        Kg_norm = AMSFormulas.standardizeTensor(Kg);
        tensor_stack_norm(:, :, j) = Kg_norm;
    end
    
    % --- [步驟 2-2] 做平均 (Average) ---
    Mean_Kg = AMSFormulas.calculateMeanTensor(tensor_stack_norm);
    
    % 將平均張量存入資料庫 (供步驟 4 使用)
    tensor_db(curr_id) = Mean_Kg;
    
    % --- [額外計算] AMS 特徵提取 (Eigen-analysis) ---
    % 這裡提取的是「平均磁感率張量」的方向與主值
    ams_analysis = AMSFormulas.performEigenAnalysis(Mean_Kg);
    
    % --- [步驟 3] 轉應變 (Convert to Strain) ---
    % 利用 AMS 的特徵值 (K_vals) 透過 Slate 公式算出應變主軸 (Er)
    K_vals = ams_analysis.vals;
    Er_mean = AMSFormulas.calculateEr(K_vals(1), K_vals(2), K_vals(3), SLATE_A, SLATE_B);
    
    % 將應變主軸轉回地理座標 (Eg) - 使用 AMS 的方向矩陣
    Eg_mean = AMSFormulas.calculateEg(Er_mean, ams_analysis.V);
    
    % 對應變張量做特徵提取 (為了得到應變軸 X, Y, Z 的大小)
    strain_analysis = AMSFormulas.performEigenAnalysis(Eg_mean);
    
    
    %% --- 儲存結果 ---
    result_count = result_count + 1;
    results_struct(result_count).ID = curr_id;
    results_struct(result_count).N  = num_subs;
    
    % A. AMS 主軸特徵 (方向與數值)
    results_struct(result_count).AMS_K1_Val = ams_analysis.vals(1);
    results_struct(result_count).AMS_K1_Dec = ams_analysis.K1_dir(1);
    results_struct(result_count).AMS_K1_Inc = ams_analysis.K1_dir(2);
    
    results_struct(result_count).AMS_K2_Val = ams_analysis.vals(2);
    results_struct(result_count).AMS_K2_Dec = ams_analysis.K2_dir(1);
    results_struct(result_count).AMS_K2_Inc = ams_analysis.K2_dir(2);
    
    results_struct(result_count).AMS_K3_Val = ams_analysis.vals(3);
    results_struct(result_count).AMS_K3_Dec = ams_analysis.K3_dir(1);
    results_struct(result_count).AMS_K3_Inc = ams_analysis.K3_dir(2);
    
    % B. AMS 形狀因子
    results_struct(result_count).L = ams_analysis.L;
    results_struct(result_count).F = ams_analysis.F;
    results_struct(result_count).P = ams_analysis.P;
    
    % C. 應變軸數值 (Strain Magnitude)
    results_struct(result_count).Strain_X = strain_analysis.vals(1);
    results_struct(result_count).Strain_Y = strain_analysis.vals(2);
    results_struct(result_count).Strain_Z = strain_analysis.vals(3);
    
    fprintf('組別 %d 完成。 AMS K1: %.1f/%.1f, Strain X: %.3f\n', ...
        curr_id, ams_analysis.K1_dir(1), ams_analysis.K1_dir(2), strain_analysis.vals(1));
end

%% ===== 4. 輸出表格 =====
Final_Table = struct2table(results_struct);
disp(' ');
disp('==================== 群組分析結果 (AMS & Strain) ====================');
disp(Final_Table);
writetable(Final_Table, 'Step3_Strain_Results.xlsx');


%% ===== [步驟 4] 算應變增量 (Incremental Analysis) =====
disp(' ');
disp('============================================================');
fprintf('           [步驟 4] 應變增量分析: No.%d -> No.%d\n', ID_INITIAL, ID_FINAL);
disp('============================================================');

if isKey(tensor_db, ID_INITIAL) && isKey(tensor_db, ID_FINAL)
    % 取出平均 AMS 張量
    Kg_init = tensor_db(ID_INITIAL);
    Kg_final = tensor_db(ID_FINAL);
    
    % 計算增量張量 FK (Magnetic Fabric Change)
    % 公式: FK = Kg_final * inv(Kg_init)
    FK = AMSFormulas.calculateFK(Kg_init, Kg_final);
    
    % 分析 FK (特徵提取)
    FK_analysis = AMSFormulas.performEigenAnalysis(FK);
    
    disp('增量張量特徵值 (Incremental Eigenvalues):');
    disp(FK_analysis.vals');
    
    disp('增量張量主軸方向 (Incremental Directions):');
    fprintf('  Max Axis: Trend %.1f, Plunge %.1f\n', FK_analysis.K1_dir(1), FK_analysis.K1_dir(2));
    fprintf('  Min Axis: Trend %.1f, Plunge %.1f\n', FK_analysis.K3_dir(1), FK_analysis.K3_dir(2));
    
    % --- 繪圖 (3D 橢球) ---
    figure('Name', 'Incremental Analysis', 'Position', [100, 100, 1200, 400]);
    
    % 取得初始與最終的 AMS 特徵 (為了繪圖)
    An_Init = AMSFormulas.performEigenAnalysis(Kg_init);
    An_Final = AMSFormulas.performEigenAnalysis(Kg_final);
    
    subplot(1,3,1); 
    plotEllipsoidSub(An_Init.vals, An_Init.V, sprintf('Initial AMS (No.%d)', ID_INITIAL));
    
    subplot(1,3,2); 
    plotEllipsoidSub(FK_analysis.vals, FK_analysis.V, 'Incremental Tensor (FK)');
    
    subplot(1,3,3); 
    plotEllipsoidSub(An_Final.vals, An_Final.V, sprintf('Final AMS (No.%d)', ID_FINAL));
    
else
    error('指定的 ID 不存在，請檢查 ID_INITIAL 與 ID_FINAL 設定。');
end


%% ===== 繪圖輔助函式 =====
function plotEllipsoidSub(vals, V, titleText)
    [u,v] = meshgrid(linspace(0,2*pi,30), linspace(0,pi,15));
    x = vals(1)*cos(u).*sin(v); y = vals(2)*sin(u).*sin(v); z = vals(3)*cos(v);
    
    pts = V * [x(:)'; y(:)'; z(:)'];
    x_r = reshape(pts(1,:), size(x)); y_r = reshape(pts(2,:), size(y)); z_r = reshape(pts(3,:), size(z));
    
    surf(x_r, y_r, z_r, 'FaceAlpha',0.6, 'EdgeColor','none', 'FaceColor','cyan');
    hold on; axis equal; grid on; set(gca,'XDir','reverse');
    
    m = max(vals)*1.5;
    plot3([-m m],[0 0],[0 0],'k--'); plot3([0 0],[-m m],[0 0],'k--'); plot3([0 0],[0 0],[-m m],'k--');
    text(m,0,0,'N'); text(0,m,0,'E'); text(0,0,m,'D');
    
    quiver3(0,0,0, V(1,1)*vals(1)*1.2, V(2,1)*vals(1)*1.2, V(3,1)*vals(1)*1.2, 'r', 'LineWidth',2);
    quiver3(0,0,0, V(1,3)*vals(3)*1.2, V(2,3)*vals(3)*1.2, V(3,3)*vals(3)*1.2, 'b', 'LineWidth',2);
    title(titleText); view(135,30); camlight; lighting gouraud; hold off;
end