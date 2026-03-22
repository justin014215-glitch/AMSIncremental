%% 腳本：多樣本群組分析 (雙重統計版: Jelinek + Bootstrap)
%
% 執行流程：
% 1. [前處理] 轉座標 (Geo) -> Jelinek 標準化 -> 平均
% 2. [統計] 同時計算 Jelinek (信心橢圓) 與 Bootstrap (信心角/誤差)
% 3. [應變] 轉換為應變張量 (Slate)
% 4. [增量] 計算指定群組間的應變增量 (FK)
%
% 輸出：包含 AMS 主軸、Jelinek 參數、Bootstrap 誤差、應變軸的完整總表

clear; clc; close all;

%% ===== 1. 參數設定 =====
FILENAME = 'test0206.xlsx';
SLATE_A = 6.897; 
SLATE_B = 0.007;

% 統計設定
NUM_BOOTSTRAPS = 1000; % Bootstrap 重抽樣次數

% 增量分析設定
ID_INITIAL = 1;
ID_FINAL   = 2;

%% ===== 2. 讀取與前處理 =====
fprintf('正在讀取檔案...\n');
try
    raw_table = readtable(FILENAME, 'VariableNamingRule', 'preserve');
catch
    raw_table = readtable(FILENAME); 
end

% 自動偵測樣本編號
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

results_struct = struct();
tensor_db = containers.Map('KeyType', 'double', 'ValueType', 'any');
result_count = 0;

fprintf('共 %d 個群組，開始分析...\n', length(unique_ids));

%% ===== 3. 迴圈處理 (計算與統計) =====
for i = 1:length(unique_ids)
    curr_id = unique_ids(i);
    sub_table = raw_table(raw_table.(col_ID) == curr_id, :);
    num_subs = height(sub_table);
    
    tensor_stack_norm = zeros(3, 3, num_subs);
    
    % --- Step 1: 轉座標 & 標準化 ---
    for j = 1:num_subs
        row = sub_table(j, :);
        K1 = row.K1; K2 = row.K2; K3 = row.K3;
        v1 = [cosd(row.iK1geo)*cosd(row.dK1geo); cosd(row.iK1geo)*sind(row.dK1geo); sind(row.iK1geo)];
        v2 = [cosd(row.iK2geo)*cosd(row.dK2geo); cosd(row.iK2geo)*sind(row.dK2geo); sind(row.iK2geo)];
        v3 = [cosd(row.iK3geo)*cosd(row.dK3geo); cosd(row.iK3geo)*sind(row.dK3geo); sind(row.iK3geo)];
        V = [v1, v2, v3];
        
        Kg = AMSFormulas.calculateKg(K1, K2, K3, V);
        Kg_norm = AMSFormulas.standardizeTensor(Kg);
        tensor_stack_norm(:, :, j) = Kg_norm;
    end
    
    % --- Step 2: 平均與 AMS 特徵提取 ---
    Mean_Kg = AMSFormulas.calculateMeanTensor(tensor_stack_norm);
    tensor_db(curr_id) = Mean_Kg;
    ams = AMSFormulas.performEigenAnalysis(Mean_Kg);
    
    % --- Step 3: 雙重統計 (Dual Statistics) ---
    % A. Jelinek Statistics
    jel_stats = AMSFormulas.calculateJelinekStats(tensor_stack_norm);
    
    % B. Bootstrap Statistics
    if num_subs >= 3
        boot_stats = AMSFormulas.calculateBootstrapStats(tensor_stack_norm, NUM_BOOTSTRAPS);
    else
        boot_stats = struct('K1_std',NaN, 'K3_std',NaN, 'V1_conf',NaN, 'V3_conf',NaN);
    end
    
    % --- Step 4: 轉應變 (Strain) ---
    K_vals = ams.vals;
    Er_mean = AMSFormulas.calculateEr(K_vals(1), K_vals(2), K_vals(3), SLATE_A, SLATE_B);
    Eg_mean = AMSFormulas.calculateEg(Er_mean, ams.V);
    strain = AMSFormulas.performEigenAnalysis(Eg_mean);
    
    % --- Step 5: 儲存結果 ---
    result_count = result_count + 1;
    results_struct(result_count).ID = curr_id;
    results_struct(result_count).N  = num_subs;
    
    % 1. AMS 數值 (K)
    results_struct(result_count).K1 = ams.vals(1);
    results_struct(result_count).K2 = ams.vals(2);
    results_struct(result_count).K3 = ams.vals(3);
    
    % 2. 數值誤差 (Bootstrap Std Dev vs Jelinek Std Err)
    results_struct(result_count).K1_Err_Boot = boot_stats.K1_std;
    results_struct(result_count).K1_Err_Jel  = jel_stats.K1_err;
    results_struct(result_count).K3_Err_Boot = boot_stats.K3_std;
    results_struct(result_count).K3_Err_Jel  = jel_stats.K3_err;
    
    % 3. AMS 方向 (Trend/Plunge)
    results_struct(result_count).K1_Dec = ams.K1_dir(1); results_struct(result_count).K1_Inc = ams.K1_dir(2);
    results_struct(result_count).K3_Dec = ams.K3_dir(1); results_struct(result_count).K3_Inc = ams.K3_dir(2);
    
    % 4. 方向誤差 (Bootstrap Angle vs Jelinek Ellipses)
    % Bootstrap: 單一數值 (圓錐半徑)
    results_struct(result_count).K1_Conf_Boot = boot_stats.V1_conf;
    results_struct(result_count).K3_Conf_Boot = boot_stats.V3_conf;
    
    % Jelinek: 橢圓半軸 (E12, E13)
    results_struct(result_count).E12_Jel = jel_stats.E12; % K1 在 1-2 面上的誤差
    results_struct(result_count).E13_Jel = jel_stats.E13; % K1 在 1-3 面上的誤差
    results_struct(result_count).E23_Jel = jel_stats.E23; % K2/K3 在 2-3 面上的誤差
    
    % 5. 應變與形狀
    results_struct(result_count).L = ams.L;
    results_struct(result_count).F = ams.F;
    results_struct(result_count).Strain_X = strain.vals(1);
    results_struct(result_count).Strain_Z = strain.vals(3);
    
    fprintf('組別 %d 完成。 Bootstrap K1誤差: %.4f | Jelinek E12: %.1f\n', ...
        curr_id, boot_stats.K1_std, jel_stats.E12);
end

%% ===== 4. 輸出結果 =====
Final_Table = struct2table(results_struct);
disp(' ');
disp('==================== 綜合分析結果 (部分預覽) ====================');
% 挑選最具代表性的欄位顯示
preview_cols = {'ID', 'N', 'K1', 'K1_Err_Boot', 'K1_Conf_Boot', 'E12_Jel', 'E13_Jel', 'Strain_X'};
if all(ismember(preview_cols, Final_Table.Properties.VariableNames))
    disp(Final_Table(:, preview_cols));
else
    disp(Final_Table);
end
writetable(Final_Table, 'Result_DualStats.xlsx');
fprintf('\n完整結果已儲存至 Result_DualStats.xlsx\n');


%% ===== 5. 增量分析 (與之前相同) =====
disp(' ');
disp('---------------- 增量分析 (Incremental) ----------------');
if isKey(tensor_db, ID_INITIAL) && isKey(tensor_db, ID_FINAL)
    Kg_init = tensor_db(ID_INITIAL);
    Kg_final = tensor_db(ID_FINAL);
    FK = AMSFormulas.calculateFK(Kg_init, Kg_final);
    FK_res = AMSFormulas.performEigenAnalysis(FK);
    
    fprintf('增量 (No.%d -> No.%d) 主軸方向:\n', ID_INITIAL, ID_FINAL);
    fprintf('Max: %.1f / %.1f\n', FK_res.K1_dir(1), FK_res.K1_dir(2));
    
    % 繪圖部分 (如需繪圖請保留 v6 的繪圖程式碼，呼叫 plotEllipsoidSub)
    figure('Name', 'Incremental Analysis', 'Position', [100, 100, 1000, 350]);
    An_Init = AMSFormulas.performEigenAnalysis(Kg_init);
    An_Final = AMSFormulas.performEigenAnalysis(Kg_final);
    
    subplot(1,3,1); plotEllipsoidSub(An_Init.vals, An_Init.V, sprintf('Initial (No.%d)', ID_INITIAL));
    subplot(1,3,2); plotEllipsoidSub(FK_res.vals, FK_res.V, 'Incremental FK');
    subplot(1,3,3); plotEllipsoidSub(An_Final.vals, An_Final.V, sprintf('Final (No.%d)', ID_FINAL));
else
    warning('指定的增量 ID 不存在。');
end

function plotEllipsoidSub(vals, V, titleText)
    [u,v] = meshgrid(linspace(0,2*pi,30), linspace(0,pi,15));
    x = vals(1)*cos(u).*sin(v); y = vals(2)*sin(u).*sin(v); z = vals(3)*cos(v);
    pts = V * [x(:)'; y(:)'; z(:)'];
    x_r = reshape(pts(1,:), size(x)); y_r = reshape(pts(2,:), size(y)); z_r = reshape(pts(3,:), size(z));
    surf(x_r, y_r, z_r, 'FaceAlpha',0.6, 'EdgeColor','none', 'FaceColor','cyan');
    hold on; axis equal; grid on; set(gca,'XDir','reverse');
    m = max(vals)*1.2;
    plot3([-m m],[0 0],[0 0],'k--'); plot3([0 0],[-m m],[0 0],'k--'); plot3([0 0],[0 0],[-m m],'k--');
    quiver3(0,0,0, V(1,1)*vals(1)*1.1, V(2,1)*vals(1)*1.1, V(3,1)*vals(1)*1.1, 'r', 'LineWidth',2);
    quiver3(0,0,0, V(1,3)*vals(3)*1.1, V(2,3)*vals(3)*1.1, V(3,3)*vals(3)*1.1, 'b', 'LineWidth',2);
    title(titleText); view(135,30); camlight; lighting gouraud; hold off;
end