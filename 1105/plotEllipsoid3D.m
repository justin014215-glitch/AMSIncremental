function plotEllipsoid3D(a, b, c, V, sampleName)
% plotEllipsoid3D(a, b, c, V, sampleName)
% 繪製 3D 橢球並顯示加強版的座標軸
%
% 修改說明：已設定 X 軸方向反轉 (XDir = reverse)
%
% 輸入:
%   a, b, c    : 三軸長度 (K1, K2, K3)
%   V          : 3x3 特徵向量矩陣 (旋轉矩陣)
%   sampleName : 樣本名稱 (選填)

    if nargin < 4 || isempty(V)
        V = eye(3);  % 預設為未旋轉
    end
    if nargin < 5
        sampleName = '';
    end
    
    % 1. 建立橢球表面數據
    [u, v] = meshgrid(linspace(0, 2*pi, 60), linspace(0, pi, 30));
    x = a * cos(u) .* sin(v);
    y = b * sin(u) .* sin(v);
    z = c * cos(v);
    
    % 2. 旋轉點雲
    pts = V * [x(:)'; y(:)'; z(:)'];
    x_rot = reshape(pts(1, :), size(x));
    y_rot = reshape(pts(2, :), size(y));
    z_rot = reshape(pts(3, :), size(z));
    
    % 3. 計算繪圖範圍 (讓軸長度大於橢球最大半徑)
    max_radius = max([a, b, c]);
    axis_len = max_radius * 1.5; % 設定軸長度為最大半徑的 1.5 倍
    
    % 4. 開始繪圖
    figure;
    hold on;
    grid on;
    axis equal;
    
    % --- [修改點] 反轉 X 軸方向 ---
    set(gca, 'XDir', 'reverse'); 
    
    % --- (A) 繪製更明顯的世界座標 XYZ 軸 (黑色長線) ---
    % 畫出穿過原點的 X, Y, Z 軸
    plot3([-axis_len axis_len], [0 0], [0 0], 'k--', 'LineWidth', 1); % X 軸
    plot3([0 0], [-axis_len axis_len], [0 0], 'k--', 'LineWidth', 1); % Y 軸
    plot3([0 0], [0 0], [-axis_len axis_len], 'k--', 'LineWidth', 1); % Z 軸
    
    % 在軸的末端加上文字標籤
    % 注意：因為 X 軸反轉了，(axis_len, 0, 0) 會顯示在圖的另一邊，但數值邏輯不變
    text(axis_len, 0, 0, ' X', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
    text(0, axis_len, 0, ' Y', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
    text(0, 0, axis_len, ' Z', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
    
    % --- (B) 繪製橢球表面 ---
    surf(x_rot, y_rot, z_rot, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'cyan');
    
    % --- (C) 繪製主軸方向 (紅/綠/藍箭頭) ---
    % 設定放大係數，讓箭頭稍微突出橢球表面
    arrow_scale = 1.2; 
    
    % 使用 Quiver3 繪製向量，並關閉 AutoScale 以便精確控制長度
    q1 = quiver3(0, 0, 0, V(1,1)*a*arrow_scale, V(2,1)*a*arrow_scale, V(3,1)*a*arrow_scale, ...
        'r', 'LineWidth', 3, 'AutoScale', 'off', 'MaxHeadSize', 0.5);
        
    q2 = quiver3(0, 0, 0, V(1,2)*b*arrow_scale, V(2,2)*b*arrow_scale, V(3,2)*b*arrow_scale, ...
        'g', 'LineWidth', 3, 'AutoScale', 'off', 'MaxHeadSize', 0.5);
        
    q3 = quiver3(0, 0, 0, V(1,3)*c*arrow_scale, V(2,3)*c*arrow_scale, V(3,3)*c*arrow_scale, ...
        'b', 'LineWidth', 3, 'AutoScale', 'off', 'MaxHeadSize', 0.5);
    
    % --- (D) 視覺美化設定 ---
    title(['3D Ellipsoid: ', sampleName], 'FontSize', 12);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    
    % 設定視圖範圍
    xlim([-axis_len axis_len]);
    ylim([-axis_len axis_len]);
    zlim([-axis_len axis_len]);
    
    % 燈光與視角
    camlight; 
    lighting gouraud;
    view(120, 30); % 設定初始視角
    
    legend([q1, q2, q3], {'K1 (Max)', 'K2 (Int)', 'K3 (Min)'}, 'Location', 'best');
    hold off;
end