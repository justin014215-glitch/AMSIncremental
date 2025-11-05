
%{
% 使用方法

K1 = 1; K2 = 1.2; K3 = 0.8;% 主軸長度
[V, ~] = eig([K1 0 0; 0 K2 0; 0 0 K3]);% 方向矩陣（可從 eig 計算）
plotEllipsoid3D(K1, K2, K3, V, 'Sample NXH25');% 呼叫繪圖函數

%}



function plotEllipsoid3D(a, b, c, V, sampleName)
% plotEllipsoid3D(a, b, c, V, sampleName)
% - a, b, c: 三軸長（可為應變或磁化率主軸）
% - V: 3x3 的特徵向量矩陣（optional，預設為單位軸）
% - sampleName: 樣本名稱字串（optional，用於圖上顯示）

    if nargin < 4 || isempty(V)
        V = eye(3);  % 預設為未旋轉橢球
    end
    if nargin < 5
        sampleName = '';
    end

    % 建立橢球表面參數
    [u, v] = meshgrid(linspace(0, 2*pi, 60), linspace(0, pi, 30));
    x = a * cos(u) .* sin(v);
    y = b * sin(u) .* sin(v);
    z = c * cos(v);

    % 展平成點雲並旋轉
    pts = V * [x(:)'; y(:)'; z(:)'];
    x_rot = reshape(pts(1, :), size(x));
    y_rot = reshape(pts(2, :), size(y));
    z_rot = reshape(pts(3, :), size(z));

    % 繪圖
    figure;
    surf(x_rot, y_rot, z_rot, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(['3D Ellipsoid: ', sampleName]);
    camlight; lighting gouraud;
    
    % 繪製主軸方向線
    hold on;
    quiver3(0, 0, 0, V(1,1)*a, V(2,1)*a, V(3,1)*a, 'r', 'LineWidth', 2);
    quiver3(0, 0, 0, V(1,2)*b, V(2,2)*b, V(3,2)*b, 'g', 'LineWidth', 2);
    quiver3(0, 0, 0, V(1,3)*c, V(2,3)*c, V(3,3)*c, 'b', 'LineWidth', 2);
    legend('Ellipsoid','Axis 1','Axis 2','Axis 3');

end
