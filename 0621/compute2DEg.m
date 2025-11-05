function Eg2D = compute2DEg(K1, K2, dK1geo, iK1geo, dK2geo, iK2geo, A, B)
    % 輸入：
    % K1, K2 - 主磁化率
    % dK1geo, iK1geo - K1的方向（度）
    % dK2geo, iK2geo - K2的方向（度）
    % A, B - Slate Coefficients (例如 6.897, 0.007)
    [filename, pathname] = uigetfile('*.xlsx', '請選擇 AMS 數據檔');
        if isequal(filename, 0)
            disp('使用者取消選擇。');
            return;
        end
    % 先把角度轉成弧度
    dK1 = deg2rad(dK1geo);
    iK1 = deg2rad(iK1geo);
    dK2 = deg2rad(dK2geo);
    iK2 = deg2rad(iK2geo);

    % 計算K0幾何平均
    K0 = (K1 * K2)^(1/2);

    % 計算有限應變指標e1, e2
    ln1pe1 = A * (K1/K0 - 1) - B;
    ln1pe2 = A * (K2/K0 - 1) - B;
    e1 = exp(ln1pe1) - 1;
    e2 = exp(ln1pe2) - 1;
    omega = (1 + e1)*(1 + e2);

    % 計算Er (2x2對角矩陣)
    Er = omega^2 * diag([(1+e1)^(-2), (1+e2)^(-2)]);

    % 計算方向矩陣V (2x2)
    % 這裡只考慮水平平面（X-Y平面）投影，忽略垂直方向
    % 以方位角與俯仰角轉換成方向向量並取前兩維
    v1 = [cos(iK1)*cos(dK1); cos(iK1)*sin(dK1)];
    v2 = [cos(iK2)*cos(dK2); cos(iK2)*sin(dK2)];
    V = [v1, v2];

    % 計算Eg
    Eg2D = V * Er * V';

end
