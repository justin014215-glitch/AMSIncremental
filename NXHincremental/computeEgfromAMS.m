clc
function Eg = compute_Eg_from_AMS(filename)
    % 計算磁感率資料轉換後的 Eg（有限應變橢球）
   

    data = readtable(filename, 'VariableNamingRule', 'preserve');

    K1 = data.K1; K2 = data.K2; K3 = data.K3;
    dK1 = deg2rad(data.dK1geo); iK1 = deg2rad(data.iK1geo);
    dK2 = deg2rad(data.dK2geo); iK2 = deg2rad(data.iK2geo);
    dK3 = deg2rad(data.dK3geo); iK3 = deg2rad(data.iK3geo);

    a = 6.897; b = 0.007;
    num_samples = length(K1);
    K0 = (K1 .* K2 .* K3).^(1/3);
    
    ln1pe1 = a .* ((K1 ./ K0) - 1) - b;
    ln1pe2 = a .* ((K2 ./ K0) - 1) - b;
    ln1pe3 = a .* ((K3 ./ K0) - 1) - b;

    e1 = exp(ln1pe1) - 1;
    e2 = exp(ln1pe2) - 1;
    e3 = exp(ln1pe3) - 1;
    omega = (1 + e1) .* (1 + e2) .* (1 + e3);

    Er = zeros(3, 3, num_samples);
    for i = 1:num_samples
        Er(:,:,i) = omega(i)^2 * diag([(1+e1(i))^(-2), (1+e2(i))^(-2), (1+e3(i))^(-2)]);
    end

    V = zeros(3,3,num_samples);
    for i = 1:num_samples
        V(:,:,i) = [cos(iK1(i))*cos(dK1(i)), cos(iK2(i))*cos(dK2(i)), cos(iK3(i))*cos(dK3(i));
                    cos(iK1(i))*sin(dK1(i)), cos(iK2(i))*sin(dK2(i)), cos(iK3(i))*sin(dK3(i));
                    sin(iK1(i)),             sin(iK2(i)),             sin(iK3(i))];
    end

    Eg = zeros(3,3,num_samples);
    for i = 1:num_samples
        Eg(:,:,i) = V(:,:,i)' * Er(:,:,i) * V(:,:,i);
    end
end
