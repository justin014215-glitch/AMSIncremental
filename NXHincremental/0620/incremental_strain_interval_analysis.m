function [Einc, U, R, eigvals_U, eigvecs_U, Eg_avg, sample_ranges] = ...
    incremental_strain_interval_analysis(Eg, interval_count, samples_per_interval)
% 執行 Eg 區間平均 → 應變增量計算 → polar decomposition

    % Step 1: Eg 區間平均
    [Eg_avg, sample_ranges] = calculate_average_Eg(Eg, interval_count, samples_per_interval);

    % 預先配置空間
    Einc = zeros(3,3,interval_count-1);
    U = zeros(3,3,interval_count-1);
    R = zeros(3,3,interval_count-1);
    eigvals_U = zeros(3, interval_count-1);
    eigvecs_U = zeros(3, 3, interval_count-1);

    for i = 1:interval_count-1
        E1 = Eg_avg(:,:,i);
        E2 = Eg_avg(:,:,i+1);
        
        % Step 2: 應變增量張量 F = E2 * inv(E1)
        F = E2 * inv(E1);
        Einc(:,:,i) = F;

        % Step 3: polar decomposition → F = R * U
        C = F' * F;
        U(:,:,i) = sqrtm(C);                   % 對稱 stretch tensor
        R(:,:,i) = F * inv(U(:,:,i));          % 旋轉張量

        % Step 4: 主值主軸分解（形變率張量）
        [eigvec, eigval] = eig(U(:,:,i));
        [sorted_vals, idx] = sort(diag(eigval), 'descend');
        eigvecs_U(:,:,i) = eigvec(:,idx);
        eigvals_U(:,i) = sorted_vals;
    end
end
