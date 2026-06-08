function results = read_RAN(filename)
% READ_RAN  讀取 Anisoft/AGICO 的 .RAN 二進位檔案並輸出 AMS 報表
%
% 用法:
%   results = read_RAN('RUEI2024.RAN')
%   results = read_RAN()          % 彈出檔案選擇視窗
%
% 輸出 (struct array，每個元素 = 一個樣品):
%   .name        樣品名稱
%   .k_bulk_SI   平均磁化率 (SI，無量綱)
%   .k1 k2 k3    Normalized 主磁化率 (k1>=k2>=k3, sum=3)
%   .T_norm      3×3 normalized susceptibility tensor
%   .eigvec      特徵向量矩陣 (各欄對應 k1, k2, k3 的方向)
%   .L           Lineation  = k1/k2
%   .F           Foliation  = k2/k3
%   .P           Degree of anisotropy = k1/k3
%   .Pprime      Jelinek P' = exp(sqrt(2*sum(eta_i^2)))
%   .T_shape     Shape parameter = (2lnk2-lnk1-lnk3)/(lnk1-lnk3)
%   .U           Lisle's U  = (2k2-k1-k3)/(k1-k3)
%   .dirs        [Dec Inc] × 3 軸方向 (specimen 座標系，deg)
%   .instr_code  儀器型號碼 (e.g., 'SST'=KLY-3, 'SSN'=KLY-4/5)
%
% 亦可呼叫 plot_RAN(results) 繪製 Flinn 圖等統計圖形。
%
% .RAN 二進位格式 (64 bytes / record, little-endian float32):
%   Header (bytes 1–64):  byte 1 = 測量型別; bytes 49–58 = O.P. 字串
%   每筆樣品 (64 bytes):
%     bytes  1–12  : 樣品名稱 (ASCII, zero/space padded)
%     bytes 13–16  : k_bulk  (float32) × 1e-6 = SI
%     bytes 17–20  : k_bulk  (重複)
%     bytes 21–24  : k11  (normalized tensor 對角線)
%     bytes 25–28  : k22
%     bytes 29–32  : k33
%     bytes 33–36  : k12  (off-diagonal; k11+k22+k33 = 3)
%     bytes 37–40  : k23
%     bytes 41–44  : k13
%     bytes 45–47  : 儀器代碼 (ASCII 3 chars)
%     bytes 48–64  : 保留

    if nargin < 1
        [file, path] = uigetfile('*.RAN', '選擇 .RAN 檔案');
        if isequal(file, 0), error('未選擇檔案'); end
        filename = fullfile(path, file);
    end

    %% ── 讀取二進位檔 ──────────────────────────────────────────────────
    fid = fopen(filename, 'rb', 'l');   % 'l' = little-endian
    if fid == -1, error('無法開啟檔案: %s', filename); end
    raw = fread(fid, inf, 'uint8=>uint8');
    fclose(fid);

    RECORD_SIZE = 64;
    num_recs    = floor((length(raw) - RECORD_SIZE) / RECORD_SIZE);

    % 解析 Header record
    op_str    = strtrim(char(raw(49:58)'));   % O.P. 字串
    meas_type = char(raw(1));                 % 測量型別碼

    %% ── 逐筆解析 ──────────────────────────────────────────────────────
    results = struct([]);
    count   = 0;

    for i = 1:num_recs
        base = RECORD_SIZE + (i-1)*RECORD_SIZE;        % 0-indexed
        rec  = raw(base+1 : base+RECORD_SIZE);         % 1-indexed slice

        % 樣品名稱 (bytes 1–12)
        name = strtrim(char(rec(1:12)'));
        if isempty(name) || all(rec(1:12) == 0), continue; end

        % 8 個 float32，bytes 13–44
        vals = double(typecast(rec(13:44), 'single'));
        k_bulk_raw = vals(1);                  % 儀器原始值
        k11 = vals(3);  k22 = vals(4);  k33 = vals(5);
        k12 = vals(6);  k23 = vals(7);  k13 = vals(8);

        instr_code = strtrim(char(rec(45:47)'));

        % 建立對稱 normalized tensor
        Tn = [k11 k12 k13;
              k12 k22 k23;
              k13 k23 k33];

        % 特徵值分解；由大到小排序 k1>=k2>=k3
        [V, D] = eig(Tn, 'vector');
        [ev, ix] = sort(D, 'descend');
        V = V(:, ix);
        k1 = ev(1);  k2 = ev(2);  k3 = ev(3);

        % AMS 參數 (Jelinek 1981)
        L       = k1 / k2;
        F       = k2 / k3;
        P       = k1 / k3;
        T_shape = (2*log(k2) - log(k1) - log(k3)) / (log(k1) - log(k3));
        U       = (2*k2 - k1 - k3) / (k1 - k3);
        eta     = log(ev) - mean(log(ev));
        Pprime  = exp(sqrt(2 * sum(eta.^2)));
        E       = L / F;

        % 主軸方向 → Dec / Inc (specimen system, z 朝下)
        dirs = zeros(3, 2);
        for j = 1:3
            v = V(:, j);
            if v(3) < 0, v = -v; end               % 統一到下半球
            dirs(j, 1) = mod(atan2d(v(2), v(1)), 360);
            dirs(j, 2) = asind(v(3));
        end

        % 存入
        count = count + 1;
        results(count).name        = name;
        results(count).k_bulk_raw  = k_bulk_raw;
        results(count).k_bulk_SI   = k_bulk_raw * 1e-6;
        results(count).T_norm      = Tn;
        results(count).eigvec      = V;
        results(count).k1          = k1;
        results(count).k2          = k2;
        results(count).k3          = k3;
        results(count).L           = L;
        results(count).F           = F;
        results(count).P           = P;
        results(count).Pprime      = Pprime;
        results(count).T_shape     = T_shape;
        results(count).U           = U;
        results(count).E           = E;
        results(count).dirs        = dirs;
        results(count).instr_code  = instr_code;
    end

    %% ── 印出報表 ──────────────────────────────────────────────────────
    print_ams_report(results, filename, meas_type, op_str);

end  % read_RAN


%% ═══════════════════════════════════════════════════════════════════════════
function print_ams_report(results, filename, meas_type, op_str)

    [~, fname, ext] = fileparts(filename);
    n   = length(results);
    S   = repmat('-', 1, 118);

    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('  %s  —  ANISOTROPY OF SUSCEPTIBILITY\n', [fname ext]);
    fprintf('  O.P.: %-20s  Type: %s   Samples: %d\n', op_str, meas_type, n);
    fprintf('%s\n\n', repmat('=', 1, 60));

    % 表 1: AMS 參數
    fprintf('%s\n', S);
    fprintf('%-12s  %10s  %7s  %7s  %7s  %6s  %6s  %7s  %7s  %7s  %6s\n', ...
        'Sample', 'k_bulk(SI)', 'k1', 'k2', 'k3', 'L', 'F', 'P', "P'", 'T', 'U');
    fprintf('%s\n', S);
    for i = 1:n
        r = results(i);
        fprintf('%-12s  %10.3e  %7.4f  %7.4f  %7.4f  %6.3f  %6.3f  %7.3f  %7.3f  %7.3f  %6.3f\n', ...
            r.name, r.k_bulk_SI, r.k1, r.k2, r.k3, ...
            r.L, r.F, r.P, r.Pprime, r.T_shape, r.U);
    end
    fprintf('%s\n\n', S);

    % 表 2: 主軸方向
    fprintf('%s\n', S);
    fprintf('%-12s   %8s %8s   %8s %8s   %8s %8s\n', ...
        'Sample', 'Dec_k1', 'Inc_k1', 'Dec_k2', 'Inc_k2', 'Dec_k3', 'Inc_k3');
    fprintf('%s\n', S);
    for i = 1:n
        r = results(i);
        fprintf('%-12s   %8.1f %8.1f   %8.1f %8.1f   %8.1f %8.1f\n', ...
            r.name, ...
            r.dirs(1,1), r.dirs(1,2), ...
            r.dirs(2,1), r.dirs(2,2), ...
            r.dirs(3,1), r.dirs(3,2));
    end
    fprintf('%s\n\n', S);

    % 統計摘要
    collect = @(f) [results.(f)];
    fprintf('  %-20s  %9s  %9s  %9s  %9s\n', '參數', 'Mean', 'Median', 'Min', 'Max');
    fprintf('  %s\n', repmat('-', 1, 65));
    prt = @(lbl, v) fprintf('  %-20s  %9.4g  %9.4g  %9.4g  %9.4g\n', ...
        lbl, mean(v), median(v), min(v), max(v));
    prt('k_bulk (SI)',       collect('k_bulk_SI'));
    prt('k1 (normalized)',   collect('k1'));
    prt('k3 (normalized)',   collect('k3'));
    prt('L  (Lineation)',    collect('L'));
    prt('F  (Foliation)',    collect('F'));
    prt('P  (Anisotropy)',   collect('P'));
    prt("P' (Jelinek)",      collect('Pprime'));
    prt('T  (Shape)',        collect('T_shape'));
    fprintf('\n');

end  % print_ams_report


%% ═══════════════════════════════════════════════════════════════════════════
function plot_RAN(results)
% PLOT_RAN  為 read_RAN 輸出繪製四種常用 AMS 統計圖
%
% 用法:
%   r = read_RAN('RUEI2024.RAN');
%   plot_RAN(r);

    collect = @(f) [results.(f)];
    n    = length(results);
    L    = collect('L');   F    = collect('F');
    P    = collect('P');   T    = collect('T_shape');
    Pp   = collect('Pprime');
    kb   = collect('k_bulk_SI');
    names = {results.name};

    figure('Name', 'AMS Summary Plots', 'NumberTitle', 'off', ...
           'Position', [80 80 1100 800]);

    %── 1. Flinn 圖 (L vs F) ────────────────────────────────────────────
    ax1 = subplot(2, 2, 1);
    scatter(L, F, 45, kb, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.4);
    cb = colorbar;  cb.Label.String = 'k_{bulk} (SI)';
    colormap(ax1, 'jet');
    hold on;
    lim = max([L F]) * 1.08;
    plot([1 lim], [1 lim], 'k--', 'LineWidth', 0.8);   % 對角線 (L=F)
    xlim([0.95 max(L)*1.05]);  ylim([0.95 max(F)*1.05]);
    xlabel('L = k_1/k_2  (Lineation)');
    ylabel('F = k_2/k_3  (Foliation)');
    title('Flinn 圖');
    grid minor;  hold off;

    %── 2. Jelinek 圖 (T vs P') ─────────────────────────────────────────
    ax2 = subplot(2, 2, 2);
    scatter(T, Pp, 45, kb, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.4);
    cb2 = colorbar;  cb2.Label.String = 'k_{bulk} (SI)';
    colormap(ax2, 'jet');
    hold on;
    xline(0, 'k--', 'LineWidth', 0.8);
    xlim([-1.1 1.1]);  ylim([1 max(Pp)*1.1]);
    xlabel('T  (Shape:  +1 = oblate,  −1 = prolate)');
    ylabel("P'  (Jelinek corrected P)");
    title("Jelinek 圖  (T  vs  P')");
    grid minor;  hold off;

    %── 3. k_bulk 直方圖 ────────────────────────────────────────────────
    subplot(2, 2, 3);
    histogram(log10(kb), 20, 'FaceColor', [0.25 0.55 0.80], 'EdgeColor', 'k');
    xlabel('log_{10}(k_{bulk})  [SI]');
    ylabel('Frequency');
    title('平均磁化率分布');
    grid on;

    %── 4. P vs k_bulk，顏色 = T ──────────────────────────────────────
    ax4 = subplot(2, 2, 4);
    scatter(kb, P, 45, T, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.4);
    cb4 = colorbar;  cb4.Label.String = 'T (Shape)';
    colormap(ax4, 'cool');
    caxis([-1 1]);
    hold on;
    set(gca, 'XScale', 'log');
    xlabel('k_{bulk}  (SI, log scale)');
    ylabel('P  (Degree of anisotropy)');
    title('P  vs  k_{bulk}  (顏色 = T)');
    grid on;

    % 標記 P 值異常大的樣品 (> mean + 2σ)
    thr = mean(P) + 2*std(P);
    idx = find(P > thr);
    if ~isempty(idx)
        scatter(ax4, kb(idx), P(idx), 80, 'r', 'x', 'LineWidth', 2);
        for k = idx(:)'
            text(ax4, kb(k), P(k) + 0.05, names{k}, ...
                'FontSize', 7, 'Color', 'r', 'HorizontalAlignment', 'center');
        end
    end
    hold off;

    sgtitle(sprintf('AMS 資料摘要  (n = %d)', n), ...
        'FontSize', 14, 'FontWeight', 'bold');

end  % plot_RAN