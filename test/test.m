% 讀取 Excel 檔案，保留原始欄位名稱
filename = 'NCIH.xlsx';
data = readtable(filename, 'VariableNamingRule', 'preserve');

% 顯示 Excel 檔案的欄位名稱
disp('欄位名稱：');
disp(data.Properties.VariableNames);

% 假設欄位名稱為 K1, K2, K3, Pj，直接讀取這些欄位
K1 = data.K1;   % 最大主軸磁感率
K2 = data.K2;   % 中間主軸磁感率
K3 = data.K3;   % 最小主軸磁感率
Km = data.Km;   % 平均磁感率
F = data.F;
L = data.L;
Int = data.Int;
T = data.T;
Pj = data.Pj;   % 校正異向性 (Pj)
K0 = (K1 .* K2 .* K3).^(1/3);
data.K0 = K0;
a = 6.897;
b = 0.007;

% 顯示已讀取的數據（選擇性）
disp('讀取的數據：');
disp(data);



