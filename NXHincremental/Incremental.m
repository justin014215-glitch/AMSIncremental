
clear all
close all
clc

%%
% 讀取 Excel 檔案，保留原始欄位名稱
filename = 'NXHAMStest.xlsx';
data = readtable(filename, 'VariableNamingRule', 'preserve');

%%
%{
顯示 Excel 檔案的欄位名稱
disp('欄位名稱：');
disp(data.Properties.VariableNames);
%}

% 假設欄位名稱為 K1, K2, K3, Pj，直接讀取這些欄位
K1 = data.K1;   % 最大主軸磁感率
K2 = data.K2;   % 中間主軸磁感率
K3 = data.K3;   % 最小主軸磁感率
Km = data.Km;   % 平均磁感率
F = data.F;
L = data.L;
T = data.T;
P = data.P;
Pj = data.Pj;   % 校正異向性 (Pj)
K0 = (K1 .* K2 .* K3).^(1/3);
Int = (sqrt((L-1).^2+(F-1).^2));
data.K0 = K0;
data.Int = Int;
a = 6.897;
b = 0.007;

%顯示已讀取的數據（選擇性）
disp('讀取的數據：');
disp(data);

%%
clc
dK1 = data.dK1geo;  % K1 方位角
iK1 = data.iK1geo;  % K1 傾角
dK2 = data.dK2geo;  % K2 方位角
iK2 = data.iK2geo;  % K2 傾角
dK3 = data.dK3geo;  % K3 方位角
iK3 = data.iK3geo;  % K3 傾角

% 轉換角度為弧度
dK1 = deg2rad(dK1);
iK1 = deg2rad(iK1);
dK2 = deg2rad(dK2);
iK2 = deg2rad(iK2);
dK3 = deg2rad(dK3);
iK3 = deg2rad(iK3);
num_samples = length(K1);

V  = zeros(3,3,num_samples);

l1_matrix = zeros(num_samples, 1);
m1_matrix = zeros(num_samples, 1);
n1_matrix = zeros(num_samples, 1);
l2_matrix = zeros(num_samples, 1);
m2_matrix = zeros(num_samples, 1);
n2_matrix = zeros(num_samples, 1);
l3_matrix = zeros(num_samples, 1);
m3_matrix = zeros(num_samples, 1);
n3_matrix = zeros(num_samples, 1);

for i = 1:num_samples
    % 計算方向餘弦(特徵向量) (l, m, n)
    l1 = cos(iK1(i)) .* cos(dK1(i));
    m1 = cos(iK1(i)) .* sin(dK1(i));
    n1 = sin(iK1(i));

    l2 = cos(iK2(i)) .* cos(dK2(i));
    m2 = cos(iK2(i)) .* sin(dK2(i));
    n2 = sin(iK2(i));

    l3 = cos(iK3(i)) .* cos(dK3(i));
    m3 = cos(iK3(i)) .* sin(dK3(i));
    n3 = sin(iK3(i));


% 建立特徵向量矩陣 V (每一組樣本對應一個 3x3 矩陣)


    l1_matrix(i) = l1;
    m1_matrix(i) = m1;
    n1_matrix(i) = n1;
    l2_matrix(i) = l2;
    m2_matrix(i) = m2;
    n2_matrix(i) = n2;
    l3_matrix(i) = l3;
    m3_matrix(i) = m3;
    n3_matrix(i) = n3;
    V(:,:,i) = [l1_matrix(i), l2_matrix(i), l3_matrix(i);
                m1_matrix(i), m2_matrix(i), m3_matrix(i); 
                n1_matrix(i), n2_matrix(i), n3_matrix(i)];

end 


disp('前 5 組樣本的旋轉矩陣：');
for i = 1:min(5, num_samples)
    fprintf('樣本 %d:\n', i);
    disp(V(:,:,i));
    
end





%%

% 計算每個主軸對應的 ln(1+e)
ln1pe1 = a.*((K1./K0)-1)-b;
ln1pe2 = a.*((K2./K0)-1)-b;
ln1pe3 = a.*((K3./K0)-1)-b;


disp('前五組樣本ln1pe1')
for i = 1:min(5, num_samples)
    fprintf('樣本 %d:\n', i);
    disp(ln1pe1);
end

%%

% 求出有限應變 e_i
e1 = exp(ln1pe1) - 1;
e2 = exp(ln1pe2) - 1;
e3 = exp(ln1pe3) - 1;


%計算w
w = (1+e1).*(1+e2).*(1+e3);


%計算偽應變橢球Er
Er = zeros(3,3,num_samples);  % 初始化 3×3×N 矩陣 (N 為樣本數)

E11_matrix = zeros(num_samples, 1);
E22_matrix = zeros(num_samples, 1);
E33_matrix = zeros(num_samples, 1);

%%
for i = 1:num_samples
%{    
    E11 = ln1pe1(i)^(-2);  % 取得單一樣本的 E11
    E22 = ln1pe2(i)^(-2);  % 取得單一樣本的 E22
    E33 = ln1pe3(i)^(-2);  % 取得單一樣本的 E33
%}
    E11 =(ln1pe2(i)^2 )*(ln1pe3(i)^2);  % 取得單一樣本的 E11
    E22 = (ln1pe1(i)^2 )*(ln1pe3(i)^2);  % 取得單一樣本的 E22
    E33 = (ln1pe1(i)^2 )*(ln1pe2(i)^2);  % 取得單一樣本的 E33

    E11_matrix(i) = E11;  % 取得單一樣本的 E11
    E22_matrix(i) = E22;  % 取得單一樣本的 E22
    E33_matrix(i) = E33;  % 取得單一樣本的 E33

    % 計算 w^2
    w_squared = w(i)^2;  % 計算 w 的平方

    % 將結果儲存到 Er(:,:,i)
    Er(:,:,i) = w_squared * [E11, 0, 0;
                             0, E22, 0;
                             0, 0, E33];
end

disp('w^2')
disp(w_squared)

disp('前五組樣本E11');
for i = 1:min(5, num_samples)
    fprintf('樣本 %d:\n', i);
    disp(E11_matrix);
end

disp('前五組樣本E22');
for i = 1:min(5, num_samples)
    fprintf('樣本 %d:\n', i);
    disp(E22_matrix);
end

disp('前五組樣本E33');
for i = 1:min(5, num_samples)
    fprintf('樣本 %d:\n', i);
    disp(E33_matrix);
end



disp('前 5 組樣本的 Er 矩陣：');
for i = 1:min(5, num_samples)
    fprintf('樣本 %d:\n', i);
    disp(Er(:,:,i));
end

%%
clc

% 初始化 Eg 矩陣，與 Er 結構相同 Eg原位有限應變橢球
Eg = zeros(3, 3, num_samples);

for i = 1:num_samples
    % 取得對應樣本的旋轉矩陣 V 和應變矩陣 Er
    Vi = V(:,:,i);  % 取出第 i 個旋轉矩陣
    Ei = Er(:,:,i); % 取出第 i 個應變張量
    VT = Vi';
    % 計算 E_g = V^T * E_r * V
    Eg(:,:,i) = (VT) * Ei * (Vi);
end



%%
% 顯示前 5 組 Eg 矩陣
disp('前 5 組樣本的 Eg 矩陣：');
for i = 1:min(5, num_samples)
    fprintf('樣本 %d:\n', i);
    disp(Eg(:,:,i));
end
%%
clc
for  i = 1:3

Eg1 = Eg(:,:,i); % 第一顆樣本 (初始應變)
Eg3 = Eg(:,:,i+1); % 第二顆樣本 (變形後應變)
Eg2 = Eg3 * (Eg1)'; 

end

disp(Eg2)
%%
clc
L = Eg2;
U = zeros(3,3,num_samples);
R = zeros(3,3,num_samples);


% 計算 D_ij 和 W_ij
for i = 1
    Lij = L(:,:,i);  % 取出第 i 組 L 矩陣
    C = Lij' * Lij;
    % 計算對稱部分 D_ij 和反對稱部分 W_ij
    U(:,:,i) = C ^ 0.5 ;  
    R(:,:,i) = Lij * U(:,:,i) ^(-1) ;
end


disp('對稱部分 U：');
for i = 1
    fprintf('第一段 %d:\n', i);
    disp(U(:,:,i));
end 


disp('反對稱部分 R：');
for i = 1
    fprintf('第一段 %d:\n', i);
    disp(R(:,:,i));
end
%%
clc
Eb = VT.^(-1) .* U .* V.^(-1);

disp(Eb)