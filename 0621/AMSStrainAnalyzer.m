classdef AMSStrainAnalyzer < handle
    properties
        filename
        data
        Eg
        ErList
        VList
        config
        Er_raw      % 原始磁化率構成的張量 Er_raw = diag(K1, K2, K3)
        EgRaw       % 由 Er_raw 推得的 Eg_raw = V' * Er_raw * V
        
        % 新增2D結果
        EgD
        ErList2D
        VList2D
    end

    methods
        
        %輸入檔案，單獨使用方法：
        %[filename, pathname] = uigetfile({'*.xlsx;*.csv;*.txt'}, '請選擇 AMS 數據檔');
        
        function obj = AMSStrainAnalyzer(filename, varargin)
            obj.filename = filename;
            obj.config = obj.parseConfig(varargin{:});
            obj.loadData();
        end

        %輸入磁感率轉應變之參數A,B
        
        function config = parseConfig(~, varargin)
            p = inputParser;
            addParameter(p, 'SlateCoeffA', 6.897, @isnumeric);
            addParameter(p, 'SlateCoeffB', 0.007, @isnumeric);
            addParameter(p, 'Verbose', true, @islogical);
            parse(p, varargin{:});
            config = p.Results;
        end

        %讀取資料
        function loadData(obj)
            try
                obj.data = readtable(obj.filename, 'VariableNamingRule', 'preserve');
                if obj.config.Verbose
                    fprintf('成功載入 %d 筆 AMS 數據\n', height(obj.data));
                end
            catch ME
                error('載入數據失敗：%s', ME.message);
            end

            required = {'K1','K2','K3','dK1geo','iK1geo','dK2geo','iK2geo','dK3geo','iK3geo'};
            missing = setdiff(required, obj.data.Properties.VariableNames);
            if ~isempty(missing)
                error('缺少必要欄位：%s', strjoin(missing, ', '));
            end
        end

        %設定參數(K1, K2, K3)
        %計算Er, Eg(計算程式在下兩個function)


        function computeFiniteStrainTensors(obj)
            K1 = obj.data.K1; K2 = obj.data.K2; K3 = obj.data.K3;
            dK1 = deg2rad(obj.data.dK1geo); iK1 = deg2rad(obj.data.iK1geo);
            dK2 = deg2rad(obj.data.dK2geo); iK2 = deg2rad(obj.data.iK2geo);
            dK3 = deg2rad(obj.data.dK3geo); iK3 = deg2rad(obj.data.iK3geo);

            n = length(K1);
            obj.Eg = zeros(3,3,n);
            obj.ErList = zeros(3,3,n);
            obj.VList = zeros(3,3,n);

            for i = 1:n
                Er = obj.computeEr(K1(i), K2(i), K3(i));
                V = obj.computeV(dK1(i), iK1(i), dK2(i), iK2(i), dK3(i), iK3(i));
                Eg = obj.computeEg(Er, V);
                obj.ErList(:,:,i) = Er;
                obj.VList(:,:,i) = V;
                obj.Eg(:,:,i) = Eg;
            end
        end

        
        %計算Er
        % computeEr - 計算磁化率應變張量
        %
        % 輸入參數：
        %   K1, K2, K3 - 主磁化率值
        %
        % 輸出參數：
        %   Er - 3x3 應變張量
        %
        % 計算方法：
        %   基於 Slate 係數的指數變換
        %   Er = ω² * diag([(1+e₁)⁻², (1+e₂)⁻², (1+e₃)⁻²])
        %
        function Er = computeEr(obj, K1, K2, K3)
            K0 = (K1 * K2 * K3)^(1/3);
            a = obj.config.SlateCoeffA;
            b = obj.config.SlateCoeffB;
            e1 = exp(a * ((K1/K0)-1) - b) - 1;
            e2 = exp(a * ((K2/K0)-1) - b) - 1;
            e3 = exp(a * ((K3/K0)-1) - b) - 1;
            omega = (1+e1)*(1+e2)*(1+e3);
            Er = omega^2 * diag([(1+e1)^-2, (1+e2)^-2, (1+e3)^-2]);
        end
        

        %computeV - 計算eigen vector
        %
        %輸入參數：trend & plunge(dk1geo, dk2geo, dk3geo, id1geo, ik2geo, ik3geo)
        %
        %輸出參數:
        %   V- eigen vector
        %
        %計算方法：V = [cos(i1)*cos(d1), cos(i2)*cos(d2), cos(i3)*cos(d3); %N = X
        %             cos(i1)*sin(d1), cos(i2)*sin(d2), cos(i3)*sin(d3); %E = Y
        %             sin(i1),         sin(i2),         sin(i3)];        %D = Z
        function V = computeV(~, d1, i1, d2, i2, d3, i3)
            V = [cosd(i1)*cosd(d1), cosd(i2)*cosd(d2), cosd(i3)*cosd(d3); %N = X
                 cosd(i1)*sind(d1), cosd(i2)*sind(d2), cosd(i3)*sin(d3); %E = Y
                 sind(i1),         sind(i2),         sind(i3)];        %D = Z
            threshold = 1e-12;
            V(abs(V) < threshold) = 0;
        end
        

        %computeEg - 計算轉回地理座標系統之橢球
        %
        %輸入參數：V,Er
        %
        %計算方法：Eg = V'(V^T) * Er * V;
        %
        
     
        function Eg = computeEg(~, Er, V)
            Eg = V' * Er * V;
        end
        %驗證model，不做磁感率轉應變
        %
        function computeErRawAll(obj)
            K1 = obj.data.K1; K2 = obj.data.K2; K3 = obj.data.K3;
            n = length(K1);
            obj.Er_raw = zeros(3,3,n);

            for i = 1:n
                Er = diag([K1(i), K2(i), K3(i)]);
                obj.Er_raw(:,:,i) = Er;
                assignin('base', sprintf('Er_raw_%d', i), Er);
                if obj.config.Verbose
                    fprintf('Er_raw_%d = diag([%.4e, %.4e, %.4e])\n', i, K1(i), K2(i), K3(i));
                end
            end
        end
        %驗證model，不做磁感率轉應變
        %
        function computeEgFromErRaw(obj)
            dK1 = deg2rad(obj.data.dK1geo); iK1 = deg2rad(obj.data.iK1geo);
            dK2 = deg2rad(obj.data.dK2geo); iK2 = deg2rad(obj.data.iK2geo);
            dK3 = deg2rad(obj.data.dK3geo); iK3 = deg2rad(obj.data.iK3geo);
            n = size(obj.Er_raw, 3);
            obj.EgRaw = zeros(3,3,n);

            for i = 1:n
                V = obj.computeV(dK1(i), iK1(i), dK2(i), iK2(i), dK3(i), iK3(i));
                Eg = V' * obj.Er_raw(:,:,i) * V;
                obj.EgRaw(:,:,i) = Eg;
                assignin('base', sprintf('Eg_raw_%d', i), Eg);
                if obj.config.Verbose
                    fprintf('Eg_raw_%d 已完成\n', i);
                end
            end
        end
        
        function [V_sorted, D_sorted] = eigSorted(~, A)
            % eigSorted - 對矩陣 A 做特徵分解並依特徵值由大到小排序
            %
            % 輸入：
            %   A - 方陣
            %
            % 輸出：
            %   V_sorted - 排序後的特徵向量矩陣（每一列是特徵向量）
            %   D_sorted - 排序後的特徵值對角矩陣（由大到小排序）
            
             [V, D] = eig(A);
             eigVals = diag(D);
             [sortedEigVals, idx] = sort(eigVals, 'descend');  % 由大到小排序
             D_sorted = diag(sortedEigVals);
             V_sorted = V(:, idx);
        end
        

    end
end


%2D分析
%{
        function computeFiniteStrainTensors2D(obj)
            K1 = obj.data.K1; K2 = obj.data.K2;
            dK1 = deg2rad(obj.data.dK1geo); iK1 = deg2rad(obj.data.iK1geo);
            dK2 = deg2rad(obj.data.dK2geo); iK2 = deg2rad(obj.data.iK2geo);
            a = obj.config.SlateCoeffA; b = obj.config.SlateCoeffB;

            n = length(K1);
            obj.EgD = zeros(2,2,n);
            obj.ErList2D = zeros(2,2,n);
            obj.VList2D = zeros(2,2,n);

            for i = 1:n
                K0 = sqrt(K1(i)*K2(i));
                e1 = exp(a*((K1(i)/K0)-1) - b) - 1;
                e2 = exp(a*((K2(i)/K0)-1) - b) - 1;
                omega = (1+e1)*(1+e2);
                Er2D = omega^2 * diag([(1+e1)^-2, (1+e2)^-2]);

                v1 = [cos(iK1(i))*cos(dK1(i)); cos(iK1(i))*sin(dK1(i))];
                v2 = [cos(iK2(i))*cos(dK2(i)); cos(iK2(i))*sin(dK2(i))];
                V2D = [v1, v2];

                obj.ErList2D(:,:,i) = Er2D;
                obj.VList2D(:,:,i) = V2D;
                obj.EgD(:,:,i) = V2D' * Er2D * V2D;
            end
        end
%}