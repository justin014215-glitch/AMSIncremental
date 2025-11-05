function [V_sorted, D_sorted] = eigSorted(A)
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
