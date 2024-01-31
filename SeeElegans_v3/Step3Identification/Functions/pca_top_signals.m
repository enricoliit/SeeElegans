function top_signals = pca_top_signals(neurons,top_sel,zs,final_time)
if nargin < 3
    zs = 0;
end
    if zs
        signal_data = [];
        for i_neu= 1:numel(neurons)
            signal_data(:,i_neu) = zscore(neurons(i_neu).coords(10:final_time,5));
        end
        
        [coeff, ~, ~, ~, explained] = pca(signal_data);
        pc1 = coeff(:,1);
        abs_weights = abs(pc1);
        [~, sorted_indices] = sort(abs_weights, 'descend');
        top_signals = sorted_indices(1:top_sel);
        
    else
        signal_data = [];
        for i_neu= 1:numel(neurons)
            signal_data(:,i_neu) = neurons(i_neu).coords(10:final_time,5);
        end
        
        [coeff, ~, ~, ~, explained] = pca(signal_data);
        pc1 = coeff(:,1);
        abs_weights = abs(pc1);
        [~, sorted_indices] = sort(abs_weights, 'descend');
        top_signals = sorted_indices(1:top_sel);
        
    end
        
   
end