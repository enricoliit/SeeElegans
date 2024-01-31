% function distmatrix_s=makeDistMatrix_s(coords_t,maxtime_s)
function distmatrix_s = makeDistMatrix_s(coords_t)
    num_neurons = size(coords_t,3);
    distmatrix_s = zeros(num_neurons, num_neurons, 1, 3);
    h = waitbar(0,'Building distance matrix ...');
    for i = 1:num_neurons
        waitbar(i/num_neurons,h);
        for j = i+1:num_neurons
            distmatrix_s(i, j, 1, :) = mean(coords_t(:, 1:3,i) - coords_t(:, 1:3,j), 1, 'omitnan');
            distmatrix_s(j, i, 1, :) = -distmatrix_s(i, j, 1, :);
        end
    end
    close(h)
end