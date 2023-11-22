function thresholded_matrix = proportional_threshold_top(matrix, percentile)
    threshold_value = prctile(abs(matrix(:)), 100 - percentile);
    thresholded_matrix = matrix;
    thresholded_matrix(abs(matrix) < threshold_value) = 0;
end