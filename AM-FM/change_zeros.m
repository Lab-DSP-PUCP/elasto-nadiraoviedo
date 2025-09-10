function matrix_no_zero = change_zeros(matrix)
    % Copia la matriz original
    matrix_no_zero = matrix;
    
    % Encuentra los índices de los valores que son cero
    [rows, cols] = find(matrix == 0);
    
    % Recorre cada valor cero
    for k = 1:length(rows)
        row = rows(k);
        col = cols(k);
        
        % Busca el valor más cercano diferente de cero en la fila
        row_values = matrix(row, :);
        row_values_diff_zero = row_values(row_values ~= 0);
        
        % Busca el valor más cercano diferente de cero en la columna
        col_values = matrix(:, col);
        col_values_diff_zero = col_values(col_values ~= 0);
        
        % Si ambos valores existen, promedialos. Si solo uno existe, úsalo.
        if ~isempty(row_values_diff_zero) && ~isempty(col_values_diff_zero)
            new_value = (row_values_diff_zero(1) + col_values_diff_zero(1)) / 2;
        elseif ~isempty(row_values_diff_zero)
            new_value = row_values_diff_zero(1);
        elseif ~isempty(col_values_diff_zero)
            new_value = col_values_diff_zero(1);
        else
            new_value = 0; % Si no hay valores no cero alrededor
        end
        
        % Reemplaza el valor cero en la matriz
        matrix_no_zero(row, col) = new_value;
    end
end
