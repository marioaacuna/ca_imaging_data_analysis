function metadata = process_experimental_condition(metadata, criteria, criteria_to_match)

if isempty(metadata), return, end
if isempty(criteria_to_match), return, end

global GC

if isfield(GC, 'experiment_name') && ~isempty(GC.experiment_name)
    % Make a new table with the conditions split
    if ~isempty(GC.columns_experimental_condition.(GC.experiment_name))
        conds = translate_experimental_condition(metadata.experimental_condition);
        add_column_names = GC.columns_experimental_condition.(GC.experiment_name);

        if size(add_column_names, 1) == 1
            conds = cell2table(conds(:, 1:length(add_column_names)), 'VariableNames',add_column_names);

        else
            all_columns_to_add = add_column_names(:, 2);
            all_columns_to_add = unique([all_columns_to_add{:}]);
            all_conds = cell(size(conds, 1), length(all_columns_to_add));
            columns_to_keep = {};

            % Add all columns
            n_rules = size(add_column_names, 1);
            for i_rule = 1:n_rules
                rule_to_match = criteria_to_match{i_rule};
                value_to_match = criteria{i_rule}{2};
                rows_this_criterion = ismember(metadata{:, rule_to_match}, value_to_match);
                if ~any(rows_this_criterion), continue, end

                columns_this_criterion = add_column_names{i_rule, 2};
                columns_to_keep = [columns_to_keep; columns_this_criterion(:)];
                values = conds(rows_this_criterion, :);
                for i_col = 1:length(columns_this_criterion)
                    col_idx = ismember(all_columns_to_add, columns_this_criterion{i_col});
                    all_conds(rows_this_criterion, col_idx) = values(:, i_col);
                end
            end
            % Discard unused columns
            columns_to_keep_idx = ismember(all_columns_to_add, unique(columns_to_keep));
            conds = cell2table(all_conds(:, columns_to_keep_idx), 'VariableNames',all_columns_to_add(columns_to_keep_idx));
        end

        try
            % Insert conds table next to experimental condition
            column_names = metadata.Properties.VariableNames;
            idx_exp_cond = find(ismember(column_names, 'experimental_condition'));
            metadata = [metadata(:, 1:idx_exp_cond), conds, metadata(:, idx_exp_cond+1:end)];
        end
    end
end


%#ok<*AGROW>
