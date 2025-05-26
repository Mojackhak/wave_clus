function key = findFirstKey(list, regexPattern)
% helper: return first fieldname that matches regex, or [] if none
    idx = find(~cellfun(@isempty, regexp(list, regexPattern, 'once')), 1);
    if isempty(idx)
        key = [];
    else
        key = list{idx};
    end
end
