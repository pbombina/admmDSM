%Convert Adjacency List into Adjacency Matrix%
lists = strsplit(fileread('jazz.txt'), {'\r', '\n'});
if isempty(lists{end}), lists = lists(1:end-1); end
lists = cellfun(@(v) sscanf(v, '%d'), lists, 'UniformOutput', false);
adj = sparse(repelem(1:numel(lists), cellfun('length', lists)), vertcat(lists{:}), 1)

full(adj)
