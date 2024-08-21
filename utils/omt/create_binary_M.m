function M = create_binary_M(M_in, best_pos)
% takes the matrix M and changes it s.t. it is a proper omt plan
% if best_pos, the best candiates from ind vec, is sent in, it is ensured 
% that only these activated candidates are used

if nargin<2 || isempty(best_pos)
    best_pos = 1:size(M_in,2);
end

% create index matrix saying which columns that can be used
ok_cols = zeros(1,size(M_in,2));
ok_cols(best_pos) = 1;

% set all M values 
M_in(:,~ok_cols) = 0;

[~,I] = max(M_in,[],2);

M = zeros(size(M_in));
M(sub2ind(size(M),1:size(M,1),I')) = 1;

end