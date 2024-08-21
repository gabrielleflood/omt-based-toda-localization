function new_pos = remove_duplicates(old_pos,thresh)
% Compare the different N 3D positions in the 3xN matrix old_pos and
% removes any duplicates. Two positions are considered a duplicate if their
% Euclidean distance is smaller than thresh

if nargin < 2 || isempty(thresh); thresh = 0.1; end

% find the distances
dist_matrix = pdist2(old_pos',old_pos');

% only look at upper triangular distances, so set the rest to NaN
dist_matrix = dist_matrix + tril(inf(size(old_pos,2)));

[row,col] = find(dist_matrix < thresh);

remove_ind = row';
keep_ind = 1:size(old_pos,2);
keep_ind(remove_ind) = [];

new_pos = old_pos(:,keep_ind);

% if length(remove_ind)>0
%     old_pos
%     new_pos
%     keyboard();
% end