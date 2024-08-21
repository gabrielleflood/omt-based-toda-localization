function cost_mat = setup_cost_matrix(tdoas_listed,candidates,r,trash_prct,trash_cost)
% compute the cost matrix, here computed as 
% c = (||x-r_k|| - ||x-r_l|| - tau)^2
% with x being candidate position, r_k and r_l the corresponding receiver
% positions, ||.|| representing 2-norm and tau being measured tdoa
% and the thrash bin cost is a constant, here set to be the trash_prct 
% percentile of all values in the rest of the cost matrix, if trash_prct is
% input argument (not empty) and the constant trash_cost if this variable 
% is not empty. 

if nargin<4 || nargin <5 || (isempty(trash_prct) && isempty(trash_cost))
    error('trash_prct or trash_cost has to be used, both need to be sent in')
end
if ~isempty(trash_prct) && ~isempty(trash_cost)
    error('Only trash_prct or trash_cost can be used')
end

cost_mat = zeros(length(tdoas_listed), size(candidates,2)+1);
for i = 1:length(tdoas_listed)
    % compute cost for every candidate in row i, corresponding to tdoa i
    % exclude dustbin cost for now
    candidates_tdoa = sqrt(sum((candidates-repmat(r(:,tdoas_listed(i).r2),[1 size(candidates,2)])).^2,1)) -...
                        sqrt(sum((candidates-repmat(r(:,tdoas_listed(i).r1),[1 size(candidates,2)])).^2,1));
%     candidates_tdoa = vecnorm( candidates-repmat(r(:,tdoas_listed(i).r1),[1 size(candidates,2)]) ) - ...
%                         vecnorm( candidates-repmat(r(:,tdoas_listed(i).r2),[1 size(candidates,2)]) );
    cost_mat(i,1:size(candidates,2)) = ( candidates_tdoa - tdoas_listed(i).tdoas ).^2;

end

% set the cost for the trash bin. Now, set the cost at all tdoas to be the
% 30th percentile of all above computed costs.
if ~isempty(trash_prct)
    all_costs = cost_mat(i,1:size(candidates,2));
    trash_cost = prctile(all_costs(:),trash_prct);
end   
cost_mat(:,end) = trash_cost;

