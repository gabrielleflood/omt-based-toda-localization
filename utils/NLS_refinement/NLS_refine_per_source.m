function s_est = NLS_refine_per_source(temp_tdoa_rec_mat,r,s0)

nbr_tdoas = size(temp_tdoa_rec_mat,1);
R1_mat = zeros(3,nbr_tdoas);
R2_mat = zeros(3,nbr_tdoas);
tdoa_vec = temp_tdoa_rec_mat(:,1)';
for kk = 1:nbr_tdoas
    r1_ind = temp_tdoa_rec_mat(kk,2);
    r2_ind = temp_tdoa_rec_mat(kk,3);
    r1 = r(:,r1_ind);
    r2 = r(:,r2_ind);
    R1_mat(:,kk) = r1;
    R2_mat(:,kk) = r2;
end

f = @(x) sum((sqrt(sum((R2_mat-x).^2))-sqrt(sum((R1_mat-x).^2))-tdoa_vec).^2);

options = optimset('TolX',1e-10,'TolFun',1e-10);
s_est = fminsearch(f,s0,options);

%%%% OBS! with bounds %%%%%
%LB = s0 - ones(size(s0));
%UB = s0 + ones(size(s0));
%s_est = fminsearchbnd(f,s0,LB,UB,options);
%s_est = fminsearch(f,s0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end