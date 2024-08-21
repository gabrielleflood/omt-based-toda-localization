function u_out = min_logsumexp_l1_ball(v_in,eta,N_cum_vec)
%
% Solves the problem
%
% min_u exp(-u+v)'*1 , s.t ||u||_1 \leq \eta
%
% where it is assumed that all entries of the vector v are strictly
% positive. This is without loss of generality as if v_k = 0 for some k,
% then the corresponding optimal u_k = 0, and these elements can be
% excluded from the optimization procedure.
%
% Solves the problem by Lagrange relaxation. Solution found by simple
% sorting.
%
%
%%%%%%%%%
% Date of creation: 19 October 2023. Filip Elvander
%%%%%%%%
% This version: 22 November 2023. Filip Elvander
%%%%%%%%

v = sort(v_in,'ascend');
N = length(v);

if nargin<3
    N_cum_vec = (N-(1:N)');
end

%%%% Check if all variables are active
if sum(v)-N*v(1)<eta
    add_const = sum(v)-eta;
    c = N;
    opt_lagrange_multiplier = add_const/c;
else

    %%%% Assume at least one variable is zero in optimum %%%.
    c_v = cumsum(v);
    u_tilde = c_v(end)-c_v-v.*N_cum_vec;

    upper_index = find(u_tilde<eta,1,'first');

    %%% this is the number of non-zero variables and also slope of function in
    %%% this interval.
    c = N-upper_index+1;

    % Off-set in this interval
    add_const = (c_v(end)-c_v(upper_index-1))-eta;

    % Find zero of linear function. This is the Lagrange multiplier
    opt_lagrange_multiplier = add_const/c;
end

u_out = max(v_in-opt_lagrange_multiplier,0);
end







