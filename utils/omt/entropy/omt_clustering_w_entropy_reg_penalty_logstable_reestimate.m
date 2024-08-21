
function [M,trash_mass] = omt_clustering_w_entropy_reg_penalty_logstable_reestimate(C,trash_cost,epsilon,R)
%
%
% Routine for reestimating the associations. 
%
% Specifically, the problem solved is
%
% min_{M,m_trash} trace(C^T M) + eta*norm(M)_{\infty,1} + epsilon*D(M)
%               + c_trash*1_\Tau^T m_trash + epsilon*D(m_trash)
%
% subject to
%           M 1_\Omega + m_trash = 1_\Tau , M^T 1_\Tau <= R(R-1)/2 1_\Omega
%
% where \Omega is the set of candidate targets and \Tau is the set of
% TDoAs, with 1_\Omega and 1_\Tau being vectors all-ones with length equal
% to cardinalities of \Omega and \Tau, respectively.
%
% D(M) = \sum_{k,ell} m_{k,ell}log(m_{k,ell} - m_{k,ell}
% is the discrete entropy. This automatically gives non-negative solutions.
%
% The first constraint ensures that the elements are <= 1 in the optimal
% point. Together with D(M), this gives optimal M with elements in [0,1].
%
% The penalty norm(M)_{\infty,1} is the sum of the infinity norms of the
% columns of M, i.e., norm(M)_{\infty,1} = \sum_{ell} max(M(:,ell)). This
% promotes the clustering.
%
% For stability, the computations are carried out in log-domain.
%
%%%%% Additional routines needed %%%%%%
% min_logsumexp_l1_ball     -       sub-routine required for implementing
%                                   the sparsity penalty.
% stable_log_sum_exp        -       function for computing log-sum-exp in a
%                                   stable way. Required for log-stable
%                                   computations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%% Input %%%%%%%%%%%
% C             -       cost matrix, should have size |\Tau|\times|\Omega|.
%                       Note: this SHOULD NOT include the trash cost.
% trash_cost    -       scalar, cost for moving to trash state. Could be
%                       vector also.
% epsilon       -       strictly positive scalar. This is the entropy
%                       regularization parameter. The smaller this
%                       parameter, the closer to the original OMT problem.
%                       Smaller value will effect convergence speed, but
%                       should not affect numerical stability.
% R             -       number of receivers.
%
%%%%%%% Output %%%%%%%%%%
% M             -       transport plan, matrix of size
%                       |\Tau|\times|\Omega|. Indicates which TDoA is
%                       associated with which target. Elements of transport
%                       plan should be between zero and one in optimum.
%                       number of TDoAs. Note: cluster_vec = M^T 1_\Tau.
% trash_mass    -       vector of length |\Tau| containing the mass
%                       transported to the trash state.
%
%%%%%%%%%%%
% Date of first version: 13 February 2024. Filip Elvander.
%%%%%%%%%%%%
% This version: 13 February 2024. Filip Elvander.
%%%%%%%%%%%

sol_tol = 1e-5;
max_iter = 1e4;

save_dual_val = 0;
if save_dual_val
    dual_hist = zeros(max_iter,1);
    primal_hist = zeros(max_iter,1);
end


%%%% Constraint violation tolerances %%%%%%
box_tol = 1e-5;
mass_balance_tol = 1e-2;
upper_bound_tol = 1e-5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N_tau,N_grid] = size(C);

%%%%% Make trash cost into vector %%%%
if length(trash_cost)==1
    trash_cost = trash_cost*ones(N_tau,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Dual variables %%%%%%
lambda = zeros(N_tau,1);
mu = zeros(N_grid,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Pre-computed quantities %%%%%%
scaled_C = 1/epsilon*C;
ones_mu = ones(N_grid,1);
ones_mu_plus_trash = ones(N_grid+1,1);
ones_lambda = ones(N_tau,1);
trash_scaled = -1/epsilon*trash_cost;%*ones(N_tau,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Other constants %%%%%%
R_tilde = R*(R-1)/2;
mu_constant = -log(R_tilde)*ones(N_grid,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k_iter = 1:max_iter

    if save_dual_val
        [pp,dd] = eval_primal_dual_entropy_w_trash_reest(C,trash_cost,epsilon,R,lambda,mu);
        dual_hist(k_iter) = dd;
        primal_hist(k_iter) = pp;
    end

    %%%%% Update lambda %%%%%%
    log_xi_lambda = stable_log_sum_exp([-scaled_C-1/epsilon*(ones_lambda*mu'),trash_scaled],ones_mu_plus_trash);
    lambda = -epsilon*log_xi_lambda;
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% Update mu %%%%%%%%%%
    log_xi_mu = stable_log_sum_exp(-scaled_C'+1/epsilon*(ones_mu*lambda'),ones_lambda);
    mu = epsilon*max(0,mu_constant+log_xi_mu);
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%% Diagnostics based on duality gap %%%%%%%%%%%%
    if k_iter>2000 && mod(k_iter,105)==0
        [pp,dd,duality_gap] = eval_primal_dual_entropy_w_trash_reest(C,trash_cost,epsilon,R,lambda,mu);
        duality_gap = duality_gap/dd;
        fprintf('Iteration %d , duality gap %.2E\n',k_iter,duality_gap)
        if duality_gap<sol_tol
            % Check feasibility %
            log_M = (-scaled_C+(1/epsilon*lambda-1/epsilon*mu'));
            M = exp(log_M);
            mass_balance_constraint = max(abs(sum(M,2)-1));
            if mass_balance_constraint<1e-6
                fprintf('Convergence after %d iterations\n',k_iter)
                break
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%% Compute transport plan %%%%%%%%
% The entries of M should be in [0,1] in the optimal point. log(M) can be
% computed in a stable manner, whereafter M can be computed.
log_M = (-scaled_C+(1/epsilon*lambda-1/epsilon*mu'));
M = exp(log_M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Compute mass transported to trash state %%%%%
trash_mass = exp(-1/epsilon*trash_cost+1/epsilon*lambda);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Check feasibility of solution %%%%%
box_constraint = max(max(max(M))-1,0); % Is M constained on [0,1]?
mass_balance_constraint = max(abs(sum(M,2)+trash_mass-1)); % Is the "mass" of the measurements accounted for?
upper_bound_constraint = max(max(sum(M,1)'-R_tilde,0)); % Is a target assigned more than R(R-1)/2 measurements?

fprintf('\n\n')
if box_constraint>box_tol
    fprintf('Box constraint violation: %.2E\n',box_constraint)
end
if mass_balance_constraint>mass_balance_tol
    fprintf('Mass balance constraint violation: %.2E\n',mass_balance_constraint)
end
if upper_bound_constraint>upper_bound_tol
    fprintf('Upper bound constraint violation: %.2E\n',upper_bound_constraint)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



