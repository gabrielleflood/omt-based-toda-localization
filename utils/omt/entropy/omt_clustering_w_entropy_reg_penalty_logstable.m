
function [M,cluster_vec] = omt_clustering_w_entropy_reg_penalty_logstable(C,epsilon,eta,R)
%
%
% Solve clustering problem where we do not specify number of targets but
% instead promote sparse solution.
%
% Specifically, the problem solved is
%
% min_{M} trace(C^T M) + eta*norm(M)_{\infty,1} + epsilon*D(M)
%
% subject to
%           M 1_\Omega = 1_\Tau , M^T 1_\Tau <= R(R-1)/2 1_\Omega
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
% epsilon       -       strictly positive scalar. This is the entropy
%                       regularization parameter. The smaller this
%                       parameter, the closer to the original OMT problem.
%                       Smaller value will effect convergence speed, but
%                       should not affect numerical stability.
% eta           -       sparsity penalty parameter. The larger, the fewer
%                       columns of the optimal transport plans will be
%                       non-zero.
% R             -       number of receivers.
%
%%%%%%% Output %%%%%%%%%%
% M             -       transport plan, matrix of size
%                       |\Tau|\times|\Omega|. Indicates which TDoA is
%                       associated with which target. Elements of transport
%                       plan should be between zero and one in optimum.
% cluster_vec   -       vector of length |\Omega|. This is the summation of
%                       M along columns, saying how much TDoA "mass" is
%                       assigned to each candidate target. Each element
%                       is constrained to be bounded from above by R(R-1)/2
%                       as each target can be associated with at most that
%                       number of TDoAs. Note: cluster_vec = M^T 1_\Tau.
%
%%%%%%%%%%%
% Date of first version: 22 November 2023. Filip Elvander.
%%%%%%%%%%%%
% This version: 12 February 2024. Filip Elvander.
%%%%%%%%%%%

sol_tol = 1e-8;
max_iter = 3e4;
Psi_update_cycle = 20;

%%%% Constraint violation tolerances %%%%%%
box_tol = 1e-5;
mass_balance_tol = 1e-2;
upper_bound_tol = 1e-5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N_tau,N_grid] = size(C);

%%%%% Dual variables %%%%%%
lambda = zeros(N_tau,1);
mu = zeros(N_grid,1);
Psi = zeros(N_tau,N_grid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Pre-computed quantities %%%%%%
scaled_C = 1/epsilon*C;
ones_mu = ones(N_grid,1);
ones_lambda = ones(N_tau,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Other constants %%%%%%
R_tilde = R*(R-1)/2;
mu_constant = -log(R_tilde)*ones(N_grid,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Psi_minus_C = 1/epsilon*Psi-scaled_C;

Psi_old = Psi;
for k_iter = 1:max_iter
    
    %%%%% Update lambda %%%%%%
    log_xi_lambda = stable_log_sum_exp(Psi_minus_C-1/epsilon*(ones_lambda*mu'),ones_mu);
    lambda = -epsilon*log_xi_lambda;
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% Update mu %%%%%%%%%%
    log_xi_mu = stable_log_sum_exp(Psi_minus_C'+1/epsilon*(ones_mu*lambda'),ones_lambda);
    mu = epsilon*max(0,mu_constant+log_xi_mu);
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Update Psi %%%%%%%%%%
    if mod(k_iter,Psi_update_cycle)==0 || k_iter ==1
        log_Xi = -scaled_C + (1/epsilon*lambda-1/epsilon*mu');
        for ell = 1:N_grid
            Psi(:,ell) = -epsilon*min_logsumexp_l1_ball(log_Xi(:,ell),eta/epsilon);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    Psi_minus_C = 1/epsilon*Psi-scaled_C;

    diagnostic_val = norm(Psi_old-Psi,"fro")/norm(Psi_old,'fro');
    if mod(k_iter,Psi_update_cycle)==0 && diagnostic_val<sol_tol
        fprintf('Convergence after %d iterations\n',k_iter)
        break
    end
    
    if mod(k_iter,100)==0
        fprintf('Iteration %d , rel. norm change %.2E\n',k_iter,diagnostic_val)
    end
    Psi_old = Psi;
end

%%%% Final iteration for consistency %%%%
%%%%% Update lambda %%%%%%
log_xi_lambda = stable_log_sum_exp(Psi_minus_C-1/epsilon*(ones_lambda*mu'),ones_mu);
lambda = -epsilon*log_xi_lambda;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Update mu %%%%%%%%%%
log_xi_mu = stable_log_sum_exp(Psi_minus_C'+1/epsilon*(ones_mu*lambda'),ones_lambda);
mu = epsilon*max(0,mu_constant+log_xi_mu);
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Compute transport plan %%%%%%%%
% The entries of M should be in [0,1] in the optimal point. log(M) can be
% computed in a stable manner, whereafter M can be computed.
log_M = (-scaled_C+1/epsilon*Psi+(1/epsilon*lambda-1/epsilon*mu'));
M = exp(log_M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Compute indicator vector %%%%%%
cluster_vec = sum(M,1)'/R_tilde;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Check feasibility of solution %%%%%
box_constraint = max(max(max(M))-1,0); % Is M constained on [0,1]?
mass_balance_constraint = max(abs(sum(M,2)-1)); % Is the "mass" of the measurements accounted for?
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



