function [CRLB_source_pos,CRLB_mat,tot_FIM] = get_CRLB_source_pos(source_pos,rec1_rec2_pair_cell,prop_speed,noise_var_vec)

% Computes the Cramer-Rao lower bound for the estimate of a source position
% in R^d (for d = 2 or 3) given a set of time difference of arrival (TDoA)
% measurements and known receiver positions.
%
% INPUT
% source_pos                -       source position coordinates, d-vector.
% rec1_rec2_pair_cell       -       cell array containing all pairs of
%                                   receivers that correspond to a TDoA
%                                   measurement. Each cell element should
%                                   be a d-times-2 matrix where the columns
%                                   are the coordinates of the respective
%                                   receiver positions for the
%                                   corresponding pair.
% prop_speed                -       propagation speed of the medium, scalar
%
% noise_var_vec             -       variance of the measurement noise. This
%                                   can be given as a vector if the noise
%                                   is assumed to have different variance
%                                   for different TDoA measurements.
%                                   In this case, the length of the vector 
%                                   should be the same as that of the cell
%                                   array rec1_rec2_pair_cell. If the noise
%                                   is isotropic, then noise_var_vec can be
%                                   given as a scalar.
%
% OUTPUT
% CRLB_source_pos           -       the CRLB for the source position
%                                   estimate, i.e., a lower bound on the
%                                   variance. For an unbiased estimator,
%                                   this serves as a bound on the
%                                   expectation of || s - \hat{s} ||_2^2
%                                   where s \in R^d is the source position
%                                   and \hat{s} is an estimator. Scalar
%                                   value.
% CRLB_mat                  -       CRLB matrix, serving as a matrix lower
%                                   bound on the covariance matrix of the
%                                   estimated coordinate vector, i.e, 
%                                   CRLB_mat-E((s-\hat{s})(s-\hat{s})^T)>=0
%                                   This is a d-times-d matrix.
% FIM                       -       Fisher information matrix for the
%                                   measurements.
%                                   Note: FIM^{-1} = CRLB_mat.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version: 4 August 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbr_TDoA_pairs = length(rec1_rec2_pair_cell);
if length(noise_var_vec) ==1
    noise_var_vec = noise_var_vec*ones(nbr_TDoA_pairs,1);
elseif length(noise_var_vec) ~=nbr_TDoA_pairs
    error('noise_var_vec must either scalar or vector of consistent length')
end

% Compute Fisher information matrix
tot_FIM = 0;
for k_TDoA = 1:nbr_TDoA_pairs
    rec_pos_mat = rec1_rec2_pair_cell{k_TDoA};
    rec1_pos = rec_pos_mat(:,1);
    rec2_pos = rec_pos_mat(:,2);
    % Compute FIM for this particular receiver-receiver pair. 
    partial_FIM = get_FIM_TDoA(source_pos,rec1_pos,rec2_pos,prop_speed,noise_var_vec(k_TDoA));

    tot_FIM = tot_FIM + partial_FIM;
end

CRLB_mat = inv(tot_FIM);
CRLB_source_pos = trace(CRLB_mat);

end
