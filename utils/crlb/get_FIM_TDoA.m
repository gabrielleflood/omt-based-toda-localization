function partial_FIM = get_FIM_TDoA(source_pos,rec1_pos,rec2_pos,prop_speed,noise_var)
%
% Computes the Fisher Information matrix for the localization of one source
% from a time difference of arrival(TDOA) measurement, i.e., the
% information contained about the source coordinates. It is assumed that
% the TDoA measurement is a Gaussian random variable.
% 
% Ambient dimension is R^d (for d = 2 or 3).
%
% INPUT
% source_pos            -           source position coordinates, d-vector.
% rec1_pos              -           coordinates for first receiver,
%                                   d-vector.
% rec2_pos              -           coordiantes for second receiver,
%                                   d-vector.
% prop_speed            -           propagation speed of the medium,
%                                   scalar.
% noise_var             -           variance of the measurement noise,
%                                   scalar.
%
% OUTPUT
% partial_FIM           -           the Fisher Information matrix, 
%                                   d-times-d matrix of rank 1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version: 4 August 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = source_pos;
x1 = rec1_pos;
x2 = rec2_pos;
rho = prop_speed;
sigma2 = noise_var;

vec1 = s-x1;
vec2 = s-x2;
norm1 = norm(vec1);
norm2 = norm(vec2);

vvec = vec1/norm1-vec2/norm2;

partial_FIM = (vvec*vvec')/(sigma2*rho^2);

end
