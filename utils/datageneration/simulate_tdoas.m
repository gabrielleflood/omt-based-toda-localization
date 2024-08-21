function [tdoas_measured,tdoas_true,r,s,nbr_missing, nbr_extra] = simulate_tdoas(R,S,sigma,r_bounds, s_bounds,P_missing,P_extra,all_tdoas)
% Input: R>4 is nbr of receivers, S nbr of senders, r_bounds and s_bounds 
% are 3xn matrices with bounds for simulations of r and s, respectively.
% P_missing is the probability that there is a TDOA measurement missing and
% P_extra is the probability that there is an extra, erroneous TDOA
% measurement. sigma is noise of std
% Output: tdoa_measurements are the TDOAs, given in some form (TBD). r are
% the receiver positions and s are the sender positions. 
% If R and/or S is not an integer but a matrix, it shoukd be a 3xR or 3xS
% matrix with the positions of the receivers/senders.
% all_tdoas is 1 if all possible tdoas from all toas should be computed and
% 0 if we assume that the correct tdoa peaks are found (only the ones with
% same sender are computed)

% set values that might not be sent as input
if nargin < 3 || isempty(sigma); sigma = 0; end
if nargin < 4 || isempty(r_bounds); r_bounds = [0 10; 0 10; 0 2]; end
if nargin < 5 || isempty(s_bounds); s_bounds = [0 10; 0 10; 0 2]; end
if nargin < 6 || isempty(P_missing); P_missing = 0; end
if nargin < 7 || isempty(P_extra); P_extra = 0; end
if nargin < 8 || isempty(all_tdoas); all_tdoas = 0; end

% define the dimension
dim = 3; 

% simulate receivers if not included as argument, otherwise, save in
% correct variables
if numel(R) == 1
    % simulate receiver positions using rbounds
    r = rand(3,R);
    r(1,:) = r_bounds(1,1) + (r_bounds(1,2)-r_bounds(1,1)).*r(1,:);
    r(2,:) = r_bounds(2,1) + (r_bounds(2,2)-r_bounds(2,1)).*r(2,:);
    r(3,:) = r_bounds(3,1) + (r_bounds(3,2)-r_bounds(3,1)).*r(3,:);
elseif size(R,1) ~= dim
    error('Input receivers must be given as a 3xR matrix')
else
    r = R;
    R = size(r,2);
end


% simulate sender if not included as argument, otherwise, save in
% correct variables
if numel(S) == 1
    % simulate sender positions using s_bounds
    s = rand(3,S);
    s(1,:) = s_bounds(1,1) + (s_bounds(1,2)-s_bounds(1,1)).*s(1,:);
    s(2,:) = s_bounds(2,1) + (s_bounds(2,2)-s_bounds(2,1)).*s(2,:);
    s(3,:) = s_bounds(3,1) + (s_bounds(3,2)-s_bounds(3,1)).*s(3,:);
elseif size(S,1) ~= dim
    error('Input sender must be given as a 3xS matrix')
else
    s = S;
    S = size(s,2);
end


%% compute true TDOAs

% first, compute toas
toas_true = pdist2(r',s'); % this gives r:s on rows, s:s on columns. 

if ~all_tdoas
% Compute tdoas from toas. Here, we assume that we have the correct tdoas, 
% s.t. only toas for the same microphone are subtracted. 
% OBS! Silly solution, improve with matrix computations later.
tdoas_true = cell(R,R);
for i = 1:R
    for j = i+1:R
        tdoas_true{i,j}.r1 = i;
        tdoas_true{i,j}.r2 = j;
        tdoas_true{i,j}.s = 1:S;
        tdoas_true{i,j}.N = S;
        tdoas_true{i,j}.tdoas = toas_true(j,:)-toas_true(i,:);
        % do we want to sort the tdoas? Does it matter?
%         [tdoas_true{i,j}.tdoas,I] = sort(tdoas_true{i,j}.tdoas);
%         tdoas_true{i,j}.s =  tdoas_true{i,j}.s(I);
    end
end
tdoas_measured = tdoas_true;
end

if all_tdoas
% here, we assume that we do not have the correct tdoas, but that there are
% erroneous tdoas that are coming from matching of receiver peaks that
% originate from different senders.
tdoas_all = cell(R,R);
for i = 1:R
    for j = i+1:R
        tdoas_all{i,j}.r1 = i;
        tdoas_all{i,j}.r2 = j;
        tdoas_all{i,j}.tdoas = repmat(toas_true(j,:),[S 1]) ...
                                -repmat(toas_true(i,:)',[1 S]);
% OBS! Silly solution, improve with matrix computations later.
        tdoas_all{i,j}.tdoas = tdoas_all{i,j}.tdoas(:)';
        
        tdoas_all{i,j}.s1 = repmat(1:S,[S 1]);
        tdoas_all{i,j}.s1 = tdoas_all{i,j}.s1(:);
        tdoas_all{i,j}.s2 = repmat((1:S)',[1 S]);
        tdoas_all{i,j}.s2 = tdoas_all{i,j}.s2(:);
        tdoas_all{i,j}.N = length(tdoas_all{i,j}.tdoas);
        % do we want to sort the tdoas? Does it matter?
%         [tdoas_true{i,j}.tdoas, I] = sort(tdoas_true{i,j}.tdoas);
%         tdoas_all{i,j}.s1 =  tdoas_all{i,j}.s1(I);
%         tdoas_all{i,j}.s2 =  tdoas_all{i,j}.s2(I);
    end
end

tdoas_measured = tdoas_all;
tdoas_true = tdoas_all;
end

%% compute measured TDOAs

% create a variable for measured tdoas

nbr_missing = 0;
nbr_extra = 0;
% add or remove measurements according to P_missing and P_extra
if P_missing>0 || P_extra>0
    for i = 1:R
        for j = i+1:R
            if rand(1) < P_missing
                % decide which tdoa measurement to remove using randperm
                pos = randperm(length(tdoas_measured{i,j}.tdoas),1);
                tdoas_measured{i,j}.tdoas(pos) = [];
                if ~all_tdoas
                    tdoas_measured{i,j}.s(pos) = [];
                else
                    tdoas_measured{i,j}.s1(pos) = [];
                    tdoas_measured{i,j}.s2(pos) = [];
                end
                tdoas_measured{i,j}.N = tdoas_measured{i,j}.N-1;
                nbr_missing = nbr_missing+1;
            end
            if rand(1) < P_extra
                % append an extra measurement
                val = mean(tdoas_measured{i,j}.tdoas) + std(tdoas_measured{i,j}.tdoas)*randn(1); 
                tdoas_measured{i,j}.tdoas = [tdoas_measured{i,j}.tdoas val];
                if ~all_tdoas
                    tdoas_measured{i,j}.s(end+1) = NaN;
                else
                    tdoas_measured{i,j}.s1(end+1) = NaN;
                    tdoas_measured{i,j}.s2(end+1) = NaN;
                end
                tdoas_measured{i,j}.N = tdoas_measured{i,j}.N+1;
                nbr_extra = nbr_extra+1;
            end
            % do we want to sort the tdoas? Does it matter?
            % tdoas_true{i,j}.tdoas = sort(tdoas_true{i,j}.tdoas);
        end
    end
end


% add noise using sigma
for i = 1:R
    for j = i+1:R
        tdoas_measured{i,j}.tdoas = tdoas_measured{i,j}.tdoas + sigma*randn(size(tdoas_measured{i,j}.tdoas));
    end
end

