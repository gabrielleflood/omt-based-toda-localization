function CRLB_mean = get_CRLB_mean(r,s,R,S,sigma)

% Create cell with all receiver-receiver pairs
rec1_rec2_pair_cell = cell(R*(R-1)/2,1);
counter = 1;
for k1 = 1:R-1
    for k2 = k1+1:R
        rec1_rec2_pair_cell{counter} = [r(:,k1),r(:,k2)];
        counter = counter+1;
    end
end

prop_speed_nominal = 1;
sigma2_noise_nominal = sigma^2;
for k_source = 1:S
    [CRLB_source_pos,CRLB_mat,tot_FIM] = ...
        get_CRLB_source_pos(s(:,k_source),rec1_rec2_pair_cell,prop_speed_nominal,sigma2_noise_nominal);
    CRLB_simulation_mat(k_source) = sqrt(trace(CRLB_mat));
    
    if sqrt(trace(CRLB_mat))>1000
        %'stop'
    end
end
CRLB_mean = mean(CRLB_simulation_mat(:));