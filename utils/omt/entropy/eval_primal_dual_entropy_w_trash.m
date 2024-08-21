function [val_primal,val_dual,duality_gap] = eval_primal_dual_entropy_w_trash(C,trash_cost,epsilon,eta,R,lambda,mu,Psi)

log_M = (-1/epsilon*C+1/epsilon*Psi+(1/epsilon*lambda-1/epsilon*mu'));
M = exp(log_M);
log_m = -1/epsilon*trash_cost+1/epsilon*lambda;
m = exp(log_m);

val_dual = -epsilon*(sum(M(:))+sum(m))-...
    R*(R-1)/2*sum(mu)+sum(lambda);

val_primal = C(:)'*M(:) + trash_cost(:)'*m(:) + ...
    epsilon*(disc_ent(M,log_M)+disc_ent(m,log_m))+...
    eta*sum(max(M));

duality_gap = val_primal-val_dual;

end

function val_out = disc_ent(M,log_M)
    val_out = M(:)'*log_M(:)-sum(M(:));
end


