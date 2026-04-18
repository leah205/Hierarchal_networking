% generates vector of ranks with ith rank being the rank of ith node in
% matrix (higher rank = more dominant)
function [ranks] = springRank(m)
    
    %row sums (# interactions where i dominates)
    out_deg = sum(m, 2);
    %col sum(# interactions where i is dominated)
    in_deg = sum(m, 1)';
   
    % left hands side of system
    B = diag(out_deg) + diag(in_deg) - (m + m');
    
    %right hand side of system
    b = out_deg - in_deg;

    % solve system for ranks
    ranks = lsqminnorm(B, b);
    ranks = ranks - mean(ranks);
end



