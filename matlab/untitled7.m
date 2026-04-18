
function [ranks] = springRank(m)
    if(~issparse(m))
        m = sparse(m);
    end
    sz = size(m);
    %row sums (# interactions where i dominates)
    out_deg = sum(A, 2);
    %col sum(# interactions where i is dominated)
    in_deg = sum(B, 2);

    %create diagonal matrix of out degrees for each node except last(fixed
    %rank)
    Diag_out = diag(out_deg(1: sz - 1));
    %create diagonal matrix of in degrees for each node except last(fixed
    %rank)
    Diag_in = diag(in_deg(1: sz - 1));

    % get last node in and out degree
    out_last = out_deg(sz);
    in_last = in_deg(sz);
    
    %fills matrix with last node dominance interactions
    last_node_dom = repmat(A(sz, 1: sz - 1), sz - 1, 1);
    %fills matrix up with last node column
    last_node_sub = repmat(A(1: sz - 1, sz)', sz - 1, 1);
    
    %mtrix without fixed last node
    last_excl_m =  A(1: sz - 1, 1: sz - 1)
    
    lhs = Diag_out + Diag_in - last_excl_m - last_excl_m' - last_node_dom - last_node_sub
    rhs = Diag_out - Diag_in + out_last - in_last;

    %use iterative solver
    [s, ~] = bicgstab(lhs, rhs, 1e-12, 200)
    %set ranks to iterative solver solution, last rank fixed at zero
    ranks = [s;0]

    fixed_m = [];



end


springRank(A);