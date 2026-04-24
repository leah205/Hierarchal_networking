% read in adjacency matrices for each time period

A = readmatrix("../data/matrix_a.csv");
B = readmatrix("../data/matrix_b.csv");
C = readmatrix("../data/matrix_c.csv");

tA = readtable("../data/matrix_a.csv");
tB = readtable("../data/matrix_b.csv");
tC = readtable("../data/matrix_c.csv");

% Process matrix for analysis
function [cleaned] = clean_input(m)
    %remove names column
    m(:, 1) = [];
    sz = size(m, 1);
    % Nan on diagonals converted to zero
    m(1:sz+1:end) = 0;
    cleaned = m;
end

% convert adjacency matrix and ranks vectors to tables
function [table, ranks] = convert_to_table(sorted_m, indices, prev_table, ranks)
    names = prev_table.Properties.VariableNames
    names(:, 1) = []
    reordered_names = names(indices);
    T_sorted = array2table(sorted_m);
    T_sorted.Properties.VariableNames = reordered_names;
    T_sorted = addvars(T_sorted, reordered_names', 'Before', 1);

    ranks = array2table(ranks);
    ranks = addvars(ranks, reordered_names', 'Before', 1);
    table = T_sorted;
end



A = clean_input(A);
B = clean_input(B);
C = clean_input(C);

%disp("total energy for time period 1")
%disp("ranks for matrix A")
%ranks = rank(A);
%disp("Energy per edge Matrix A:")

% find ranks for each time period 
%sort adjacency matrix and ranks vector by ranks descending
ranksA = rank(A);
[energy, p_val] = energy_test(A, ranksA, 1000);
[rASorted, Ind] = sort(ranksA, 'descend');
ASorted = A(Ind, Ind);
[new_table_A, new_ranks_A] = convert_to_table(ASorted, Ind, tA, rASorted);

ranksB = rank(B);
[rBSorted, Ind] = sort(ranksB, 'descend');
BSorted = B(Ind, Ind);
[new_table_B, new_ranks_B] = convert_to_table(BSorted, Ind, tB, rBSorted);

ranksC = rank(C);
[rCSorted, Ind] = sort(ranksC, 'descend');
CSorted = C(Ind, Ind);
[new_table_C, new_ranks_C] = convert_to_table(CSorted, Ind, tC, rCSorted);




%T_sorted
writetable(new_table_A, "../data/matrix_a_sorted.csv");
writetable(new_ranks_A, "../data/ranks_a.csv");

writetable(new_table_B, "../data/matrix_b_sorted.csv");
writetable(new_ranks_B, "../data/ranks_b.csv");

writetable(new_table_C, "../data/matrix_c_sorted.csv");
writetable(new_ranks_C, "../data/ranks_c.csv");



%disp("total energy for time period 2")
%B = clean_input(B)
%ranks = springRank(B);
%get_tot_energy(ranks, B)


%disp("total energy for time period 3")
%C = clean_input(C)
%ranks = springRank(C);
%get_tot_energy(ranks, C)