A = readmatrix("../data/matrix_a.csv");
B = readmatrix("../data/matrix_b.csv");
C = readmatrix("../data/matrix_c.csv");

tA = readtable("../data/matrix_a.csv");
tB = readtable("../data/matrix_b.csv");
tC = readtable("../data/matrix_c.csv");

function [cleaned] = clean_input(m)
    m(:, 1) = [];
    sz = size(m, 1);
    m(1:sz+1:end) = 0;
    cleaned = m;
end

function [table] = convert_to_table(sorted_m, indices, prev_table)
    names = prev_table.Properties.VariableNames
    names(:, 1) = []
    reordered_names = names(indices)
    T_sorted = array2table(sorted_m);
    T_sorted.Properties.VariableNames = reordered_names;
    T_sorted = addvars(T_sorted, reordered_names', 'Before', 1)
    table = T_sorted
end



A = clean_input(A);

%disp("total energy for time period 1")
%disp("ranks for matrix A")
%ranks = rank(A);
%[energy, p_val] = energy_test(A, ranks, 1000)
%disp("Energy per edge Matrix A:")


ranksA = rank(A)
[XSorted, Ind] = sort(ranksA, 'descend');
ASorted = A(Ind, Ind);
new_table_A = convert_to_table(ASorted, Ind, tA)

ranksB = rank(B)
[XSorted, Ind] = sort(ranksB, 'descend');
BSorted = B(Ind, Ind);
new_table_B = convert_to_table(BSorted, Ind, tB)




%T_sorted
writetable(new_table_A, "../data/matrix_a_sorted.csv")
writematrix(ranks, "../data/ranks_a.csv")



%disp("total energy for time period 2")
%B = clean_input(B)
%ranks = springRank(B);
%get_tot_energy(ranks, B)


%disp("total energy for time period 3")
%C = clean_input(C)
%ranks = springRank(C);
%get_tot_energy(ranks, C)