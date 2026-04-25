% simplified version of sigma_a from the paper
function [sigma] = accuracy_simp(test_matrix, train_matrix)
    s = rank(train_matrix);
    % s = s(randperm(length(s))); % can be used for evaluating statistical significance of rank()
    n = length(s);
    number_correct = 0; % the number of edge directions correctly predicted
    for i = 1 : n
        for j = 1 : n
            if (test_matrix(i,j) > 0) && (s(i) > s(j)) 
                number_correct = number_correct + test_matrix(i,j); 
            end
        end     
    end
    sigma = number_correct / sum(sum(test_matrix));
end

% full version of sigma_a from the paper
function [sigma] = accuracy_full(test_matrix, train_matrix)
    s = rank(train_matrix);
    n = length(s);
    
    function [a] = negacc(b)
        t = sum(sum(train_matrix));
        y = 0;
        for i = 1:n
            for j = 1:n
                p = (1 + exp(-2*b*(s(i)-s(j))))^(-1);
                y = y + abs(train_matrix(i, j) - (train_matrix(i,j) + ...
                    train_matrix(j, i))*p);
            end
        end
        a = y/t - 1;
    end

    b = fminbnd(@negacc, 1e-6, 1000);
    
    t = sum(sum(test_matrix));
    y = 0;
    for i_ = 1:n
        for j_ = 1:n
            p = (1 + exp(-2*b*(s(i_)-s(j_))))^(-1);
            y = y + abs(test_matrix(i_,j_) - (test_matrix(i_,j_) + test_matrix(j_,i_))*p);
        end     
    end
    sigma = 1 - 0.5*y/t;
end

% returns two arrays with each trial's accuracy
function [sigma_s, sigma_f] = cross_validate(A, n_reps, n_folds)
    if ~issparse(A)
        A = sparse(A);
    end

	% initialize arrays of results
    sigma_s = zeros(n_reps*n_folds, 1);
    sigma_f = zeros(n_reps*n_folds, 1);

    % row and column subscripts of each interaction
	% and a vector v containing each interaction
    [row, col, v] = find(triu(A + transpose(A)));
    fold_size = floor(length(v)/n_folds); 

    for r = 1:n_reps
		% randomly shuffle the indices of the entries of v
        indices = randperm(length(v));

        % split A into n_folds number of folds. 
        for f = 1:n_folds - 1
            fold{f} = indices((f-1)*fold_size + 1 : f*fold_size);
        end
        %  store the remainder in the final entry
        fold{n_folds} = indices((n_folds - 1)* fold_size + 1 : end);

        % train and test the model n_folds number of times, cycling through 
        % each of the folds taking the role of the test set
        for f = 1:n_folds
            test_i = row(fold{f});
            test_j = col(fold{f});
            test_indices = sub2ind(size(A),test_i, test_j);
            test_transpose_indices = sub2ind(size(A),test_j,test_i);

            train_matrix = A;
            train_matrix(test_indices) = 0;
            train_matrix(test_transpose_indices) = 0;

            test_matrix = A - train_matrix;

            sigma_s((r-1)*n_folds + f,1) = accuracy_simp(test_matrix, train_matrix);
            sigma_f((r-1)*n_folds + f,1) = accuracy_full(test_matrix, train_matrix);
        end    
    end    
end

function [upset_rate] = upset(A, s) 
    total_interactions = 0;
    upsets = 0;
    n = length(s);

    for i = 1:n
        for j = 1:n
            if A(i,j) > 0  % interactions where i dominated j
                if s(i) < s(j)  % i is lower ranked than j
                    upsets = upsets + A(i,j);
                end
                total_interactions = total_interactions + A(i,j);
            end
        end
    end
    
    upset_rate = upsets / total_interactions;
end


A = readmatrix("../data/matrix_a_sorted.csv");
A(:,1) = [];
[A_simp, A_full] = cross_validate(A, 50, 5);
fprintf("Simplified accuracy for matrix A with 5 folds: %.4f \n", mean(A_simp));
fprintf("Full accuracy for matrix A with 5 folds: %.4f \n", mean(A_full));

% test accuracy, this should be close to 1
accuracy(A, A);

% confirm that accuracy is lower with a lower fold size and higher with a higher fold size
[A_simp, A_full] = cross_validate(A, 10, 2);
fprintf("Simplified accuracy for matrix A with 2 folds: %.4f \n", mean(A_simp));
fprintf("Full accuracy for matrix A with 2 folds: %.4f \n", mean(A_full));

[A_simp, A_full] = cross_validate(A, 10, 20);
fprintf("Simplified accuracy for matrix A with 20 folds: %.4f \n", mean(A_simp));
fprintf("Full accuracy for matrix A with 20 folds: %.4f \n", mean(A_full));

B = readmatrix("../data/matrix_b_sorted.csv");
B(:,1) = [];
[B_simp, B_full] = cross_validate(B, 50, 5);
fprintf("Average simplified accuracy for matrix B: %.4f \n", mean(B_simp));
fprintf("Average full accuracy for matrix B: %.4f \n", mean(B_full));

C = readmatrix("../data/matrix_c_sorted.csv");
C(:,1) = [];
[C_simp, C_full] = cross_validate(C, 50, 5);
fprintf("Average simplified accuracy for matrix C: %.4f \n", mean(C_simp));
fprintf("Average full accuracy for matrix C: %.4f \n", mean(C_full));

% matrix of males at time period A
mA = A;
mA([5 7 10:15 17:18 20:22], :) = [];
mA(:, [5 7 10:15 17:18 20:22]) = [];
[mA_simp, ~] = cross_validate(mA, 50, 5);
fprintf("Average accuracy for males at time period A: %.4f \n", mean(mA_simp));

% matrix of females at time period A
fA = A;
fA([1:4 6 8:9 16 19], :) = [];
fA(:, [1:4 6 8:9 16 19]) = [];
[fA_simp, ~] = cross_validate(fA, 50, 5);
fprintf("Average accuracy for females at time period A: %.4f \n", mean(fA_simp));


s = rank(A);
s([5 7 10:15 17:18 20:22], :) = [];
sM = rank(mA);
[rM, pM] = corr(sM,s,"Type", "Spearman");

%{
figure('Position',[200 200 800 500]);
scatter(sM, s)
xlabel("Male-only Rank")
ylabel("Full Rank")
title("Correlation of Male Ranks")
lsline()
%}

s = rank(A);
s([1:4 6 8:9 16 19], :) = [];
sF = rank(fA);
[rF, pF] = corr(sF,s,"Type", "Spearman");

%{
figure('Position',[200 200 800 500]);
scatter(sF, s)
xlabel("Female-only Rank")
ylabel("Full Rank")
title("Correlation of Female Ranks")
lsline()
%}

fprintf("Correlation between male-only hierarchy ranks and full hierarchy ranks is %.4f with a p-value of %.4f\n",rM,pM);
fprintf("Correlation between female-only hierarchy ranks and full hierarchy ranks is %.4f with a p-value of %.4f\n",rF,pF);

fprintf("Upset rate among males: %.4f\n", upset(mA, sM));
fprintf("Upset rate among females: %.4f\n", upset(fA, sF));
fprintf("Total upset rate: %.4f\n", upset(A, rank(A)));


% Use A as the training set to predict B
% remove GS which is not in A and reorder
B2 = readmatrix("../data/matrix_b_reordered.csv");
B2(:,1) = [];
rAB = accuracy_simp(B2, A);
fprintf("%.4f\n", mean(mean(rAB)));

%{
% scatter plots 

figure('Position',[200 200 800 500]);
scatter(A_simp, A_full)
lsline
title('Accuracy of A rankings')
xlabel('Simplified \sigma_a')
ylabel('\sigma_a')

figure('Position',[200 200 800 500]);
scatter(B_simp, B_full)
lsline
title('Accuracy of B Rankings')
xlabel('Simplified \sigma_a')
ylabel('\sigma_a')

figure('Position',[200 200 800 500]);
scatter(C_simp, C_full)
lsline
title('Accuracy of C Rankings')
xlabel('Simplified \sigma_a')
ylabel('\sigma_a')
%}
