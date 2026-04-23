% simplified version of sigma_a from the paper
function [sigma] = accuracy(test_matrix, train_matrix)
    s = rank(train_matrix);
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

% returns a matrix with each trial's accuracy
function [results] = cross_validate(A, n_reps, n_folds)
    if ~issparse(A)
        A = sparse(A);
    end

	% initialize matrix of results
    results = zeros(n_reps, n_folds);

    % row and column subscripts of each interaction
	% and a vector v containing each interaction
    [row, col, v] = find(triu(A + transpose(A)));
    fold_size = floor(length(v)/n_folds); 

    for r = 1:n_reps
		% randomly shuffle the entries of v
        rv = v(randperm(length(v)));

		% split A into n_folds number of folds
        for f = 1:n_folds - 1
            fold{f} = rv((f-1)*fold_size + 1 : f*fold_size);
        end
		% store the remainder in the final entry
        fold{n_folds} = rv((n_folds - 1)* fold_size + 1 : end);

        % train and test the model n_folds number of times, cycling through 
        % each of the folds take the role of the test set
        for f = 1:n_folds
            test_i = row(fold{f});
            test_j = col(fold{f});
            test_indices = sub2ind(size(A),test_i, test_j);
            test_transpose_indices = sub2ind(size(A),test_j,test_i);

            train_matrix = A;
            train_matrix(test_indices) = 0;
            train_matrix(test_transpose_indices) = 0;

            test_matrix = A - train_matrix;

            results(r,f) = accuracy(test_matrix, train_matrix);
        end    
    end    
end

% test accuracy, this should be close to 1
accuracy(A, A);

A = readmatrix("matrix_a_sorted.csv");
A(:,1) = [];
rA = cross_validate(A, 50, 5);
fprintf("%.2f\n", mean(mean(rA)));

B = readmatrix("matrix_b_sorted.csv");
B(:,1) = [];
rB = cross_validate(B, 50, 5);
fprintf("%.2f\n", mean(mean(rB)));

C = readmatrix("matrix_c_sorted.csv");
C(:,1) = [];
rC = cross_validate(C, 50, 5);
fprintf("%.2f\n", mean(mean(rC)));
