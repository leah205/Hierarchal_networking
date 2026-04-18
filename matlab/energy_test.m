
function [energy, p_val] = test(A, ranks, n)  
    energy = get_tot_energy(ranks, A) / sum(sum(A));
    random_energies = get_random_energies(A, 1000);
    p_val = sum(random_energies <= energy) + 1 / 10001;
  
end

function [random] = generate_randomized(m)
    sz = size(m, 1);
    random = zeros(sz);

    for i = 1:sz
        for j = 1:sz
            interactions = m(i, j) + m(j, i);
            random_i_dom = randi([0, interactions]);
            random(i, j) = random_i_dom;
            random(j, i) = interactions - random_i_dom;
        end
    end
end

function [avg_random_energy_per_edge] = get_random_energies(m, n)
    sum_energies = 0;
    for i = 1:n
        random = generate_randomized(m);
        random_ranks = rank(random);
        sum_energies = sum_energies + get_tot_energy(random_ranks, random) / sum(sum(random));
    end
    avg_random_energy_per_edge = sum_energies / n;
    disp("Average random energy per edge for matrix A:")
    disp(avg_random_energy_per_edge);
end

%A(i, j) has # of fearful signals from j when interacting with i
%( A(i, j) -> i dominates j

function [energy] = get_tot_energy(ranks, m)
H = zeros(size(m, 1)); % Initialize the H matrix    
%create H matrix storing energy of each spring
    for(i = 1:size(m, 1))
        for(j = 1:size(m, 1))
            H(i, j) = 0.5 * (ranks(i) - ranks(j) - 1)^2;
        end
    end
    energy = 0;

    %adds up spring energies to find total energy
    for(i = 1:size(m, 1))
        for(j = 1:size(m, 1))
            energy = energy + H(i, j) * m(i, j);
        end
    end
end



