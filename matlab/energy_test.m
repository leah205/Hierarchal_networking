
function [energy, p_val] = test(A, ranks, n)  
    tot_energy = get_tot_energy(ranks, A);
    energy = tot_energy / sum(sum(A));
    disp("total energy for A");
    tot_energy;
    [avg_random, random_energies] = get_random_energies(A, 1000);

    figure;
    histogram(random_energies, 15);
    hold on;

    xline(tot_energy, 'r', 'LineWidth', 2);
    xlabel('Ground State Energy');
    ylabel('Frequency');
    title('Ground State Energy of hierarchy compared with Randomized Networks');
    legend('Randomised Network Distribution', 'Computed Energy');

 
    p_val = (sum(random_energies <= energy) + 1) / (n + 1);
  
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

function [avg_random_energy_per_edge, random_energies] = get_random_energies(m, n)
    random_energies = []
    sum_energies = 0;
    sum_tot_energies = 0;
    for i = 1:n
        random = generate_randomized(m);
        random_ranks = rank(random);
        random_energy = get_tot_energy(random_ranks, random);
        sum_energies = sum_energies +  random_energy / sum(sum(random));
        sum_tot_energies = sum_tot_energies + random_energy;
        random_energies(end + 1) = random_energy;
    end
    avg_random_energy_per_edge = sum_energies / n;
    disp("Average random energy per edge for matrix A:")
    disp(avg_random_energy_per_edge);
     disp("Average random energy  for matrix A:")
    disp(sum_tot_energies / n);
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



