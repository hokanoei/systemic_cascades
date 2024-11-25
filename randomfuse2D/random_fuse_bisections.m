format longEng

% 12.14 on 07.07.23

Rfuse = 1;
Rsmall = 10^9;
eps = 10^-10;
sigma = 1; 
threshold = 1e-3; 
Nsamples = 5; 

% Set the desired accuracy
desiredAccuracy = 1e-3;

L_vec = [500]; 
dist_vals = [1e-2, 1e-3, 1e-4]; 

for L = L_vec
    
    I_th_vec = []; 

    for runs = 1 : Nsamples 
        
        disp([L runs])
       
        % Fixing the disordered seed at this value of L
        N = L * L + 2;
        [Ic_zero1, Aij1] = set_Ic_zero(L, N, sigma);  
        orderParameter = @(Ib) (computeOrderParameter(Rfuse, Rsmall, eps, Ib, L, Ic_zero1, Aij1)); 
    
        % Define the function that calculates the order parameter R for a given Ib
        %orderParameter = @(Ib) (computeOrderParameter(Rfuse, Rsmall, eps, Ib, L, sigma)); 
    
        % Define the initial interval [a, b] within which to search for the critical threshold
        Ith = 0.5 * L / log(L); 
        a = 0.5 * Ith;
        b = 2 * Ith;
    
        % Perform the bisection search
        while (b - a) > desiredAccuracy
            
            c = (a + b) / 2; % Calculate the midpoint
            [G_c, ~, ~] = orderParameter(c);
    
            % Check if the order parameter crosses the threshold
            if G_c < threshold
                b = c;
            else
                a = c;
            end
        end
    
        % The critical threshold is given by the midpoint of the final interval
        criticalThreshold = (a + b) / 2;
        I_th_vec = [I_th_vec criticalThreshold];
        
        %disp([criticalThreshold])
        %disp('Finished the bisection seach -- Evaluating the evolution')
        
        
        %figure; 
        %hold on; 

        %i = 0; 
        for dist = dist_vals
            
            [G_evo, iter_evo, count_evo] = computeOrderParameterEvolution(Rfuse, Rsmall, eps, criticalThreshold, dist, L, Ic_zero1, Aij1); 

            if length(G_evo) > 0 && G_evo(end) < 0.3

                count_evo = full(count_evo); 
                S1fuse = count_evo(2 : end) - count_evo(1 : end - 1);
                filename = ['/Users/ivan_bonamassa/Desktop/random_fuse/fuse_cascades/S_t_L', num2str(L), 'Ic_delta', num2str(dist), '_test', num2str(runs), '.mat'];
                save(filename, 'S1fuse');
            
            end
            
            %i = i + 1; 
            %c = i / length(dist_vals); 

            %colormap(jet);  % Set the colormap to "jet"
            %cMap = colormap;  % Get the colormap values
            %color = cMap(round(c*size(cMap,1)), :);  % Get the color from the colormap based on the scaled index

            %label = ['dist ' num2str(dist)];  % Generate a label based on the loop variable
            %plot(iter_evo, G_evo, '-o', 'Color', color, 'DisplayName', label);
            %ylim([0.0, 1.0]);

        end

        %legend('show');
        %legend('Location', 'best');

        % Save every time a run is finished
        filename = ['/Users/ivan_bonamassa/Desktop/random_fuse/random_fuse_Ith_L', num2str(L), '.mat'];
        save(filename, 'I_th_vec'); 

    end
    
end


% Function to calculate the order parameter R for a given Ib and L

function [G, state, NOIs] = computeOrderParameter(Rfuse, Rsmall, eps, Ib, L, Ic_zero1, Aij1)
%function [G, state, NOIs] = computeOrderParameter(Rfuse, Rsmall, eps, Ib, L, sigma)

    N = L * L + 2;
    iterNOI = N;
    W1 = zeros(N,1);
    W1(N) = 0;
    Ib_vec1 = zeros(N-1,1);

    % Comment if you want to keep the same random seed
    %[Ic_zero1, Aij1] = set_Ic_zero(L, N, sigma);  
    Ic_T1 = Ic_zero1; 
    R_ij1 = Rfuse * Aij1;  

    flag1 = 0;

    R1 = Rfuse;
    Ib_vec1(1) = Ib;

    for iter = 1 : iterNOI
        
        if flag1 == 0

            [G1, state1] = set_G(W1, L, N, R_ij1, Rsmall, Ic_T1, Rfuse); % random fuse model
            
            G21 = G1(1:N-1, 1:N-1);
            W_next1 = zeros(N,1);
            W_next1(1:N-1) = G21 \ Ib_vec1;
            W_next1(N) = 0;
            R1 = W_next1(1) / Ib;
            W_pre1 = W1;
            R1pre = W1(1) / Ib;
            W1 = W_next1;
    
            epsilon = abs(R1 - R1pre) / R1;
    
            if epsilon < eps
                break;
            end

        end
        
        if 1 / R1 < 10^-9 && flag1 == 0
            flag1 = 1;
        end

        if flag1 
            break
        end

        if any(isnan(W_next1))
            break;
        end

    end

    G = 1 / R1;
    state = state1; 
    NOIs = iter; 

end


function [G_evo, iter_evo, count_evo] = computeOrderParameterEvolution(Rfuse, Rsmall, eps, Ib, dist, L, Ic_zero1, Aij1)
%function [G, state, NOIs] = computeOrderParameter(Rfuse, Rsmall, eps, Ib, L, sigma)

    G_cas_evo = []; 
    iter_cas_evo = []; 
    count_cas_evo = []; 

    N = L * L + 2;
    iterNOI = 100;
    W1 = zeros(N,1);
    W1(N) = 0;
    Ib_vec1 = zeros(N-1,1);

    % Comment if you want to keep the same random seed
    %[Ic_zero1, Aij1] = set_Ic_zero(L, N, sigma);  
    Ic_T1 = Ic_zero1; 
    R_ij1 = Rfuse * Aij1;  

    flag1 = 0;

    R1 = Rfuse;
    Ib_vec1(1) = Ib + dist;

    for iter = 1 : iterNOI
        
        if flag1 == 0

            [G1, state1] = set_G(W1, L, N, R_ij1, Rsmall, Ic_T1, Rfuse); % random fuse model
            
            G21 = G1(1:N-1, 1:N-1);
            W_next1 = zeros(N,1);
            W_next1(1:N-1) = G21 \ Ib_vec1;
            W_next1(N) = 0;
            R1 = W_next1(1) / Ib;
            W_pre1 = W1;
            R1pre = W1(1) / Ib;
            W1 = W_next1;
    
        end
        
        if 1 / R1 < 10^-3 && flag1 == 0
            flag1 = 1;
        end

        if flag1 
            break
        end

        if any(isnan(W_next1))
            break;
        end

        count = count_states(state1);

        G_cas_evo = [G_cas_evo 1 / R1]; 
        iter_cas_evo = [iter_cas_evo iter]; 
        count_cas_evo = [count_cas_evo count]; 

    end
    
    G_evo = G_cas_evo; 
    iter_evo = iter_cas_evo; 
    count_evo = count_cas_evo; 

end


function count = count_states(state)
    count = sum(sum(state == 3))/2;
end

