%% Random fuse model - 2D disordered lattices
% Iterative Kirchoff equations for the random fuse model

format longEng

Rfuse = 1; 
Rsmall = 10^9;
eps = 10^-10;

%iterNOI = 2.5*10^1;

Imin = 10^-6; 
dI = 0.1; 

sigma = 1; 
L_vec = [2^5];
I_th_vec = []; % save the thresholds

figure(); 
for L = L_vec

    %% parameters
    N = L*L + 2;
    iterNOI = L;
    W1 = zeros(N,1);
    W1(N) = 0;
    %R1 = Rfuse;
    Ib_vec1 = zeros(N-1,1);

    %% initialization
    R_heating_coupled_vec1 = [];
    Ib_vals = [];
    iter_heating_coupled_vec = [];
    
    [Ic_zero1, Aij1] = set_Ic_zero(L, N, sigma);  %% comment this if you want to load
    Ic_T1 = Ic_zero1; 
    R_ij1 = Rfuse * Aij1;    %% Fuses have all the same resistance
    flag1 = 0;
    
    %% random fuse solve
    Ith = 0.5 * L / log(L); 
    Imin = 0.5 * Ith; 
    Imax = 2 * Ith; 
    
    for Ib = Imin : dI : Imax

        R1 = Rfuse; 
        Ib_vec1(1) = Ib;

        for iter = 1 : iterNOI

            if flag1 == 0
                [L Ib iter R1]
                    
                [G1, state1] = set_G(W1, L, N, R_ij1, Rsmall, Ic_T1, Rfuse); % random fuse model
                G21 = G1(1:N-1, 1:N-1);
                W_next1 = zeros(N,1);
                W_next1(1:N-1) = G21 \ Ib_vec1;
                W_next1(N) = 0;
                R1 = W_next1(1) / Ib;
                W_pre1 = W1;
                R1pre = W1(1) / Ib;
                W1 = W_next1;
                  
                %epsilon = (norm(W_next1 - W_pre1)) / (norm(W_next1));
                epsilon = abs(R1 - R1pre) / R1;
                   
                if epsilon < eps
                    break;
                end

            end

            if 1 / R1 < 10^-7 && flag1 == 0
                flag1 = 1;
            end

            if flag1 
                break
            end

            if any(isnan(W_next1))
                break;
            end

        end

        Ib_vals = [Ib_vals Ib];
        R_heating_coupled_vec1 = [R_heating_coupled_vec1 1 / R1];
        iter_heating_coupled_vec = [iter_heating_coupled_vec iter];

        if flag1 
            break
        end

    end

    hold on; 
    plot(Ib_vals, R_heating_coupled_vec1, '-o');
    plot([Ith], [0], '-o', color = 'r');
    ylim([0, 1.0]);
    xlim([Imin, Imax]);
    ylabel('G(S)');
    xlabel('I(A)');
    I_th_vec = [I_th_vec Ib]; 
filename = ['/Users/ivan_bonamassa/Desktop/random_fuse_Ith_L.mat'];
save(filename, 'L_vec', 'I_th_vec'); 

hold off; 

end
