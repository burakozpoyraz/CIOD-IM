%==========================================================================
% Burak Özpoyraz, 2020
% Ergodic Secrecy Rate of Efficient Alamouti Scheme

% ARGUMENTS
% 1-) num_iterations: Number of iterations for Monte Carlo simulation
% 2-) M: Constellation size
% 3-) P_tot_des: Desired total power to be transmitted during an iteration
% 4-) theta: Ratio of the power allocated to information symbols
% 5-) SNRdB: Signal-to-noise ratio in dB scale

% OUTPUTS
% 1-) Rb: Ergodic rate of the legitimate receiver
% 2-) Re: Ergodic rate of the eavesdropper
% 3-) Rs: Ergodic secrecy rate
%==========================================================================

%% MAIN FUNCTION
function [Rb, Re, Rs] = EfficientAlamoutiSecrecy(num_iterations, M, P_tot_des, alpha, SNRdB)
    %% PARAMETERS
    P_X_des = P_tot_des * alpha; % Desired power of the Alamouti matrix
    P_symbol_des = P_X_des / 4; % Desired power of a single symbol in the constellation
    P_Z_des = P_tot_des * (1 - alpha); % Desired power of the AN matrix
    
    if M > 8
        ss = qammod(0 : M-1, M, 'Gray');
    else
        ss = pskmod(0 : M-1, M, 0, 'Gray');
    end
    Pc = ss * ss'/ M; % Average power of a symbol in the constellation
    ss = ss * sqrt(P_symbol_des / Pc);  % Consellation power normalization: E[|xk|^2] = 1
    X = AllPossibleMatrices(M, ss);
    
    N0 = P_symbol_des / 10^(SNRdB / 10); % The noise power
    
    Rb_array = zeros(1, num_iterations); % Array of Ergodic Rates of Bob
    Re_array = zeros(1, num_iterations); % Array of Ergodic Rates of Eve
    Rs_array = zeros(1, num_iterations); % Array of Ergodic Secrecy Rates
    
    %% SIMULATION
    for iter_index = 1 : num_iterations
        fprintf("Iteration Index: %d\n", iter_index);
        
        % Channel & Noise//////////////////////////////////////////////////
        h = (randn(1, 3) + 1i * randn(1, 3)) / sqrt(2); % Channel of the legitimate receiver
        g = (randn(1, 3) + 1i * randn(1, 3)) / sqrt(2); % Channel of the eavesdropper
        
        nb = sqrt(N0 / 2) * (randn(1, 2) + 1i * randn(1, 2)); % Noise of the legitimate receiver
        ne = sqrt(N0 / 2) * (randn(1, 2) + 1i * randn(1, 2)); % Noise of the eavesdropper
        % /////////////////////////////////////////////////////////////////
        
        % Ergodic Rates and Secrecy Rate///////////////////////////////////
        i_sum_bob = 0;
        i_sum_eve = 0;
        for i = 1 : 2
            hi = [h(i) h(3)];
            gi = [g(i) g(3)];
            
            n_sum_bob = 0;
            n_sum_eve = 0;
            for n = 1 : M^2
                Xn = X(:, :, n);
                
                i2_sum_bob = 0;
                i2_sum_eve = 0;
                for i2 = 1 : 2
                    hi2 = [h(i2) h(3)];
                    gi2 = [g(i2) g(3)];
                    
                    n2_sum_bob = 0;
                    n2_sum_eve = 0;
                    for n2 = 1 : M^2
                        Xn2 = X(:, :, n2);
                        
                        % Artificial Noise/////////////////////////////////
                        A1 = i;
                        Z = AN_Matrix(h, A1);
                        P_Z = sum(sum(abs(Z).^2));
                        Z_N = Z * sqrt(P_Z_des / P_Z); % Artificial noise normalization
                        % /////////////////////////////////////////////////
                        
                        % Linear Whitening Transformation//////////////////
                        ne_hat = g * Z_N + ne;
                        Ne_hat = ne_hat * ne_hat' / 2;
                        Psie = sqrt(N0) * (Ne_hat)^(-1 / 2);
                        
                        gi_tilde = Psie * gi;
                        gi2_tilde = Psie * gi2;
                        ne_tilde = Psie * ne_hat;
                        % /////////////////////////////////////////////////
                        
                        value_bob = exp(-(norm(hi * Xn - hi2 * Xn2 + nb)^2 - norm(nb)^2) / N0);
                        n2_sum_bob = n2_sum_bob + value_bob;
                        
                        value_eve = exp(-(norm(gi_tilde * Xn - gi2_tilde * Xn2 + ne_tilde)^2 - norm(ne_tilde)^2) / N0);
                        n2_sum_eve = n2_sum_eve + value_eve;
                    end
                    i2_sum_bob = i2_sum_bob + n2_sum_bob;
                    i2_sum_eve = i2_sum_eve + n2_sum_eve;
                end
                n_sum_bob = n_sum_bob + log2(i2_sum_bob);
                n_sum_eve = n_sum_eve + log2(i2_sum_eve);
            end
            i_sum_bob = i_sum_bob + n_sum_bob;
            i_sum_eve = i_sum_eve + n_sum_eve;
        end
        Rb_ins = (log2(2 * M^2) / 2) - (i_sum_bob / (4 * M^2));
        Re_ins = (log2(2 * M^2) / 2) - (i_sum_eve / (4 * M^2));
        Rb_array(iter_index) = Rb_ins;
        Re_array(iter_index) = Re_ins;
        Rs_array(iter_index) = max(0, (Rb_ins - Re_ins));
        % /////////////////////////////////////////////////////////////////
    end
    Rb = mean(Rb_array);
    Re = mean(Re_array);
    Rs = mean(Rs_array);
end

%% INNER FUNCTIONS (TOTAL OF 2)
%==========================================================================
% 1. All Possible Alamouti Matrices

% ARGUMENTS
% 1-) M: Constellation size
% 2-) ss: Constellation set

% OUTPUT
% - X: Box of all possible Alamouti matrices
%==========================================================================
function X = AllPossibleMatrices(M, ss)
    X = zeros(2, 2, M^2);
    n = 1;
    for x1_index = 1 : M
        for x2_index = 1 : M
            x1_current = ss(x1_index);
            x2_current = ss(x2_index);

            X(:, :, n) = [x1_current x2_current;
                          -x2_current' x1_current'];
            n = n + 1;
        end
    end
end
%==========================================================================


%==========================================================================
% 2. Constructing Artificial Noise Matrix

% ARGUMENTS
% 1-) h: Channel of the legitimate receiver
% 2-) A1: Index of the first activated transmit antenna

% OUTPUT
% - Z: Artificial noise matrix
%==========================================================================
function Z = AN_Matrix(h, A1)
    v = (randn(1, 1) + 1i * randn(1, 1)) / sqrt(2);
    
    z11 = -h(3) * v;
    z21 = h(A1) * v;
    z12 = -h(3)' * v;
    z22 = -h(A1)' * v;
    
    Z = zeros(3, 2);
    Z(A1, :) = [z11 -z12'];
    Z(3, :) = [z21 z22'];
end
%==========================================================================