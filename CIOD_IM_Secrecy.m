%==========================================================================
% Burak Özpoyraz, 2020
% Ergodic Secrecy Rate of CIOD-IM Scheme

% ARGUMENTS
% 1-) num_iterations: Number of iterations for Monte Carlo simulation
% 2-) N: Number of transmit antennas at Alice
% 3-) M: Constellation size
% 4-) P_tot_des: Desired total power to be transmitted during an iteration
% 5-) alpha: Ratio of the power allocated to CIOD matrix
% 6-) sigma2: Power of the erroneous estimation of the channel
% 7-) SNRdB: Signal-to-noise ratio in dB scale
% 8-) mod_type: Constellation scheme (PSK or QAM)

% OUTPUTS
% 1-) Rb: Ergodic rate of the legitimate receiver
% 2-) Re: Ergodic rate of the eavesdropper
% 3-) Rs: Ergodic secrecy rate
%==========================================================================

%% MAIN FUNCTION
function [Rb, Re, Rs] = CIOD_IM_Secrecy(num_iterations, N, M, P_tot_des, alpha, sigma2, SNRdB, mod_type)
    %% PARAMETERS
    P_S_des = P_tot_des * alpha; % Desired power of the CIOD matrix
    P_symbol_des = P_S_des / 8; % Desired power of a single symbol in the constellation
    P_Z_des = P_tot_des * (1 - alpha); % Desired power of the AN matrix

    rotation_angle = FindAngle(M, mod_type); % Rotation angle of CIOD constellation in degree
    rotation_theta = rotation_angle * pi / 180; % Rotation angle of CIOD constellation in radian
    if mod_type == "QAM"
        ss = qammod(0 : M-1, M, "Gray") * exp(1i * rotation_theta); % Constellation set
    else
        ss = pskmod(0 : M-1, M, 0, "Gray") * exp(1i * rotation_theta); % Constellation set
    end
    Pc = ss * ss'/ M; % Average power of a symbol in the constellation
    ss = ss * sqrt(P_symbol_des / Pc);  % Consellation power normalization: E[|xk|^2] = 1
    
    A = N * M^4 / 2; % Number of all possible CIOD matrices
    X = AllPossibleMatrices(N, M, ss); % Set of all CIOD matrices

    N0 = P_symbol_des / (10^(SNRdB / 10)); % The noise power
    
    Rb_array = zeros(1, num_iterations); % Array of ergodic rates of Bob
    Re_array = zeros(1, num_iterations); % Array of ergodic rates of Eve
    Rs_array = zeros(1, num_iterations); % Array of ergodic secrecy rates   

    %% SIMULATION
    for iter_index = 1 : num_iterations
        fprintf("Iteration Index: %d\n", iter_index);
        
        % Channel & Noise//////////////////////////////////////////////////
        h_est = (randn(1, N) + 1i * randn(1, N)) / sqrt(2); % Estimated channel of the legitimate receiver
        h_err = (randn(1, N) + 1i * randn(1, N)) / sqrt(2); % Erroneous channel of the legitimate receiver

        g_est = (randn(1, N) + 1i * randn(1, N)) / sqrt(2); % Estimated channel of the eavesdropper        
        g_err = (randn(1, N) + 1i * randn(1, N)) / sqrt(2); % Erroneous channel of the eavesdropper

        nb = sqrt(N0 / 2) * (randn(1, 4) + 1i * randn(1, 4)); % Noise of the legitimate receiver
        ne = sqrt(N0 / 2) * (randn(1, 4) + 1i * randn(1, 4)); % Noise of the eavesdropper
        % /////////////////////////////////////////////////////////////////

        % Ergodic Rates and Secrecy Rate///////////////////////////////////
        n_sum_bob = 0;
        n_sum_eve = 0;
        for n = 1 : A
            % CIOD Matrix//////////////////////////////////////////////////
            Xn = X(:, :, n);
            % /////////////////////////////////////////////////////////////
            
            n2_sum_bob = 0;
            n2_sum_eve = 0;
            for n2 = 1 : A
                % CIOD Matrix//////////////////////////////////////////////
                Xn2 = X(:, :, n2);
                % /////////////////////////////////////////////////////////
                
                % Artificial Noise/////////////////////////////////////////
                CIOD_index = floor((n - 1) / M^4);
                Z = AN_Matrix(h_est, CIOD_index, N);
                P_Z = sum(sum(abs(Z).^2));
                Z_N = Z * sqrt(P_Z_des / P_Z); % Artificial noise normalization
                % /////////////////////////////////////////////////////////
                
                % Linear Whitening Transformation//////////////////////////////////
                nb_hat = sqrt(sigma2) * h_err * (Xn + Z_N) + nb;
                Nb_hat = nb_hat * nb_hat' / 4;
                Psib = sqrt(N0) * Nb_hat^(-1 / 2);
                h_est_tilde = Psib * sqrt((1 - sigma2)) * h_est;
                nb_tilde = Psib * nb_hat;

                ne_hat = sqrt(1 - sigma2) * g_est * Z_N + sqrt(sigma2) * g_err * (Xn + Z_N) + ne;
                Ne_hat = ne_hat * ne_hat' / 4;
                Psie = sqrt(N0) * (Ne_hat)^(-1 / 2);
                g_est_tilde = Psie * sqrt((1 - sigma2)) * g_est;
                ne_tilde = ne_hat * Psie;
                % /////////////////////////////////////////////////////////////////
                
                value_bob = exp(-(norm(h_est_tilde * (Xn - Xn2) + nb_tilde)^2 - norm(nb_tilde)^2) / N0);
                n2_sum_bob = n2_sum_bob + value_bob;
                
                value_eve = exp(-(norm(g_est_tilde * (Xn - Xn2) + ne_tilde)^2 - norm(ne_tilde)^2) / N0);
                n2_sum_eve = n2_sum_eve + value_eve;
            end
            n_sum_bob = n_sum_bob + log2(n2_sum_bob);
            n_sum_eve = n_sum_eve + log2(n2_sum_eve);
        end
        Rb_ins = (log2(A) - (n_sum_bob / A)) / 4;
        Re_ins = (log2(A) - (n_sum_eve / A)) / 4;
        Rb_array(iter_index) = Rb_ins;
        Re_array(iter_index) = Re_ins;
        Rs_array(iter_index) = max(0, (Rb_ins - Re_ins));
        % /////////////////////////////////////////////////////////////////
    end
    Rb = mean(Rb_array);
    Re = mean(Re_array);
    Rs = mean(Rs_array);
end

%% INNER FUNCTIONS (TOTAL OF 5)
%==========================================================================
% 1. Rotation Angle of CIOD Constellation

% ARGUMENTS
% 1-) M: Constellation size
% 2-) mod_type: Constellation scheme (PSK or QAM)

% OUTPUT
% - rotation_angle: Rotation angle of CIOD in degree
%==========================================================================
function rotation_angle = FindAngle(M, mod_type)
    if mod_type == "QAM"
        if mod(log2(M), 2) == 0 % Square Constellations
            rotation_angle = 31.7175;
        else % Rectangular Constellations
            rotation_angle = 0; % Left 0 to be studied later
        end
    elseif mod_type == "PSK"
        if M == 4 % QPSK
            rotation_angle = 13.2825;
        else
            rotation_angle = 0; % Left 0 to be studied later
        end
    end
end
%==========================================================================


%==========================================================================
% 2. All Possible CIOD Matrices

% ARGUMENTS
% 1-) N: Number of transmit antennas at Alice
% 2-) M: Constellation size
% 3-) ss: Constellation set

% OUTPUT
% - X: Box of all possible CIOD matrices
%==========================================================================
function X = AllPossibleMatrices(N, M, ss)
    A = N * M^4 / 2;
    X = zeros(N, 4, A);
    n = 1;
    for CIOD_index = 0 : (N / 2) - 1
        for symbol1_index = 1 : M
            x1 = ss(symbol1_index);
            for symbol2_index = 1 : M
                x2 = ss(symbol2_index);
                for symbol3_index = 1 : M
                    x3 = ss(symbol3_index);
                    for symbol4_index = 1 : M
                        x4 = ss(symbol4_index);
                        x_symbol_array = [x1 x2 x3 x4];
                        S = CIOD_Matrix(x_symbol_array, CIOD_index, N);
                        X(:, :, n) = S;
                        n = n + 1;
                    end
                end
            end
        end
    end
end
%==========================================================================


%==========================================================================
% 3. Constructing CIOD Matrix 

% ARGUMENTS
% 1-) x_array: Symbol array to be transmitted ([x1 x2 x3 x4])
% 2-) CIOD_index: Index of the CIOD matrix to be selected from the index
%                 modulation map
% 3-) N: Number of transmit antennas

% OUTPUT
% - S: CIOD matrix
%==========================================================================
function S = CIOD_Matrix(x_array, CIOD_index, N)
    x1_tilde = real(x_array(1)) + 1i * imag(x_array(3));
    x2_tilde = real(x_array(2)) + 1i * imag(x_array(4));
    x3_tilde = real(x_array(3)) + 1i * imag(x_array(1));
    x4_tilde = real(x_array(4)) + 1i * imag(x_array(2));
    
    Theta1 = Theta_Matrix([x1_tilde x2_tilde]);
    Theta2 = Theta_Matrix([x3_tilde x4_tilde]);
    
    S = zeros(N, 4);
    
    row_range1 = (2 * CIOD_index + 1) : (2 * CIOD_index + 2);
    col_range1 = 1 : 2;
    S(row_range1, col_range1) = Theta1;

    row_range2 = (N - 2 * CIOD_index - 1) : (N - 2 * CIOD_index);
    col_range2 = 3 : 4;
    S(row_range2, col_range2) = Theta2;
end
%==========================================================================


%==========================================================================
% 4. Constructing Theta Matrix

% ARGUMENTS
% 1-) x_array: Symbol array to be transmitted ([x1 x2 x3 x4])

% OUTPUT
% - Theta: Theta matrix
%==========================================================================
function Theta = Theta_Matrix(x_array)
    x1 = x_array(1);
    x2 = x_array(2);
    Theta = [x1 -x2';
             x2 x1'];
end
%==========================================================================


%==========================================================================
% 5. Constructing Artificial Noise Matrix

% ARGUMENTS
% 1-) h: Channel vector of the legitimate receiver
% 2-) CIOD_index: Index of the CIOD matrix to be selected from the index
%                 modulation map
% 3-) N: Number of transmit antennas

% OUTPUT
% - Z: Artificial noise matrix
%==========================================================================
function Z = AN_Matrix(h, CIOD_index, N)
    v = (randn(1, 1) + 1i * randn(1, 1)) / sqrt(2);
    
    Z = zeros(N, 4);
    
    z11 = -h(2 * CIOD_index + 2) * v;
    z21 = h(2 * CIOD_index + 1) * v;
    z12 = -h(N - 2 * CIOD_index) * v;
    z22 = h(N - 2 * CIOD_index - 1) * v;
    
    Q1 = [z11 z11;
          z21 z21];
    Q2 = [z12 z12;
          z22 z22];

    row_range1 = (2 * CIOD_index + 1) : (2 * CIOD_index + 2);
    col_range1 = 1 : 2;
    Z(row_range1, col_range1) = Q1;

    row_range2 = (N - 2 * CIOD_index - 1) : (N - 2 * CIOD_index);
    col_range2 = 3 : 4;
    Z(row_range2, col_range2) = Q2;
end
%==========================================================================