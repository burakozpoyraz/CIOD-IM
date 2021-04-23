%==========================================================================
% Burak Özpoyraz, 2020
% Theoretical Bit Error Rate Upper Bound of CIOD-IM Scheme

% ARGUMENTS
% 1-) N: Number of transmit antennas at Alice
% 2-) M: Constellation size
% 3-) P_tot_des: Desired total power to be transmitted during an iteration
% 4-) alpha: Ratio of the power allocated to CIOD matrix
% 5-) SNRdB: Signal-to-noise ratio in dB scale
% 6-) mod_type: Constellation scheme (PSK or QAM)

% OUTPUTS
% - Pb_upper: Theoretical bit error rate upper bound of the legitimate receiver
%==========================================================================

%% MAIN FUNCTION
function Pb_upper = CIOD_IM_TheoBER(N, M, P_tot_des, alpha, SNRdB, mod_type)
    %% PARAMETERS
    P_S_des = P_tot_des * alpha; % Desired power of the CIOD matrix
    P_symbol_des = P_S_des / 8; % Desired power of a single symbol in the constellation

    rotation_angle = FindAngle(M, mod_type); % Rotation angle of CIOD constellation in degree
    rotation_theta = rotation_angle * pi / 180; % Rotation angle of CIOD constellation in radian
    if mod_type == "QAM"
        ss = qammod(0 : M-1, M, "Gray") * exp(1i * rotation_theta); % Constellation set
    else
        ss = pskmod(0 : M-1, M, 0, "Gray") * exp(1i * rotation_theta); % Constellation set
    end
    Pc = ss * ss'/ M; % Average power of a symbol in the constellation
    ss = ss * sqrt(P_symbol_des / Pc); % Consellation power normalization: E[|xk|^2] = 1
    
    m = log2(M); % Number of bits for one symbol selection
    n_IM = log2(N) - 1; % Number of bits for index modulation
    k_tot = n_IM + 4 * m; % Number of total bits transmitted during an iteration
    N0 = P_symbol_des / (10^(SNRdB / 10)); % The noise power
    
    A = N * M^4 / 2; % Number of all possible CIOD matrices
    X = AllPossibleMatrices(N, M, ss); % Set of all CIOD matrices
    
    Pb_upper = 0; % Array of theoretical BER upper bound
    
    %% CALCULATION
    for n = 1 : A
        n_bits = Dec2Bit(n - 1, k_tot);
        Xn = X(:, :, n);
        for u = 1 : A
            fprintf("SNRdB = %d, n = %d, u = %d\n", SNRdB, n, u);
            u_bits = Dec2Bit(u - 1, k_tot);
            Xu = X(:, :, u);

            e_nu = sum(xor(n_bits, u_bits));
            if e_nu == 0
                PEP = 0;
            else
                PEP = PairwiseErrorProbability(Xn, Xu, P_S_des, N0);
            end
            Pb_upper = Pb_upper + (e_nu * PEP / k_tot);
        end
    end
end

%% INNER FUNCTIONS (TOTAL OF 6)
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
% 5. Conversion from Decimal to Bit

% ARGUMENTS
% 1-) decimal: Decimal value to be converted to bit array
% 2-) len: Number of bits that the resulting bit array has

% OUTPUT
% - bit_array: Corresponding bit array of the decimal value
%==========================================================================
function bit_array = Dec2Bit(decimal, len)
    if decimal > (2^len - 1)
        bit_array = -1;
    elseif decimal == 0 || decimal == 1
        bit_array = [zeros(1, len - 1), decimal];
    else
        quotient = [];
        quotient(1) = decimal;
        remainder = [];
        remainder(1) = -1;
        counter = 2;
        while quotient(counter - 1) > 1
            quotient(counter) = floor(quotient(counter - 1) / 2);
            remainder(counter) = mod(quotient(counter - 1), 2);
            counter = counter + 1;
        end
        counter = counter - 1;
        zero_pad_length = len - counter;
        zero_pad = zeros(1, zero_pad_length);

        bit0 = quotient(length(quotient));
        bit_rest = remainder(length(remainder) : -1 : 2);
        bit_array = [bit0, bit_rest];

        bit_array = [zero_pad, bit_array];
    end
end
%==========================================================================


%==========================================================================
% 6. Pairwise Error Probability

% ARGUMENTS
% 1-) Xn: Transmitted CIOD matrix
% 2-) Xu: Erroneously detected CIOD matrix
% 3-) P_S_des: Desired power of the CIOD matrix
% 4-) N0: The noise power
%==========================================================================
function PEP = PairwiseErrorProbability(Xn, Xu, P_S_des, N0)
    Phi_nu = Xn - Xu;
    Delta = Phi_nu' * Phi_nu;
    lambda_array = eig(Delta);
    
    syms theta
    MGF = (1 + ((P_S_des * lambda_array(1)) / (2 * N0 * sin(theta)^2)))^(-1) *...
          (1 + ((P_S_des * lambda_array(2)) / (2 * N0 * sin(theta)^2)))^(-1) *...
          (1 + ((P_S_des * lambda_array(3)) / (2 * N0 * sin(theta)^2)))^(-1) *...
          (1 + ((P_S_des * lambda_array(4)) / (2 * N0 * sin(theta)^2)))^(-1);
    
    PEP = double(int(MGF, [0 (pi / 2)]) / pi);
end
%==========================================================================