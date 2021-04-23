%==========================================================================
% Burak Özpoyraz, 2020
% Bit Error Rate of CIOD-IM Scheme

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
% 1-) BER_bob: Bit error rate of the legitimate receiver
% 2-) BER_eve: Bit error rate of the eavesdropper
% 3-) error_bob: Number of bit errors of the legitimate receiver
% 4-) error_eve: Number of bit errors of the eavesdropper
%==========================================================================

%% MAIN FUNCTION
function [BER_bob, BER_eve, error_bob, error_eve] = CIOD_IM_BER(num_iterations, N, M, P_tot_des, alpha, sigma2, SNRdB, mod_type)
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
    ss = ss * sqrt(P_symbol_des / Pc); % Consellation power normalization: E[|xk|^2] = P_symbol_des
    
    m = log2(M); % Number of bits for one symbol selection
    n_IM = log2(N) - 1; % Number of bits for index modulation
    k_tot = n_IM + 4 * m; % Number of total bits transmitted during an iteration
    N0 = P_symbol_des / (10^(SNRdB / 10)); % The noise power
    
    num_bits = num_iterations * k_tot; % Number of total bits transmitted during the whole simulation
    bit_array = randi([0, 1], 1, num_bits);

    error_bob = 0;
    error_eve = 0;
    
    %% SIMULATION
    for bit_index = 1 : k_tot : num_bits
        % Transmitter//////////////////////////////////////////////////////
        CIOD_bits = bit_array(bit_index : bit_index + n_IM - 1);
        CIOD_index = Bit2Dec(CIOD_bits);
        
        x_symbol_array = zeros(1, 4);
        x_bit_array = zeros(1, (k_tot - n_IM));
        for k = 1 : 4
            x_bit = bit_array((bit_index + n_IM + (k - 1) * m) : (bit_index + n_IM + k * m - 1));
            x_symbol_array(k) = ss(Bit2Dec(x_bit) + 1);
            for j = 1 : m
                x_bit_array(j + m * (k - 1)) = x_bit(j);
            end
        end
        S = CIOD_Matrix(x_symbol_array, CIOD_index, N);
        % /////////////////////////////////////////////////////////////////
        
        % Channel & Noise//////////////////////////////////////////////////
        h_est = (randn(1, N) + 1i * randn(1, N)) / sqrt(2); % Estimated channel of the legitimate receiver
        h_err = (randn(1, N) + 1i * randn(1, N)) / sqrt(2); % Erroneous channel of the legitimate receiver
        h = sqrt(1 - sigma2) * h_est + sqrt(sigma2) * h_err; % Total channel of the legitimate receiver
        
        g_est = (randn(1, N) + 1i * randn(1, N)) / sqrt(2); % Estimated channel of the eavesdropper
        g_err = (randn(1, N) + 1i * randn(1, N)) / sqrt(2); % Erroneous channel of the eavesdropper
        g = sqrt(1 - sigma2) * g_est + sqrt(sigma2) * g_err; % Total channel of the eavesdropper

        nb = sqrt(N0 / 2) * (randn(1, 4) + 1i * randn(1, 4)); % Noise of the legitimate receiver
        ne = sqrt(N0 / 2) * (randn(1, 4) + 1i * randn(1, 4)); % Noise of the eavesdropper
        % /////////////////////////////////////////////////////////////////
        
        % Artificial Noise/////////////////////////////////////////////////
        Z = AN_Matrix(h_est, CIOD_index, N);
        P_Z = sum(sum(abs(Z).^2));
        Z_N = Z * sqrt(P_Z_des / P_Z); % Artificial noise normalization
        % /////////////////////////////////////////////////////////////////

        % Receiver/////////////////////////////////////////////////////////
        yb = h * (S + Z_N) + nb; % Received signal at the legitimate receiver        
        ye = g * (S + Z_N) + ne; % Received signal at the eavesdropper
        % /////////////////////////////////////////////////////////////////
        
        % Detector/////////////////////////////////////////////////////////
        detected_x_bit_array_bob = zeros(1, (k_tot - n_IM));
        detected_x_bit_array_eve = zeros(1, (k_tot - n_IM));

        error_metric_bob = zeros(4, M, 2^n_IM);
        error_metric_eve = zeros(4, M, 2^n_IM);
        for CIOD_i = 0 : (2^n_IM - 1)
            A = A_Matrices(CIOD_i, N);
            for k = 1 : 4
                A_first = A(:, :, (2 * k - 1));
                A_second = A(:, :, (2 * k));
                for j = 1 : M
                    xj = ss(j);
                    xjI = real(xj);
                    xjQ = imag(xj);
                    error_metric_bob(k, j, (CIOD_i + 1)) = norm(yb - h_est * (A_first * xjI + A_second * xjQ))^2;
                    error_metric_eve(k, j, (CIOD_i + 1)) = norm(ye - g_est * (A_first * xjI + A_second * xjQ))^2;
                end
            end
        end

        detected_CIOD_index_bob = DetectCIOD_Index(error_metric_bob, n_IM);
        detected_CIOD_bit_bob = Dec2Bit(detected_CIOD_index_bob, n_IM);

        detected_CIOD_index_eve = DetectCIOD_Index(error_metric_eve, n_IM);
        detected_CIOD_bit_eve = Dec2Bit(detected_CIOD_index_eve, n_IM);

        for k = 1 : 4
            [~, min_ind_bob] = min(error_metric_bob(k, :, detected_CIOD_index_bob + 1));
            detected_x_bit_array_bob(((k - 1) * m + 1) : k * m) = Dec2Bit(min_ind_bob - 1, m);

            [~, min_ind_eve] = min(error_metric_eve(k, :, detected_CIOD_index_eve + 1));
            detected_x_bit_array_eve(((k - 1) * m + 1) : k * m) = Dec2Bit(min_ind_eve - 1, m);
        end

        error_bob = error_bob + sum(xor(detected_CIOD_bit_bob, CIOD_bits));
        error_bob = error_bob + sum(xor(detected_x_bit_array_bob, x_bit_array));

        error_eve = error_eve + sum(xor(detected_CIOD_bit_eve, CIOD_bits));
        error_eve = error_eve + sum(xor(detected_x_bit_array_eve, x_bit_array));
        % /////////////////////////////////////////////////////////////////
    end
    BER_bob = error_bob / num_bits;
    BER_eve = error_eve / num_bits;
end

%% INNER FUNCTIONS (TOTAL OF 8)
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
% 2. Conversion from Bit to Decimal

% ARGUMENTS
% 1-) bit_array: Bit array to be converted to decimal value

% OUTPUT
% - decimal_value: Corresponding decimal value of the bit array
%==========================================================================
function decimal_value = Bit2Dec(bit_array)
    len = length(bit_array);
    decimal_value = bit_array * (2.^((len - 1) : -1 : 0))';
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
% 1-) h: Channel of the legitimate receiver
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


%==========================================================================
% 6. Constructing Weight Matrices for CIOD Decoding

% ARGUMENTS
% 1-) CIOD_index: Index of the CIOD matrix to be selected from the index
%                 modulation map
% 2-) N: Number of transmit antennas

% OUTPUT
% - A: Weight Matrices Box
%==========================================================================
function A = A_Matrices(CIOD_index, N)
    j = 1i;
    A = zeros(N, 4, 8);
    
    inner_matrices = zeros(2, 2, 4);
    inner_matrices(:, :, 1) = [1 0; 0 1];
    inner_matrices(:, :, 2) = [j 0; 0 -j];
    inner_matrices(:, :, 3) = [0 -1; 1 0];
    inner_matrices(:, :, 4) = [0 j; j 0];
    
    row_range1 = (2 * CIOD_index + 1) : (2 * CIOD_index + 2);
    col_range1 = 1 : 2;

    row_range2 = (N - 2 * CIOD_index - 1) : (N - 2 * CIOD_index);
    col_range2 = 3 : 4;
    for k = 1 : 8
        if k <= 4 && mod(k, 2) == 1
            A(row_range1, col_range1, k) = inner_matrices(:, :, k);
        elseif k <= 4 && mod(k, 2) == 0
            A(row_range2, col_range2, k) = inner_matrices(:, :, k);
        elseif k > 4 && mod(k, 2) == 0
            A(row_range1, col_range1, k) = inner_matrices(:, :, (k - 4));
        elseif k > 4 && mod(k, 2) == 1
            A(row_range2, col_range2, k) = inner_matrices(:, :, (k - 4));
        end
    end
end
%==========================================================================


%==========================================================================
% 7. Detection of the CIOD Index

% ARGUMENTS
% 1-) error_metric: Box of the distances between each received symbol and
%                   each symbol in the constellation for different CIOD
%                   matrices according to different index modulation bits
% 2-) n_IM: Number of bits for index modulation

% OUTPUT
% - detected_CIOD_index: Detected index of CIOD matrix from the map
%==========================================================================
function detected_CIOD_index = DetectCIOD_Index(error_metric, n_IM)
    CIOD_index_total_errors = zeros(1, 2^n_IM);
    for i = 1 : 2^n_IM
        total_error = 0;
        for k = 1 : 4
            min_error = min(error_metric(k, :, i));
            total_error = total_error + min_error;
        end
        CIOD_index_total_errors(i) = total_error;
    end
    
    [~, min_ind] = min(CIOD_index_total_errors);
    detected_CIOD_index = min_ind - 1;
end
%==========================================================================


%==========================================================================
% 8. Conversion from Decimal to Bit

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