%==========================================================================
% Burak Özpoyraz, 2020
% Bit Error Rate of Efficient Alamouti Scheme

% ARGUMENTS
% 1-) num_iterations: Number of iterations for Monte Carlo simulation
% 2-) M: Constellation size
% 3-) P_tot_des: Desired total power to be transmitted during an iteration
% 4-) alpha: Ratio of the power allocated to Alamouti matrix
% 5-) SNRdB: Signal-to-noise ratio in dB scale

% OUTPUTS
% 1-) BER_bob: Bit error rate of the legitimate receiver
% 2-) BER_eve: Bit error rate of the eavesdropper
% 3-) error_bob: Number of bit errors of the legitimate receiver
% 4-) error_eve: Number of bit errors of the eavesdropper
%==========================================================================

%% MAIN FUNCTION
function [BER_bob, BER_eve, error_bob, error_eve] = EfficientAlamoutiBER(num_iterations, M, P_tot_des, alpha, SNRdB)
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
    ss = ss * sqrt(P_symbol_des / Pc); % Consellation power normalization: E[|xk|^2] = 1

    m = log2(M); % Number of bits for one symbol selection
    k_tot = 1 + 2 * m; % Number of total bits transmitted during an iteration
    N0 = P_symbol_des / 10^(SNRdB / 10); % The noise power

    num_bits = num_iterations * k_tot; % Number of total bits transmitted during the whole simulation
    bit_array = randi([0, 1], 1, num_bits);
    
    error_bob = 0;
    error_eve = 0;
    
     %% SIMULATION
    for bit_index = 1 : k_tot : num_bits
        % Transmitter//////////////////////////////////////////////////////
        A1_bit = bit_array(bit_index);
        A1 = A1_bit + 1;

        x1_bits = zeros(1, m);
        for j = 1 : m
            x1_bits(j) = bit_array(bit_index + 1 + j - 1);
        end
        x1 = ss(Bit2Dec(x1_bits) + 1);

        x2_bits = zeros(1, m);
        for j = 1 : m
            x2_bits(j) = bit_array(bit_index + 1 + m + j - 1);
        end
        x2 = ss(Bit2Dec(x2_bits) + 1);
        
        X = zeros(3, 2);
        X(A1, :) = [x1 x2];
        X(3, :) = [-x2' x1'];
        % /////////////////////////////////////////////////////////////////
        
        % Channel & Noise//////////////////////////////////////////////////
        h = (randn(1, 3) + 1i * randn(1, 3)) / sqrt(2); % Channel of the legitimate receiver
        g = (randn(1, 3) + 1i * randn(1, 3)) / sqrt(2); % Channel of the eavesdropper
        
        nb = sqrt(N0 / 2) * (randn(1, 2) + 1i * randn(1, 2)); % Noise of the legitimate receiver
        ne = sqrt(N0 / 2) * (randn(1, 2) + 1i * randn(1, 2)); % Noise of the eavesdropper
        % /////////////////////////////////////////////////////////////////
        
        % Artificial Noise/////////////////////////////////////////////////
        Z = AN_Matrix(h, A1);
        P_Z = sum(sum(abs(Z).^2));
        Z_N = Z * sqrt(P_Z_des / P_Z); % Artificial noise normalization
        % /////////////////////////////////////////////////////////////////
        
        % Receiver/////////////////////////////////////////////////////////
        yb = h * (X + Z_N) + nb; % Received signal at the legitimate receiver
        ye = g * (X + Z_N) + ne; % Received signal at the eavesdropper
        % /////////////////////////////////////////////////////////////////
        
        % Detector/////////////////////////////////////////////////////////
        metric_bob = zeros(2, M, M);
        metric_eve = zeros(2, M, M);
        for A1_index = 1 : 2
            for x1_index = 1 : M
                for x2_index = 1 : M
                    x1_current = ss(x1_index);
                    x2_current = ss(x2_index);
                    
                    X_current = zeros(3, 2);
                    X_current(A1_index, :) = [x1_current x2_current];
                    X_current(3, :) = [-x2_current' x1_current'];
                    
                    metric_bob(A1_index, x1_index, x2_index) = norm(yb - h * X_current)^2;
                    metric_eve(A1_index, x1_index, x2_index) = norm(ye - g * X_current)^2;
                end
            end
        end
        [min_dim1_bob, min_dim2_bob, min_dim3_bob] = Min3D(metric_bob);
        detected_A1_bit_bob = min_dim1_bob - 1;
        detected_x1_bits_bob = Dec2Bit(min_dim2_bob - 1, m);
        detected_x2_bits_bob = Dec2Bit(min_dim3_bob - 1, m);

        [min_dim1_eve, min_dim2_eve, min_dim3_eve] = Min3D(metric_eve);
        detected_A1_bit_eve = min_dim1_eve - 1;
        detected_x1_bits_eve = Dec2Bit(min_dim2_eve - 1, m);
        detected_x2_bits_eve = Dec2Bit(min_dim3_eve - 1, m);

        error_bob = error_bob + sum(xor(detected_A1_bit_bob, A1_bit));
        error_bob = error_bob + sum(xor(detected_x1_bits_bob, x1_bits));
        error_bob = error_bob + sum(xor(detected_x2_bits_bob, x2_bits));

        error_eve = error_eve + sum(xor(detected_A1_bit_eve, A1_bit));
        error_eve = error_eve + sum(xor(detected_x1_bits_eve, x1_bits));
        error_eve = error_eve + sum(xor(detected_x2_bits_eve, x2_bits));
        % /////////////////////////////////////////////////////////////////
    end
    BER_bob = error_bob / num_bits;
    BER_eve = error_eve / num_bits;
end

%% INNER FUNCTIONS (TOTAL OF 5)
%==========================================================================
% 1. Conversion from Bit to Decimal

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
    
    
%==========================================================================
% 3. Minimum of 3D Matrix

% ARGUMENTS
% 1-) matrix3D: 3D matrix

% OUTPUT
% 1-) min_dim1: Dimension-1 index of minimum element
% 2-) min_dim1: Dimension-2 index of minimum element
% 3-) min_dim1: Dimension-3 index of minimum element
%==========================================================================
function [min_dim1, min_dim2, min_dim3] = Min3D(matrix3D)
    dim1 = size(matrix3D, 1);
    dim2 = size(matrix3D, 2);
    
    min_val = min(min(min(matrix3D)));
    index = find(matrix3D == min_val);
    
    [q1, r1] = QuoRem(index, (dim1 * dim2));
    if r1 > 0
        min_dim3 = q1 + 1;
    else
        min_dim3 = q1;
    end

    [q2, r2] = QuoRem(r1, dim1);
    if r2 > 0
    min_dim2 = q2 + 1;
    elseif r2 == 0 && q2 == 0
        min_dim2 = dim2;
    else
        min_dim2 = q2;
    end

    if r2 == 0
        min_dim1 = dim1;
    else
        min_dim1 = r2;
    end

    min_dim1 = double(min_dim1);
    min_dim2 = double(min_dim2);
    min_dim3 = double(min_dim3);
end
%==========================================================================    
    
    
%==========================================================================
% 4. Quotient and Remainder

% ARGUMENTS
% 1-) a: Dividend
% 2-) b: Divider

% OUTPUT
% 1-) q: Quotient
% 2-) r: Remainder
%==========================================================================
function [q, r] = QuoRem(a, b)
    q = floor(a / b);
    r = a - q * b;
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