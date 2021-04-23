%==========================================================================
% Burak Özpoyraz, 2020
% Bit Error Rate of Conventional Scheme

% ARGUMENTS
% 1-) num_iterations: Number of iterations for Monte Carlo simulation
% 2-) Na: Number of transmit antennas at Alice
% 3-) Nt: Number of selected transmit antennas at Alice
% 4-) M: Constellation size
% 5-) P_tot_des: Desired total power to be transmitted during an iteration
% 6-) alpha: Ratio of the power allocated to information symbols
% 7-) SNRdB: Signal-to-noise ratio in dB scale
% 8-) TAS_type: Transmit antenna selection scheme ("SLNR" or "Random")

% OUTPUTS
% 1-) BER_bob: Bit error rate of the legitimate receiver
% 2-) BER_eve: Bit error rate of the eavesdropper
% 3-) error_bob: Number of bit errors of the legitimate receiver
% 4-) error_eve: Number of bit errors of the eavesdropper
%==========================================================================

%% MAIN FUNCTION
function [BER_bob, BER_eve, error_bob, error_eve] = ConventionalBER(num_iterations, Na, Nt, M, P_tot_des, alpha, SNRdB, TAS_type)
    %% PARAMETERS
    P_x_des = P_tot_des * alpha; % Desired power of the information symbol
    P_z_des = P_tot_des * (1 - alpha); % Desired power of the AN matrix

    if M > 8
        ss = qammod(0 : M-1, M, 'Gray');
    else
        ss = pskmod(0 : M-1, M, 0, 'Gray');
    end
    Pc = ss * ss' / M; % Average power of a symbol in the constellation
    ss = ss * sqrt(P_x_des / Pc); % Consellation power normalization: E[|xk|^2] = 1
    
    nt = log2(Nt); % Number of bits for spatial modulation
    m = log2(M); % Number of bits for the information symbol
    k_tot = nt + m; % Number of total bits transmitted during an iteration
    N0 = P_x_des / 10^(SNRdB / 10); % The noise power

    num_bits = num_iterations * k_tot; % Number of total bits transmitted during the whole simulation
    bit_array = randi([0, 1], 1, num_bits);
    
    error_bob = 0;
    error_eve = 0;
    
    %% SIMULATION
    for bit_index = 1 : k_tot : num_bits
        % Transmitter//////////////////////////////////////////////////////
        index_bits = zeros(1, nt);
        for j = 1 : nt
            index_bits(j) = bit_array(bit_index + j - 1);
        end
        antenna_index = Bit2Dec(index_bits) + 1;

        x_bits = zeros(1, m);
        for j = 1 : m
            x_bits(j) = bit_array(bit_index + nt + j - 1);
        end
        x = ss(Bit2Dec(x_bits) + 1);

        ej = zeros(Nt, 1);
        ej(antenna_index) = 1;
        % /////////////////////////////////////////////////////////////////
        
        % Channel & Noise//////////////////////////////////////////////////
        h = (randn(1, Na) + 1i * randn(1, Na)) / sqrt(2); % Channel of the legitimate receiver
        g = (randn(1, Na) + 1i * randn(1, Na)) / sqrt(2); % Channel of the eavesdropper
        
        nb = sqrt(N0 / 2) * (randn(1, 1) + 1i * randn(1, 1)); % Noise of the legitimate receiver
        ne = sqrt(N0 / 2) * (randn(1, 1) + 1i * randn(1, 1)); % Noise of the eavesdropper
        % /////////////////////////////////////////////////////////////////
        
        % Transmit Antenna Selection///////////////////////////////////////
        [Tk, Tq] = TAS(TAS_type, Na, Nt, h, g, N0);
        % /////////////////////////////////////////////////////////////////
        
        % Artificial Noise/////////////////////////////////////////////////
        v = (randn(1, 1) + 1i * randn(1, 1)) / sqrt(2);
        
        z1 = -h * Tq * v;
        z2 = h * Tk * ej * v;
        
        z = [z1 z2];
        P_z = sum(abs(z).^2);
        z1_N = z1 * sqrt(P_z_des / P_z); % Artificial noise normalization
        z2_N = z2 * sqrt(P_z_des / P_z); % Artificial noise normalization
        % /////////////////////////////////////////////////////////////////
        
        % Receiver/////////////////////////////////////////////////////////
        yb = h * Tk * ej * (x + z1_N) + h * Tq * z2_N + nb;
        ye = g * Tk * ej * (x + z1_N) + g * Tq * z2_N + ne;
        % /////////////////////////////////////////////////////////////////
        
        % Detector/////////////////////////////////////////////////////////
        metric_bob = zeros(Nt, M);
        metric_eve = zeros(Nt, M);
        for j = 1 : Nt
            h_Tk = h * Tk;
            g_Tk = g * Tk;
            h_prime = h_Tk(j);
            g_prime = g_Tk(j);
            for k = 1 : M
                xk = ss(k);
                metric_bob(j, k) = abs(yb - h_prime * xk)^2;
                metric_eve(j, k) = abs(ye - g_prime * xk)^2;
            end
        end
        [min_of_every_col_bob, min_index_of_every_col_bob] = min(metric_bob);
        [~, col_index_of_min_value_bob] = min(min_of_every_col_bob);
        row_index_of_min_value_bob = min_index_of_every_col_bob(col_index_of_min_value_bob);
        detected_index_bits_bob = Dec2Bit(row_index_of_min_value_bob - 1, nt);
        detected_x_bits_bob = Dec2Bit(col_index_of_min_value_bob - 1, m);

        [min_of_every_col_eve, min_index_of_every_col_eve] = min(metric_eve);
        [~, col_index_of_min_value_eve] = min(min_of_every_col_eve);
        row_index_of_min_value_eve = min_index_of_every_col_eve(col_index_of_min_value_eve);
        detected_index_bits_eve = Dec2Bit(row_index_of_min_value_eve - 1, nt);
        detected_x_bits_eve = Dec2Bit(col_index_of_min_value_eve - 1, m);

        error_bob = error_bob + sum(xor(detected_index_bits_bob, index_bits));
        error_bob = error_bob + sum(xor(detected_x_bits_bob, x_bits));

        error_eve = error_eve + sum(xor(detected_index_bits_eve, index_bits));
        error_eve = error_eve + sum(xor(detected_x_bits_eve, x_bits));
        % /////////////////////////////////////////////////////////////////  
    end
    BER_bob = error_bob / num_bits;
    BER_eve = error_eve / num_bits;
end

%% INNER FUNCTIONS (TOTAL OF 3)
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
% 2. Transmit Antenna Selection

% ARGUMENTS
% 1-) TAS_type: Transmit antenna selection scheme ("SLNR" or "Random")
% 2-) Na: Number of transmit antennas at Alice
% 3-) Nt: Number of selected transmit antennas at Alice
% 4-) h: Channel of the legitimate receiver
% 5-) g: Channel of the eavesdropper
% 6-) N0: % The noise power

% OUTPUT
% - decimal_value: Corresponding decimal value of the bit array
%==========================================================================
function [Tk, Tq] = TAS(TAS_type, Na, Nt, h, g, N0)
    if TAS_type == "SLNR"
        SLNR_array = zeros(1, Na);
        for j = 1 : Na
            SLNR = abs(h(j))^2 / (abs(g(j))^2 + N0);
            SLNR_array(j) = SLNR;
        end
        [~, descending_SLNR_index_array] = sort(SLNR_array, 'descend');
        selected_antenna_index_array = descending_SLNR_index_array(1 : Nt);
        non_selected_antenna_index_array = descending_SLNR_index_array((Nt + 1) : end);

        sorted_selected_antenna_index_array = sort(selected_antenna_index_array);        

        I_Na = eye(Na);
        Tk = zeros(Na, Nt);
        for j = 1 : Nt
            index = sorted_selected_antenna_index_array(j);
            eye_col = I_Na(:, index);
            for row = 1 : Na
                Tk(row, j) = eye_col(row);
            end
        end

        random_antenna_index_of_rest = non_selected_antenna_index_array(randi([1, (Na - Nt)], 1, 1));
        Tq = I_Na(:, random_antenna_index_of_rest);
    elseif TAS_type == "Random"
        random_array = randperm(Na);

        selected_antenna_index_array = random_array(1 : Nt);
        non_selected_antenna_index_array = random_array((Nt + 1) : end);

        sorted_selected_antenna_index_array = sort(selected_antenna_index_array);        

        I_Na = eye(Na);
        Tk = zeros(Na, Nt);
        for j = 1 : Nt
            index = sorted_selected_antenna_index_array(j);
            eye_col = I_Na(:, index);
            for row = 1 : Na
                Tk(row, j) = eye_col(row);
            end
        end

        random_antenna_index_of_rest = non_selected_antenna_index_array(randi([1, (Na - Nt)], 1, 1));
        Tq = I_Na(:, random_antenna_index_of_rest);
    end
end
%==========================================================================


%==========================================================================
% 3. Conversion from Decimal to Bit

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