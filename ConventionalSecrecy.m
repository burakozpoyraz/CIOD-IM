%==========================================================================
% Burak Özpoyraz, 2020
% Ergodic Secrecy Rate of Conventional Scheme

% ARGUMENTS
% 1-) num_iterations: Number of iterations for Monte Carlo simulation
% 2-) Na: Number of transmit antennas at Alice
% 3-) Nt: Number of selected transmit antennas
% 4-) M: Constellation size
% 5-) P_tot_des: Desired total power to be transmitted during an iteration
% 6-) theta: Ratio of the power allocated to information symbols
% 7-) SNRdB: Signal-to-noise ratio in dB scale
% 8-) TAS_type: Transmit antenna selection scheme ("SLNR" or "Random")

% OUTPUTS
% 1-) Rb: Ergodic rate of the legitimate receiver
% 2-) Re: Ergodic rate of the eavesdropper
% 3-) Rs: Ergodic secrecy rate
%==========================================================================

%% MAIN FUNCTION
function [Rb, Re, Rs] = ConventionalSecrecy(num_iterations, Na, Nt, M, P_tot_des, alpha, SNRdB, TAS_type)
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
    
    N0 = P_x_des / 10^(SNRdB / 10); % The noise power
    
    Rb_array = zeros(1, num_iterations); % Array of Ergodic Rates of Bob
    Re_array = zeros(1, num_iterations); % Array of Ergodic Rates of Eve
    Rs_array = zeros(1, num_iterations); % Array of Ergodic Secrecy Rates

    %% SIMULATION
    for iter_index = 1 : num_iterations
        fprintf("Iteration Index: %d\n", iter_index);
        
        % Channel & Noise//////////////////////////////////////////////////
        h = (randn(1, Na) + 1i * randn(1, Na)) / sqrt(2); % Channel of the legitimate receiver
        g = (randn(1, Na) + 1i * randn(1, Na)) / sqrt(2); % Channel of the eavesdropper
        
        nb = sqrt(N0 / 2) * (randn(1, 1) + 1i * randn(1, 1)); % Noise of the legitimate receiver
        ne = sqrt(N0 / 2) * (randn(1, 1) + 1i * randn(1, 1)); % Noise of the eavesdropper 
        % /////////////////////////////////////////////////////////////////
        
        % Transmit Antenna Selection///////////////////////////////////////
        [Tk, Tq] = TAS(TAS_type, Na, Nt, h, g, N0);
        
        h_Tk = h * Tk;
        h_Tq = h * Tq;
        
        g_Tk = g * Tk;
        g_Tq = g * Tq;
        % /////////////////////////////////////////////////////////////////
        
        % Ergodic Rates and Secrecy Rate///////////////////////////////////
        j_sum_bob = 0;
        j_sum_eve = 0;
        for j = 1 : Nt
            hj = h_Tk(:, j);
            gj = g_Tk(:, j);
            
            m_sum_bob = 0;
            m_sum_eve = 0;
            for m = 1 : M
                xm = ss(m);
                
                j2_sum_bob = 0;
                j2_sum_eve = 0;
                for j2 = 1 : Nt
                    hj2 = h_Tk(:, j2);
                    gj2 = g_Tk(:, j2);
                    
                    m2_sum_bob = 0;
                    m2_sum_eve = 0;
                    for m2 = 1 : M
                        xm2 = ss(m2);
                        
                        % Artificial Noise/////////////////////////////////
                        v = (randn(1, 1) + 1i * randn(1, 1)) / sqrt(2);
                        z1 = h_Tq * v;
                        z2 = h_Tk(:, j) * v;
                        z = [z1 z2];
                        P_z = sum(abs(z).^2);
                        z1_N = z1 * sqrt(P_z_des / P_z); % Artificial noise normalization
                        z2_N = z2 * sqrt(P_z_des / P_z); % Artificial noise normalization
                        % /////////////////////////////////////////////////
                        
                        % Coloured Noise & Whitening Transformation////////
                        ne_hat = gj * z1_N + g_Tq * z2_N + ne;
                        We_hat = ne_hat * ne_hat';
                        Q = sqrt(N0) * (We_hat)^(-1 / 2);
                        
                        gj_tilde = Q * gj;
                        gj2_tilde = Q * gj2;
                        ne_tilde = Q * ne_hat;
                        % /////////////////////////////////////////////////
                        
                        value_bob = exp(-(norm(hj * xm - hj2 * xm2 + nb)^2 - norm(nb)^2) / N0);
                        m2_sum_bob = m2_sum_bob + value_bob;
                        
                        value_eve = exp(-(norm(gj_tilde * xm - gj2_tilde * xm2 + ne_tilde)^2 - norm(ne_tilde)^2) / N0);
                        m2_sum_eve = m2_sum_eve + value_eve;
                    end
                    j2_sum_bob = j2_sum_bob + m2_sum_bob;
                    j2_sum_eve = j2_sum_eve + m2_sum_eve;
                end
                m_sum_bob = m_sum_bob + log2(j2_sum_bob);
                m_sum_eve = m_sum_eve + log2(j2_sum_eve);
            end
            j_sum_bob = j_sum_bob + m_sum_bob;
            j_sum_eve = j_sum_eve + m_sum_eve;
        end
        Rb_ins = log2(Nt * M) - (j_sum_bob / (Nt * M));
        Re_ins = log2(Nt * M) - (j_sum_eve / (Nt * M));
        Rb_array(iter_index) = Rb_ins;
        Re_array(iter_index) = Re_ins;
        Rs_array(iter_index) = max(0, (Rb_ins - Re_ins));
        % /////////////////////////////////////////////////////////////////
    end
    Rb = mean(Rb_array);
    Re = mean(Re_array);
    Rs = mean(Rs_array);
end
%% INNER FUNCTIONS (TOTAL OF 1)
%==========================================================================
% 1. Transmit Antenna Selection

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