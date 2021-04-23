%==========================================================================
% Burak Özpoyraz, 2020
% Figures of "Index Modulation Based CIOD for Secure Communications" via
% Parallel Cores

% The rows of the BER matrices includes the following values:
% 1. Row: SNRdB Array
% 2. Row: Number of Iterations
% 3. Row: Bit Error Rate of Bob
% 4. Row: Bit Error Rate of Eve
% 5. Row: Number of Errors of Bob
% 6. Row: Number of Errors of Eve

% The rows of the secrecy matrices includes the following values:
% 1. Row: SNRdB Array
% 2. Row: Ergodic Rate of Bob
% 3. Row: Ergodic Rate of Eve
% 4. Row: Ergodic Secrecy Rate
%==========================================================================

%% PARALLELIZATION PARAMETERS
num_workers = 25;

%% COLOR CODES
dark_red = [0.6350 0.0780 0.1840];
dark_blue = [0 0.4470 0.7410];
dark_green = [0.4660 0.6740 0.1880];
dark_yellow = [0.9290 0.6940 0.1250];

%% FIGURE - 3
SNRdB_array = (0 : 5 : 30);
len_SNRdB_array = length(SNRdB_array);

P_tot_des = 1;
alpha = 0.5;
sigma2 = 0;

CIOD_N4_M4 = zeros(6, len_SNRdB_array);
CIOD_N4_M4(1, :) = SNRdB_array;
CIOD_N4_M4(2, :) = [7e2 2e3 2e4 5e5 5e6 5e7 7e2] / num_workers;

CIOD_N32_M4 = zeros(6, len_SNRdB_array);
CIOD_N32_M4(1, :) = SNRdB_array;
CIOD_N32_M4(2, :) = [7e2 7e2 7e3 3e5 6e6 2e7 3e5] / num_workers;

CIOD_N4_M16 = zeros(6, len_SNRdB_array);
CIOD_N4_M16(1, :) = SNRdB_array;
CIOD_N4_M16(2, :) = [7e2 7e2 7e2 7e3 1e5 3e6 4e7] / num_workers;

CIOD_N32_M16 = zeros(6, len_SNRdB_array);
CIOD_N32_M16(1, :) = SNRdB_array;
CIOD_N32_M16(2, :) = [1e2 2e2 4e2 1e3 5e4 9e5 7e6] / num_workers;

CIOD_N4_M4_theo_BER = zeros(1, len_SNRdB_array);
CIOD_N32_M4_theo_BER = zeros(1, len_SNRdB_array);

for SNRdB_index = 1 : len_SNRdB_array
    SNRdB = SNRdB_array(SNRdB_index);
    fprintf("SNR = %d dB has just started. \n", SNRdB);
    
    % CIOD - N=4, QPSK////////////////////////////////////////////////////
    fprintf("CIOD - N=4, QPSK \n");
    num_iterations = CIOD_N4_M4(2, SNRdB_index);
    N = 4;
    M = 4;
    mod_type = "PSK";
    parfor parfor_index = 1 : num_workers
        [CIOD_N4_M4_BER_Bob_parfor(parfor_index),...
         CIOD_N4_M4_BER_Eve_parfor(parfor_index),...
         CIOD_N4_M4_num_error_Bob_parfor(parfor_index),...
         CIOD_N4_M4_num_error_Eve_parfor(parfor_index)] =...
        CIOD_IM_BER(num_iterations, N, M, P_tot_des, alpha, sigma2, SNRdB, mod_type);
    end
    
    CIOD_N4_M4(3, SNRdB_index) = mean(CIOD_N4_M4_BER_Bob_parfor);
    CIOD_N4_M4(4, SNRdB_index) = mean(CIOD_N4_M4_BER_Eve_parfor);
    CIOD_N4_M4(5, SNRdB_index) = sum(CIOD_N4_M4_num_error_Bob_parfor);
    CIOD_N4_M4(6, SNRdB_index) = sum(CIOD_N4_M4_num_error_Eve_parfor);

    CIOD_N4_M4_theo_BER(SNRdB_index) = CIOD_IM_TheoBER(N, M, P_tot_des, alpha, SNRdB, mod_type);
    % /////////////////////////////////////////////////////////////////////
    
    % CIOD - N=32, QPSK////////////////////////////////////////////////////
    fprintf("CIOD - N=32, QPSK \n");
    num_iterations = CIOD_N32_M4(2, SNRdB_index);
    N = 32;
    M = 4;
    mod_type = "PSK";
    parfor parfor_index = 1 : num_workers
        [CIOD_N32_M4_BER_Bob_parfor(parfor_index),...
         CIOD_N32_M4_BER_Eve_parfor(parfor_index),...
         CIOD_N32_M4_num_error_Bob_parfor(parfor_index),...
         CIOD_N32_M4_num_error_Eve_parfor(parfor_index)] =...
        CIOD_IM_BER(num_iterations, N, M, P_tot_des, alpha, sigma2, SNRdB, mod_type);
    end
    
    CIOD_N32_M4(3, SNRdB_index) = mean(CIOD_N32_M4_BER_Bob_parfor);
    CIOD_N32_M4(4, SNRdB_index) = mean(CIOD_N32_M4_BER_Eve_parfor);
    CIOD_N32_M4(5, SNRdB_index) = sum(CIOD_N32_M4_num_error_Bob_parfor);
    CIOD_N32_M4(6, SNRdB_index) = sum(CIOD_N32_M4_num_error_Eve_parfor);

    CIOD_N32_M4_theo_BER(SNRdB_index) = CIOD_IM_TheoBER(N, M, P_tot_des, alpha, SNRdB, mod_type);
    % /////////////////////////////////////////////////////////////////////
    
    % CIOD - N=4, 16-QAM///////////////////////////////////////////////////
    fprintf("CIOD - N=4, 16-QAM \n");
    num_iterations = CIOD_N4_M16(2, SNRdB_index);
    N = 4;
    M = 16;
    mod_type = "QAM";
    parfor parfor_index = 1 : num_workers
        [CIOD_N4_M16_BER_Bob_parfor(parfor_index),...
         CIOD_N4_M16_BER_Eve_parfor(parfor_index),...
         CIOD_N4_M16_num_error_Bob_parfor(parfor_index),...
         CIOD_N4_M16_num_error_Eve_parfor(parfor_index)] =...
        CIOD_IM_BER(num_iterations, N, M, P_tot_des, alpha, sigma2, SNRdB, mod_type);
    end
    
    CIOD_N4_M16(3, SNRdB_index) = mean(CIOD_N4_M16_BER_Bob_parfor);
    CIOD_N4_M16(4, SNRdB_index) = mean(CIOD_N4_M16_BER_Eve_parfor);
    CIOD_N4_M16(5, SNRdB_index) = sum(CIOD_N4_M16_num_error_Bob_parfor);
    CIOD_N4_M16(6, SNRdB_index) = sum(CIOD_N4_M16_num_error_Eve_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    % CIOD - N=32, 16-QAM///////////////////////////////////////////////////
    fprintf("CIOD - N=32, 16-QAM \n");
    num_iterations = CIOD_N4_M16(2, SNRdB_index);
    N = 32;
    M = 16;
    mod_type = "QAM";
    parfor parfor_index = 1 : num_workers
        [CIOD_N32_M16_BER_Bob_parfor(parfor_index),...
         CIOD_N32_M16_BER_Eve_parfor(parfor_index),...
         CIOD_N32_M16_num_error_Bob_parfor(parfor_index),...
         CIOD_N32_M16_num_error_Eve_parfor(parfor_index)] =...
        CIOD_IM_BER(num_iterations, N, M, P_tot_des, alpha, sigma2, SNRdB, mod_type);
    end
    
    CIOD_N32_M16(3, SNRdB_index) = mean(CIOD_N32_M16_BER_Bob_parfor);
    CIOD_N32_M16(4, SNRdB_index) = mean(CIOD_N32_M16_BER_Eve_parfor);
    CIOD_N32_M16(5, SNRdB_index) = sum(CIOD_N32_M16_num_error_Bob_parfor);
    CIOD_N32_M16(6, SNRdB_index) = sum(CIOD_N32_M16_num_error_Eve_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    fprintf("SNR = %d dB has just finished. \n\n", SNRdB);
end

fig3 = figure;
semilogy(SNRdB_array, CIOD_N4_M4(3, :), "d-", "Color", dark_red,...
                                              "MarkerEdgeColor", dark_red, ...
                                              "MarkerFaceColor", "r", ...
                                              "MarkerSize", 8);
hold on;
semilogy(SNRdB_array, CIOD_N32_M4(3, :), "o-", "Color", dark_yellow,...
                                               "MarkerEdgeColor", dark_yellow, ...
                                               "MarkerFaceColor", "y", ...
                                               "MarkerSize", 8);
hold on;
semilogy(SNRdB_array, CIOD_N4_M16(3, :), "s-", "Color", dark_blue,...
                                               "MarkerEdgeColor", dark_blue, ...
                                               "MarkerFaceColor", "b", ...
                                               "MarkerSize", 8);
hold on;
semilogy(SNRdB_array, CIOD_N32_M16(3, :), "^-", "Color", dark_green,...
                                                "MarkerEdgeColor", dark_green, ...
                                                "MarkerFaceColor", "g", ...
                                                "MarkerSize", 8);
hold on;
semilogy(SNRdB_array, CIOD_N4_M4_theo_BER, "k--", "LineWidth", 1.5);
hold on;
semilogy(SNRdB_array, CIOD_N32_M4_theo_BER, "k--", "LineWidth", 1.5);
hold on;
semilogy(SNRdB_array, CIOD_N4_M4(4, :), "d--", "Color", dark_red);
hold on;
semilogy(SNRdB_array, CIOD_N32_M4(4, :), "o--", "Color", dark_yellow);
hold on;
semilogy(SNRdB_array, CIOD_N4_M16(4, :), "s--", "Color", dark_blue);
hold on;
semilogy(SNRdB_array, CIOD_N32_M16(4, :), "^--", "Color", dark_green);

xlabel("$E_s / N_0$", "Interpreter", "latex");
ylabel("BER", "Interpreter", "latex");
legend("CIOD-IM, ($N=4$, QPSK, $l=2.25$)",...
       "CIOD-IM, ($N=32$, QPSK, $l=3$)",...
       "CIOD-IM, ($N=4$, 16-QAM, $l=4.25$)",...
       "CIOD-IM, ($N=32$, 16-QAM, $l=5$)",...
       "CIOD-IM Theoretical Upper Bound",...
       "Location", "northwest", "FontSize", 11, "Interpreter", "latex");
ylim([1e-5 100]);
grid;

annotation("rectangle", [0.735, 0.605, 0.08, 0.07]);
annotation("textbox", [0.725, 0.505, 0.1, 0.1], "String", "Eve",...
                                                "FontSize", 10,...
                                                "Interpreter", "latex",...
                                                "HorizontalAlignment", "center",...
                                                "LineStyle", "none");
annotation("rectangle", [0.39, 0.37, 0.27, 0.04]);
annotation("textbox", [0.637, 0.315, 0.1, 0.1], "String", "Bob",...
                                                "FontSize", 10,...
                                                "Interpreter", "latex",...
                                                "HorizontalAlignment", "center",...
                                                "LineStyle", "none");
% Save As PDF
set(fig3, "Units", "Inches");
pos = get(fig3, "Position");
set(fig3, "PaperPositionMode", "Auto", "PaperUnits", "Inches", "PaperSize", [pos(3), pos(4)])
print(fig3, "Figure3", "-dpdf", "-r0")

%% FIGURE - 4
SNRdB_array = (0 : 3 : 30);
len_SNRdB_array = length(SNRdB_array);

CIOD_N8_M4 = zeros(6, len_SNRdB_array);
CIOD_N8_M4(1, :) = SNRdB_array;
CIOD_N8_M4(2, :) = [5e2 8e2 3e3 8e3 4e4 3e5 2e6 1e7 8e7 1e3 1e3] / num_workers;

eff_ala_M4 = zeros(6, len_SNRdB_array);
eff_ala_M4(1, :) = SNRdB_array;
eff_ala_M4(2, :) = [2e3 4e3 4e3 8e3 4e4 6e4 2e5 8e5 2e6 5e6 9e6] / num_workers;

conv_Na3_Nt2 = zeros(6, len_SNRdB_array);
conv_Na3_Nt2(1, :) = SNRdB_array;
conv_Na3_Nt2(2, :) = [2e3 4e3 5e3 8e3 3e4 4e4 5e4 1e5 4e5 6e5 9e5] / num_workers;

conv_Na6_Nt4 = zeros(6, len_SNRdB_array);
conv_Na6_Nt4(1, :) = SNRdB_array;
conv_Na6_Nt4(2, :) = [2e3 2e3 3e3 4e3 5e3 8e3 3e4 4e4 7e4 2e5 4e5] / num_workers;

for SNRdB_index = 1 : len_SNRdB_array
    SNRdB = SNRdB_array(SNRdB_index);
    fprintf("SNR = %d dB has just started. \n", SNRdB);
    
    % Conventional - Na=3, Nt=2////////////////////////////////////////////
    fprintf("Conventional - Na=3, Nt=2 \n");
    num_iterations = conv_Na3_Nt2(2, SNRdB_index);
    Na = 3;
    Nt = 2;
    M = 2;
    P_tot_des = 1;
    alpha = 0.5;
    TAS_type = "SLNR";
    parfor parfor_index = 1 : num_workers
        [conv_Na3_Nt2_BER_Bob_parfor(parfor_index),...
         conv_Na3_Nt2_BER_Eve_parfor(parfor_index),...
         conv_Na3_Nt2_num_error_Bob_parfor(parfor_index),...
         conv_Na3_Nt2_num_error_Eve_parfor(parfor_index)] =...
        ConventionalBER(num_iterations, Na, Nt, M, P_tot_des, alpha, SNRdB, TAS_type);
    end
    
    conv_Na3_Nt2(3, SNRdB_index) = mean(conv_Na3_Nt2_BER_Bob_parfor);
    conv_Na3_Nt2(4, SNRdB_index) = mean(conv_Na3_Nt2_BER_Eve_parfor);
    conv_Na3_Nt2(5, SNRdB_index) = sum(conv_Na3_Nt2_num_error_Bob_parfor);
    conv_Na3_Nt2(6, SNRdB_index) = sum(conv_Na3_Nt2_num_error_Eve_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    % Conventional - Na=6, Nt=4////////////////////////////////////////////
    fprintf("Conventional - Na=6, Nt=4 \n");
    num_iterations = conv_Na6_Nt4(2, SNRdB_index);
    Na = 6;
    Nt = 4;
    M = 2;
    P_tot_des = 1;
    alpha = 0.5;
    TAS_type = "SLNR";
    parfor parfor_index = 1 : num_workers
        [conv_Na6_Nt4_BER_Bob_parfor(parfor_index),...
         conv_Na6_Nt4_BER_Eve_parfor(parfor_index),...
         conv_Na6_Nt4_num_error_Bob_parfor(parfor_index),...
         conv_Na6_Nt4_num_error_Eve_parfor(parfor_index)] =...
        ConventionalBER(num_iterations, Na, Nt, M, P_tot_des, alpha, SNRdB, TAS_type);
    end
    
    conv_Na6_Nt4(3, SNRdB_index) = mean(conv_Na6_Nt4_BER_Bob_parfor);
    conv_Na6_Nt4(4, SNRdB_index) = mean(conv_Na6_Nt4_BER_Eve_parfor);
    conv_Na6_Nt4(5, SNRdB_index) = sum(conv_Na6_Nt4_num_error_Bob_parfor);
    conv_Na6_Nt4(6, SNRdB_index) = sum(conv_Na6_Nt4_num_error_Eve_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    % Efficient Alamouti///////////////////////////////////////////////////
    fprintf("Efficient Alamouti \n");
    num_iterations = eff_ala_M4(2, SNRdB_index);
    M = 4;
    P_tot_des = 1;
    alpha = 0.5;
    parfor parfor_index = 1 : num_workers
        [eff_ala_M4_BER_Bob_parfor(parfor_index),...
         eff_ala_M4_BER_Eve_parfor(parfor_index),...
         eff_ala_M4_num_error_Bob_parfor(parfor_index),...
         eff_ala_M4_num_error_Eve_parfor(parfor_index)] =...
        EfficientAlamoutiBER(num_iterations, M, P_tot_des, alpha, SNRdB);
    end
    
    eff_ala_M4(3, SNRdB_index) = mean(eff_ala_M4_BER_Bob_parfor);
    eff_ala_M4(4, SNRdB_index) = mean(eff_ala_M4_BER_Eve_parfor);
    eff_ala_M4(5, SNRdB_index) = sum(eff_ala_M4_num_error_Bob_parfor);
    eff_ala_M4(6, SNRdB_index) = sum(eff_ala_M4_num_error_Eve_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    % CIOD - N=8///////////////////////////////////////////////////////////
    fprintf("CIOD - N=8 \n");
    num_iterations = CIOD_N8_M4(2, SNRdB_index);
    N = 8;
    M = 4;
    P_tot_des = 1;
    alpha = 0.5;
    sigma2 = 0;
    mod_type = "PSK";
    parfor parfor_index = 1 : num_workers
        [CIOD_N8_M4_BER_Bob_parfor(parfor_index),...
         CIOD_N8_M4_BER_Eve_parfor(parfor_index),...
         CIOD_N8_M4_num_error_Bob_parfor(parfor_index),...
         CIOD_N8_M4_num_error_Eve_parfor(parfor_index)] =...
        CIOD_IM_BER(num_iterations, N, M, P_tot_des, alpha, sigma2, SNRdB, mod_type);
    end
    
    CIOD_N8_M4(3, SNRdB_index) = mean(CIOD_N8_M4_BER_Bob_parfor);
    CIOD_N8_M4(4, SNRdB_index) = mean(CIOD_N8_M4_BER_Eve_parfor);
    CIOD_N8_M4(5, SNRdB_index) = sum(CIOD_N8_M4_num_error_Bob_parfor);
    CIOD_N8_M4(6, SNRdB_index) = sum(CIOD_N8_M4_num_error_Eve_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    fprintf("SNR = %d dB has just finished. \n\n", SNRdB);
end

fig4 = figure;
semilogy(SNRdB_array, CIOD_N8_M4(3, :), "d-", "Color", dark_red,...
                                              "MarkerEdgeColor", dark_red, ...
                                              "MarkerFaceColor", "r", ...
                                              "MarkerSize", 8);
hold on;
semilogy(SNRdB_array, conv_Na3_Nt2(3, :), "s-", "Color", dark_blue,...
                                                "MarkerEdgeColor", dark_blue, ...
                                                "MarkerFaceColor", "b", ...
                                                "MarkerSize", 8);
hold on;
semilogy(SNRdB_array, conv_Na6_Nt4(3, :), "^-", "Color", dark_green,...
                                                "MarkerEdgeColor", dark_green, ...
                                                "MarkerFaceColor", "g", ...
                                                "MarkerSize", 8);
hold on;
semilogy(SNRdB_array, eff_ala_M4(3, :), "o-", "Color", dark_yellow,...
                                              "MarkerEdgeColor", dark_yellow, ...
                                              "MarkerFaceColor", "y", ...
                                              "MarkerSize", 8);
hold on;
semilogy(SNRdB_array, CIOD_N8_M4(4, :), "d--", "Color", dark_red);
hold on;
semilogy(SNRdB_array, conv_Na3_Nt2(4, :), "s--", "Color", dark_blue);
hold on;
semilogy(SNRdB_array, conv_Na6_Nt4(4, :), "^--", "Color", dark_green);
hold on;
semilogy(SNRdB_array, eff_ala_M4(4, :), "o--", "Color", dark_yellow);

xlabel("$E_s / N_0$", "Interpreter", "latex");
ylabel("BER", "Interpreter", "latex");
legend("CIOD-IM, ($N=8$, QPSK, $l=2.5$)",...
       "Conventional, ($N_a=3$, $N_t=2$, BPSK, $l=2$)",...
       "Conventional, ($N_a=6$, $N_t=4$, BPSK, $l=3$)",...
       "Efficient Alamouti, ($N_a=3$, $N_t=2$, QPSK, $l=2.5$)",...
       "Location", "southwest", "FontSize", 11, "Interpreter", "latex");
ylim([1e-5 1]);
grid;

annotation("rectangle", [0.71, 0.81, 0.08, 0.08]);
annotation("textbox", [0.7, 0.71, 0.1, 0.1], "String", "Eve",...
                                              "FontSize", 10,...
                                              "Interpreter", "latex",...
                                              "HorizontalAlignment", "center",...
                                              "LineStyle", "none");
annotation("rectangle", [0.26, 0.70, 0.27, 0.05]);
annotation("textbox", [0.51, 0.65, 0.1, 0.1], "String", "Bob",...
                                              "FontSize", 10,...
                                              "Interpreter", "latex",...
                                              "HorizontalAlignment", "center",...
                                              "LineStyle", "none");

% Save As PDF
set(fig4, "Units", "Inches");
pos = get(fig4, "Position");
set(fig4, "PaperPositionMode", "Auto", "PaperUnits", "Inches", "PaperSize", [pos(3), pos(4)])
print(fig4, "Figure4", "-dpdf", "-r0")

%% FIGURE - 5
% Figure - 5a//////////////////////////////////////////////////////////////
SNRdB_array = (0 : 5 : 30);
len_SNRdB_array = length(SNRdB_array);

N = 4;
M = 4;
P_tot_des = 1;
alpha = 0.5;
mod_type = "PSK";

CIOD_0 = zeros(6, len_SNRdB_array);
CIOD_0(1, :) = SNRdB_array;
CIOD_0(2, :) = [7e2 2e3 2e4 5e5 5e6 5e7 7e2] / num_workers;

CIOD_003 = zeros(6, len_SNRdB_array);
CIOD_003(1, :) = SNRdB_array;
CIOD_003(2, :) = [5e2 3e3 3e3 6e3 8e3 8e3 9e3] / num_workers;

CIOD_01 = zeros(6, len_SNRdB_array);
CIOD_01(1, :) = SNRdB_array;
CIOD_01(2, :) = [5e2 7e2 9e2 2e3 3e3 3e3 3e3] / num_workers;

for SNRdB_index = 1 : len_SNRdB_array
    SNRdB = SNRdB_array(SNRdB_index);
    fprintf("SNR = %d dB has just started. \n", SNRdB);
    
    % CIOD - sigma2 = 0////////////////////////////////////////////////////
    fprintf("CIOD - sigma2 = 0 \n");
    num_iterations = CIOD_0(2, SNRdB_index);
    sigma2 = 0;
    parfor parfor_index = 1 : num_workers
        [CIOD_0_BER_Bob_parfor(parfor_index),...
         CIOD_0_BER_Eve_parfor(parfor_index),...
         CIOD_0_num_error_Bob_parfor(parfor_index),...
         CIOD_0_num_error_Eve_parfor(parfor_index)] =...
        CIOD_IM_BER(num_iterations, N, M, P_tot_des, alpha, sigma2, SNRdB, mod_type);
    end
    
    CIOD_0(3, SNRdB_index) = mean(CIOD_0_BER_Bob_parfor);
    CIOD_0(4, SNRdB_index) = mean(CIOD_0_BER_Eve_parfor);
    CIOD_0(5, SNRdB_index) = sum(CIOD_0_num_error_Bob_parfor);
    CIOD_0(6, SNRdB_index) = sum(CIOD_0_num_error_Eve_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    % CIOD - sigma2 = 0.03/////////////////////////////////////////////////
    fprintf("CIOD - sigma2 = 0.03 \n");
    num_iterations = CIOD_003(2, SNRdB_index);
    sigma2 = 0.03;
    parfor parfor_index = 1 : num_workers
        [CIOD_003_BER_Bob_parfor(parfor_index),...
         CIOD_003_BER_Eve_parfor(parfor_index),...
         CIOD_003_num_error_Bob_parfor(parfor_index),...
         CIOD_003_num_error_Eve_parfor(parfor_index)] =...
        CIOD_IM_BER(num_iterations, N, M, P_tot_des, alpha, sigma2, SNRdB, mod_type);
    end
    
    CIOD_003(3, SNRdB_index) = mean(CIOD_003_BER_Bob_parfor);
    CIOD_003(4, SNRdB_index) = mean(CIOD_003_BER_Eve_parfor);
    CIOD_003(5, SNRdB_index) = sum(CIOD_003_num_error_Bob_parfor);
    CIOD_003(6, SNRdB_index) = sum(CIOD_003_num_error_Eve_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    % CIOD - sigma2 = 0.1//////////////////////////////////////////////////
    fprintf("CIOD - sigma2 = 0.1 \n");
    num_iterations = CIOD_01(2, SNRdB_index);
    sigma2 = 0.1;
    parfor parfor_index = 1 : num_workers
        [CIOD_01_BER_Bob_parfor(parfor_index),...
         CIOD_01_BER_Eve_parfor(parfor_index),...
         CIOD_01_num_error_Bob_parfor(parfor_index),...
         CIOD_01_num_error_Eve_parfor(parfor_index)] =...
        CIOD_IM_BER(num_iterations, N, M, P_tot_des, alpha, sigma2, SNRdB, mod_type);
    end
    
    CIOD_01(3, SNRdB_index) = mean(CIOD_01_BER_Bob_parfor);
    CIOD_01(4, SNRdB_index) = mean(CIOD_01_BER_Eve_parfor);
    CIOD_01(5, SNRdB_index) = sum(CIOD_01_num_error_Bob_parfor);
    CIOD_01(6, SNRdB_index) = sum(CIOD_01_num_error_Eve_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    fprintf("SNR = %d dB has just finished. \n\n", SNRdB);
end

fig5 = figure;
subplot(1, 2, 1);
semilogy(SNRdB_array, CIOD_0(3, :), "d-", "Color", dark_red,...
                                          "MarkerEdgeColor", dark_red, ...
                                          "MarkerFaceColor", "r", ...
                                          "MarkerSize", 8);
hold on;
semilogy(SNRdB_array, CIOD_003(3, :), "o-", "Color", dark_yellow,...
                                            "MarkerEdgeColor", dark_yellow, ...
                                            "MarkerFaceColor", "y", ...
                                            "MarkerSize", 8);
hold on;
semilogy(SNRdB_array, CIOD_01(3, :), "*-", "Color", dark_green,...
                                           "MarkerEdgeColor", dark_green, ...
                                           "MarkerFaceColor", "g", ...
                                           "MarkerSize", 8);
hold on;
semilogy(SNRdB_array, CIOD_0(4, :), "d--", "Color", dark_red);
hold on;
semilogy(SNRdB_array, CIOD_003(4, :), "o--", "Color", dark_yellow);
hold on;
semilogy(SNRdB_array, CIOD_01(4, :), "*--", "Color", dark_green);

xlabel(["$E_s / N_0$";"(a)"], "Interpreter", "latex");
ylabel("BER", "Interpreter", "latex");
legend("$\sigma^2=0$",...
       "$\sigma^2=0.03$",...
       "$\sigma^2=0.1$",...
       "Location", "southwest", "FontSize", 11, "Interpreter", "latex");
ylim([1e-5 1]);   
grid;

annotation("rectangle", [0.313, 0.81, 0.08, 0.05]);
annotation("textbox", [0.3, 0.81, 0.1, 0.1], "String", "Eve",...
                                              "FontSize", 10,...
                                              "Interpreter", "latex",...
                                              "HorizontalAlignment", "center",...
                                              "LineStyle", "none");
annotation("rectangle", [0.167, 0.7, 0.04, 0.12]);
annotation("textbox", [0.14, 0.61, 0.1, 0.1], "String", "Bob",...
                                             "FontSize", 10,...
                                             "Interpreter", "latex",...
                                             "HorizontalAlignment", "center",...
                                             "LineStyle", "none");
% /////////////////////////////////////////////////////////////////////////

% Figure - 5b//////////////////////////////////////////////////////////////
SNRdB_array = (0 : 3 : 30);
len_SNRdB_array = length(SNRdB_array);

num_iterations = 2e3 / num_workers;
N = 4;
M = 4;
P_tot_des = 1;
sigma2 = 0;
mod_type = "PSK";

CIOD_01_Sec = zeros(4, len_SNRdB_array);
CIOD_01_Sec(1, :) = SNRdB_array;

CIOD_05_Sec = zeros(4, len_SNRdB_array);
CIOD_05_Sec(1, :) = SNRdB_array;

CIOD_09_Sec = zeros(4, len_SNRdB_array);
CIOD_09_Sec(1, :) = SNRdB_array;

for SNRdB_index = 1 : len_SNRdB_array
    SNRdB = SNRdB_array(SNRdB_index);
    fprintf("SNR = %d dB has just started. \n", SNRdB);
    
    % Alpha = 0.1//////////////////////////////////////////////////////////
    fprintf("Alpha = 0.1 \n");
    alpha = 0.1;
    parfor parfor_index = 1 : num_workers
        [CIOD_01_Rb_parfor(parfor_index),...
         CIOD_01_Re_parfor(parfor_index),...
         CIOD_01_Rs_parfor(parfor_index)] =...
        CIOD_IM_Secrecy(num_iterations, N, M, P_tot_des, alpha, sigma2, SNRdB, mod_type);
    end
    CIOD_01_Sec(2, SNRdB_index) = mean(CIOD_01_Rb_parfor);
    CIOD_01_Sec(3, SNRdB_index) = mean(CIOD_01_Re_parfor);
    CIOD_01_Sec(4, SNRdB_index) = mean(CIOD_01_Rs_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    % Alpha = 0.5//////////////////////////////////////////////////////////
    fprintf("Alpha = 0.5 \n");
    alpha = 0.5;
    parfor parfor_index = 1 : num_workers
        [CIOD_05_Rb_parfor(parfor_index),...
         CIOD_05_Re_parfor(parfor_index),...
         CIOD_05_Rs_parfor(parfor_index)] =...
        CIOD_IM_Secrecy(num_iterations, N, M, P_tot_des, alpha, sigma2, SNRdB, mod_type);
    end
    CIOD_05_Sec(2, SNRdB_index) = mean(CIOD_05_Rb_parfor);
    CIOD_05_Sec(3, SNRdB_index) = mean(CIOD_05_Re_parfor);
    CIOD_05_Sec(4, SNRdB_index) = mean(CIOD_05_Rs_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    % Alpha = 0.9//////////////////////////////////////////////////////////
    fprintf("Alpha = 0.9 \n");
    alpha = 0.9;
    parfor parfor_index = 1 : num_workers
        [CIOD_09_Rb_parfor(parfor_index),...
         CIOD_09_Re_parfor(parfor_index),...
         CIOD_09_Rs_parfor(parfor_index)] =...
        CIOD_IM_Secrecy(num_iterations, N, M, P_tot_des, alpha, sigma2, SNRdB, mod_type);
    end
    CIOD_09_Sec(2, SNRdB_index) = mean(CIOD_09_Rb_parfor);
    CIOD_09_Sec(3, SNRdB_index) = mean(CIOD_09_Re_parfor);
    CIOD_09_Sec(4, SNRdB_index) = mean(CIOD_09_Rs_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    fprintf("SNR = %d dB has just finished. \n\n", SNRdB);
end

subplot(1, 2, 2);
plot(SNRdB_array, CIOD_01_Sec(4, :), "d-", "Color", dark_red,...
                                       "MarkerEdgeColor", dark_red, ...
                                       "MarkerFaceColor", "r", ...
                                       "MarkerSize", 8);
hold on;
plot(SNRdB_array, CIOD_05_Sec(4, :), "o-", "Color", dark_yellow,...
                                       "MarkerEdgeColor", dark_yellow, ...
                                       "MarkerFaceColor", "y", ...
                                       "MarkerSize", 8);
hold on;
plot(SNRdB_array, CIOD_09_Sec(4, :), "*-", "Color", dark_green,...
                                       "MarkerEdgeColor", dark_green, ...
                                       "MarkerFaceColor", "g", ...
                                       "MarkerSize", 8);
                                        
xlabel(["$E_s / N_0$";"(b)"], "Interpreter", "latex");
ylabel("Ergodic Secrecy Rate", "Interpreter", "latex");
legend("$\alpha=0.1$",...
       "$\alpha=0.5$",...
       "$\alpha=0.9$",...
       "Location", "northwest", "FontSize", 9, "Interpreter", "latex");
ylim([0 2.5]);   
grid;
% /////////////////////////////////////////////////////////////////////////

% Save As PDF
set(fig5, "Units", "Inches");
pos = get(fig5, "Position");
set(fig5, "PaperPositionMode", "Auto", "PaperUnits", "Inches", "PaperSize", [pos(3), pos(4)])
print(fig5, "Figure5", "-dpdf", "-r0")

%% FIGURE - 6
SNRdB_array = (0 : 3 : 30);
len_SNRdB_array = length(SNRdB_array);

CIOD_N8_M4_Sec = zeros(4, len_SNRdB_array);
CIOD_N8_M4_Sec(1, :) = SNRdB_array;

eff_ala_M4_Sec = zeros(4, len_SNRdB_array);
eff_ala_M4_Sec(1, :) = SNRdB_array;

conv_Na3_Nt2_Sec = zeros(4, len_SNRdB_array);
conv_Na3_Nt2_Sec(1, :) = SNRdB_array;

for SNRdB_index = 1 : len_SNRdB_array
    SNRdB = SNRdB_array(SNRdB_index);
    fprintf("SNR = %d dB has just started. \n", SNRdB);
    
    % Conventional - Na=3, Nt=2////////////////////////////////////////////
    fprintf("Conventional - Na=3, Nt=2 \n");
    num_iterations = 1e7 / num_workers;
    Na = 3;
    Nt = 2;
    M = 2;
    P_tot_des = 1;
    alpha = 0.5;
    TAS_type = "SLNR";
    parfor parfor_index = 1 : num_workers
        [conv_Na3_Nt2_Rb_parfor(parfor_index),...
         conv_Na3_Nt2_Re_parfor(parfor_index),...
         conv_Na3_Nt2_Rs_parfor(parfor_index)] =...
        ConventionalSecrecy(num_iterations, Na, Nt, M, P_tot_des, alpha, SNRdB, TAS_type);
    end
    conv_Na3_Nt2_Sec(2, SNRdB_index) = mean(conv_Na3_Nt2_Rb_parfor);
    conv_Na3_Nt2_Sec(3, SNRdB_index) = mean(conv_Na3_Nt2_Re_parfor);
    conv_Na3_Nt2_Sec(4, SNRdB_index) = mean(conv_Na3_Nt2_Rs_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    % Efficient Alamouti///////////////////////////////////////////////////
    fprintf("Efficient Alamouti \n");
    num_iterations = 1e5 / num_workers;
    M = 4;
    P_tot_des = 1;
    alpha = 0.5;
    parfor parfor_index = 1 : num_workers
        [eff_ala_M4_Rb_parfor(parfor_index),...
         eff_ala_M4_Re_parfor(parfor_index),...
         eff_ala_M4_Rs_parfor(parfor_index)] =...
        EfficientAlamoutiSecrecy(num_iterations, M, P_tot_des, alpha, SNRdB);
    end
    eff_ala_M4_Sec(2, SNRdB_index) = mean(eff_ala_M4_Rb_parfor);
    eff_ala_M4_Sec(3, SNRdB_index) = mean(eff_ala_M4_Re_parfor);
    eff_ala_M4_Sec(4, SNRdB_index) = mean(eff_ala_M4_Rs_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    % CIOD - N=8///////////////////////////////////////////////////////////
    fprintf("CIOD - N=8 \n");
    num_iterations = 2e3 / num_workers;
    N = 8;
    M = 4;
    P_tot_des = 1;
    alpha = 0.5;
    sigma2 = 0;
    mod_type = "PSK";
    parfor parfor_index = 1 : num_workers
        [CIOD_N8_M4_Rb_parfor(parfor_index),...
         CIOD_N8_M4_Re_parfor(parfor_index),...
         CIOD_N8_M4_Rs_parfor(parfor_index)] =...
        CIOD_IM_Secrecy(num_iterations, N, M, P_tot_des, alpha, sigma2, SNRdB, mod_type);
    end
    CIOD_N8_M4_Sec(2, SNRdB_index) = mean(CIOD_N8_M4_Rb_parfor);
    CIOD_N8_M4_Sec(3, SNRdB_index) = mean(CIOD_N8_M4_Re_parfor);
    CIOD_N8_M4_Sec(4, SNRdB_index) = mean(CIOD_N8_M4_Rs_parfor);
    % /////////////////////////////////////////////////////////////////////
    
    fprintf("SNR = %d dB has just finished. \n\n", SNRdB);
end

fig6 = figure;
plot(SNRdB_array, CIOD_N8_M4_Sec(4, :), "d-", "Color", dark_red,...
                                              "MarkerEdgeColor", dark_red, ...
                                              "MarkerFaceColor", "r", ...
                                              "MarkerSize", 8);
hold on;
plot(SNRdB_array, conv_Na3_Nt2_Sec(4, :), "s-", "Color", dark_blue,...
                                                "MarkerEdgeColor", dark_blue, ...
                                                "MarkerFaceColor", "b", ...
                                                "MarkerSize", 8);
hold on;
plot(SNRdB_array, eff_ala_M4_Sec(4, :), "o-", "Color", dark_yellow,...
                                              "MarkerEdgeColor", dark_yellow, ...
                                              "MarkerFaceColor", "y", ...
                                              "MarkerSize", 8);

xlabel("$E_s / N_0$", "Interpreter", "latex");
ylabel("Ergodic Secrecy Rate", "Interpreter", "latex");
legend("CIOD-IM, ($N=8$, QPSK, $l=2.5$)",...
       "Conventional, ($N_a=3$, $N_t=2$, BPSK, $l=2$)",...
       "Efficient Alamouti, ($N_a=3$, $N_t=2$, QPSK, $l=2.5$)",...
       "Location", "southeast", "FontSize", 11, "Interpreter", "latex");
grid;

% Save As PDF
set(fig6, "Units", "Inches");
pos = get(fig6, "Position");
set(fig6, "PaperPositionMode", "Auto", "PaperUnits", "Inches", "PaperSize", [pos(3), pos(4)])
print(fig6, "Figure6", "-dpdf", "-r0")