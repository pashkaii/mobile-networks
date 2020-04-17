clc; clear; close all

g_x = [1 0 1 0 0 1 1 0 1 0 1 1 1 1 0 0 1];
k = 8;                  % Длина сообщения
r = length(g_x) - 1;    % Длина CRC части
num_messages = 2^k;     % Количество возможных сообщений
G = [
[1, 1, 0, 1]
[1, 0, 1, 1]
[1, 0, 0, 0]
[0, 1, 1, 1]
[0, 1, 0, 0]
[0, 0, 1, 0]
[0, 0, 0, 1]
];
H = [
[1, 0, 1, 0, 1, 0, 1]
[0, 1, 1, 0, 0, 1, 1]
[0, 0, 0, 1, 1, 1, 1]
];

messages = de2bi((0:num_messages-1)');     % Все возможные сообщения
codewords = zeros(num_messages, k+r);      % Все кодовые слова с CRC
mes_6_4 = zeros(num_messages, 6, 4);     
channel_6_7 = zeros(num_messages, 6, 7); % Замодулированные BPSK (-1, 1)
                                           % без добавленных 3х нулей,
% Генерация                                  которые отправляются в канал
for j=1:num_messages
    [~, c] = gfdeconv([zeros(1, r), messages(j,:)], g_x);   % CRC часть
    codewords(j,:) = xor([zeros(1, r), messages(j,:)], ...  % Сообщение +
                         [c, zeros(1, k+r - length(c))]);   % + CRC
    mes_6_4(j, :, :) = reshape(codewords(j, :), 4, 6)';     % + нули
    for i=1:6
        channel_6_7(j, i, :) = mod(reshape(mes_6_4(j, i, :), 1, 4) ...
                                    * G', 2)*(-2) + 1; 
    end
end

% Теоретический расчет
SNRdB = -10:10;
SNR = 10.^(SNRdB./10);
Pe_bit_theor = qfunc(sqrt(2.*SNR));         
Ped_theor_up = ones(1, length(SNR)).*1/2^r; 
Ped_theor = zeros(1,length(SNRdB));
w = sum(codewords(2:end, :)');
d = min(w);
for i=1:length(SNR)
    Pe_theoretical = 0;
    for j=d:(k + r)
        Pe_theoretical = Pe_theoretical + sum(w == j) * ...
                (Pe_bit_theor(i)^j) * ((1 - Pe_bit_theor(i))^((k + r)-j));
    end
    Ped_theor(i) = Pe_theoretical;
end

% Моделирование
Pe_bit = zeros(1,length(SNRdB));  % Вероятность ошибки на бит
Pe_bit_after_Hamming = zeros(1,length(SNRdB));
Ped = zeros(1,length(SNRdB));     % Вероятность ошибки декодирования
T = zeros(1,length(SNRdB));       % Пропускная способность канала
for i=1:length(SNR)
    disp(i);
    sigma = sqrt(1/(2*SNR(i)));
    nTests = 0; nSent = 0;
    nErrDecode = 0;
    nErrBits = 0; nErrBits_H = 0;
    messages_to_send = 0;
    while messages_to_send < 10*(SNR(i)+11)
        messages_to_send = messages_to_send + 1;
        rand_ind = randi(num_messages, 1);
        c = reshape(channel_6_7(rand_ind, :, :), 6, 7);
        while 1
            nTests = nTests + 1;
            AWGN = c + sigma*randn(6, 7);
            [corrected, nErr, v_] = correctError(AWGN, H, c);
            nErrBits = nErrBits + nErr;
            nErrBits_H = nErrBits_H + v_;
            
            [~, S] = gfdeconv(corrected, g_x); 
            v = sum(xor(codewords(rand_ind, :), corrected));
            
            if (bi2de(S) == 0)
                if (v_ > 0)
                    nErrDecode = nErrDecode + 1;
                end
                nSent = nSent + 1;
                break;
            end
        end
    end
    Ped(i) = nErrDecode/nTests;
    Pe_bit(i) = nErrBits/(nTests*42);
    Pe_bit_after_Hamming(i) = nErrBits_H/(nTests*42);
    T(i) = k * nSent / (42 * nTests);
end

figure;
subplot(2, 2, 1);
semilogy(SNRdB, Ped, 'ko', ...
         SNRdB, Ped_theor_up, 'r.-', ...
         SNRdB, Ped_theor, 'b.-')
legend('Ped', 'Ped theor up', 'Ped theor');
subplot(2, 2, 3);
semilogy(SNRdB, Pe_bit, 'k', ...
         SNRdB, Pe_bit_after_Hamming, 'r.--', ...
         SNRdB, Pe_bit_theor, 'b.-');
legend('Pe bit', 'Pe bit after Hamming', 'Pe bit theor');
xlabel('E/N_0, dB')
subplot(1, 2, 2);
plot(SNRdB, T, 'b.-');
ylabel('T, пропускная способность');
xlabel('E/N_0, dB')

% Синдромное декодирование
function [concatenated, nErrBits, v] = correctError(AWGN, H, c)
	unBPSK = AWGN < 0;
    nErrBits = sum(sum(xor(c < 0, unBPSK)));
    corrected = zeros(6, 7);
    
    for part=1:6
        corr = unBPSK(part, :);     
        add = zeros(1, 7);
        S = bi2de(mod(corr*H',2));
        if S ~= 0
            add(1, S) = 1;
            corr = xor(corr, add);
        end
        corrected(part, :) = corr(1, :);
    end
    v = sum(sum(xor(corrected, c<0)));
    concatenated = [corrected(1, 3) corrected(1, 5:end) ...
                    corrected(2, 3) corrected(2, 5:end) ...
                    corrected(3, 3) corrected(3, 5:end) ...
                    corrected(4, 3) corrected(4, 5:end) ...
                    corrected(5, 3) corrected(5, 5:end) ...
                    corrected(6, 3) corrected(6, 5:end)];
end