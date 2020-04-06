clc
clear
close

E = 1;
k = 4;
r = 8;
N = 500000; % количество сообщений при моделировании

%gX_1 = x^8 + x^4 + x^3 + x^2 + 1
gX_1 = [1 0 0 0 1 1 1 0 1];
codes_1 = codeBook(k, gX_1);
codes_BPSK_1 = codes_1.*-2 + 1;
%gX_2 = x^8 + x^7 + x^6 + x^4 + x^2 + 1
gX_2 = [1 1 1 0 1 0 1 0 1];
codes_2 = codeBook(k, gX_2);
codes_BPSK_2 = codes_2.*-2 + 1;
%gX_3 = x^8 + x^5 + x^4 + 1
gX_3 = [1 0 0 1 1 0 0 0 1];
codes_3 = codeBook(k, gX_3);
codes_BPSK_3 = codes_3.*-2 + 1;

SNRdB = -20 : 2 : 0;
PeBits = zeros (3, length(SNRdB));
PEDs = zeros (3, length(SNRdB));

for i = 1 : length(SNRdB)
    disp(SNRdB(i));
    SNR = 10.^(SNRdB(i)/10);
    sigma = sqrt(E / (2*SNR));
    
    [PeBit_1, PED_1] = model(k, gX_1, codes_1, codes_BPSK_1, sigma, N);
    [PeBit_2, PED_2] = model(k, gX_2, codes_2, codes_BPSK_2, sigma, N);
    [PeBit_3, PED_3] = model(k, gX_3, codes_3, codes_BPSK_3, sigma, N);
    
    PeBits(1, i) = PeBit_1;
    PeBits(2, i) = PeBit_2;
    PeBits(3, i) = PeBit_3;
    
    PEDs(1, i) = PED_1;
    PEDs(2, i) = PED_2;
    PEDs(3, i) = PED_3;
end

% ---- Вероятности ошибки на бит для BPSK.  ----
SNRtheor = 10.^(SNRdB/10);
PeBitstheor = qfunc(sqrt(2*SNRtheor));

% ---- Вероятность ошибки декодирования CRC. ----
%Точная вероятность ошибки декодирования CRC:
n = r + k;
PEDsExact = zeros (3, length(SNRdB));
A_1 = A_func(codes_1); 
d_min_1 = min(sum(codes_1(2:end, :),2));

for i = 1 : length(SNRdB)
    for j = d_min_1 : n
        PEDsExact(1, i) = PEDsExact(1, i) + A_1(j + 1) * ...
            PeBitstheor(i)^j * (1 - PeBitstheor(i))^(n - j);
    end
end

A_2 = A_func(codes_2);
d_min_2 = min(sum(codes_2(2:end, :),2));

for i = 1 : length(SNRdB)
    for j = d_min_2 : n
        PEDsExact(2, i) = PEDsExact(2, i) + A_2(j + 1) * ...
            PeBitstheor(i)^j * (1 - PeBitstheor(i))^(n - j);
    end
end

A_3 = A_func(codes_3);
d_min_3 = min(sum(codes_3(2:end, :),2));

for i = 1 : length(SNRdB)
    for j = d_min_3 : n
        PEDsExact(3, i) = PEDsExact(3, i) + A_3(j + 1) * ...
            PeBitstheor(i)^j * (1 - PeBitstheor(i))^(n - j);
    end
end

figure();
axis('square');
semilogy(SNRdB, PeBitstheor, 'r-', SNRdB, PeBits(1, :),'b--', ...
    SNRdB, PeBits(2, :), 'y.-', SNRdB, PeBits(3, :), 'c*-');
xlabel('SNRdB'); 
ylabel('PeBit');
legend ({'Теоретическое значение вероятности ошибки на бит', ...
    'Практика для x^8 + x^4 + x^3 + x^2 + 1', ...
    'Практика для x^8 + x^7 + x^6 + x^4 + x^2 + 1', ...
    'Практика для x^8 + x^5 + x^4 + 1'}, ...
'Location','southwest')

figure();
axis('square');
semilogy(SNRdB, PEDs(1, :), 'b+-', SNRdB, PEDsExact(1, :), 'bx-', ...
    SNRdB, PEDs(2, :), 'y+-', SNRdB, PEDsExact(2, :), 'yx-', ...
    SNRdB, PEDs(3, :), 'c+-', SNRdB, PEDsExact(3, :), 'cx-');
xlabel('SNRdB'); 
ylabel('PED');

legend({'Практика для x^8 + x^4 + x^3 + x^2 + 1', ...
    'Теория для x^8 + x^4 + x^3 + x^2 + 1', ...
    'Практика для x^8 + x^7 + x^6 + x^4 + x^2 + 1', ...
    'Теория для x^8 + x^7 + x^6 + x^4 + x^2 + 1'...
    'Практика для x^8 + x^5 + x^4 + 1', ...
    'Теория для x^8 + x^5 + x^4 + 1'},...
'Location','east')