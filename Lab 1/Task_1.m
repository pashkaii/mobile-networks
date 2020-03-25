clc
clear
close

E = 1;
k = 4;
N = 1000000; % ���������� ��������� ��� �������������

%gX = x^16+x^13+x^12+x^11+x^10+x^8+x^6+x^5+x^2+1
gX = [1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1];

%gX = x^4+x^3+x^2+1
%gX = [1, 1, 1, 0, 1];
codes = codeBook(k, gX);

SNRdB = -20 : 0;
PeBits = zeros (1, length(SNRdB));
PEDs = zeros (1, length(SNRdB));

parfor i = 1 : length(SNRdB)
    disp(SNRdB(i));
    SNR = 10.^(SNRdB(i)/10);
    sigma = sqrt(E / (2*SNR));
    
    [PeBit, PED] = model(k, gX, codes, sigma, N);
    
    PeBits(1, i) = PeBit;
    PEDs(1, i) = PED;
end

% ---- ����������� ������ �� ��� ��� BPSK.  ----
SNRtheor = 10.^(SNRdB/10);
PeBitstheor = qfunc(sqrt(2*SNRtheor));


% ---- ����������� ������ ������������� CRC. ----
% ������� ������� (���������������):
r = length(gX) - 1;
PEDs_asimp = (ones(1, length(SNRdB)) ./ 2).^ r;

%������ ����������� ������ ������������� CRC:
A = A_func(codes);
d_min = min(sum(codes(2:end, :),2));
n = r + k;

PEDsExact = zeros (1, length(SNRdB));
for i = 1 : length(SNRdB)
    for j = d_min : n
        PEDsExact(1, i) = PEDsExact(1, i) + A(j + 1) * ...
            PeBitstheor(i)^j * (1 - PeBitstheor(i))^(n - j);
    end
end

%����� ������ ������� ������� ������ ������������� CRC:
PEDs_asimp_2 = (2^k - 1).*PeBitstheor.^d_min;

figure();
axis('square');
semilogy(SNRdB, PeBits, 'b.-', SNRdB, PeBitstheor, 'r-');
xlabel('SNRdB'); 
ylabel('PeBit');
legend ({'������������ �������� ����������� ������ �� ���', ...
    '������������� �������� ����������� ������ �� ���'}, ...
'Location','southwest')

figure();
axis('square');
hold on
semilogy(SNRdB, PEDs, 'b+-');
semilogy(SNRdB, PEDs_asimp, 'r.-');
semilogy(SNRdB, PEDsExact, 'cx-');
semilogy(SNRdB, PEDs_asimp_2, 'm--');
xlabel('SNRdB'); 
ylabel('PED');

legend({'������������ �������� ����������� ������ ��������', ...
    '��������������� ������� ������� ����������� ������ ��������', ...
    '������ �������� ����������� ������ ��������', ...
    '����� ������ ������� ������� ����������� ������ ��������'},...
'Location','east')