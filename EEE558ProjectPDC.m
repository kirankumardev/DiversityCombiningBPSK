clear
% Number of Bits or Symbols
N = 10^6; 
% Transmitter part
% Generating 0,1 with equal probability
ip = rand(1,N)>0.5;
% BPSK modulation 0 -> -1; 1 -> 0
s = 2*ip-1; 
% Array of Diversity M = 1 and M = 2
nRx =  [1 4];
% Multiple SNR(dB) values
Eb_N0_dB = [0:50]; 
iter = 0;
while(iter<10^6)
for jj = 1:length(nRx)
        for ii = 1:length(Eb_N0_dB)
            % White gaussian noise, 0dB variance
            n = 1/sqrt(2)*[randn(nRx(jj),N) + j*randn(nRx(jj),N)]; 
            % Channel 
            h = 1/sqrt(2)*[randn(nRx(jj),N) + j*randn(nRx(jj),N)]; 
            % Channel and Noise addition
            sD = kron(ones(nRx(jj),1),s);
            y = h.*sD + 10^(-Eb_N0_dB(ii)/20)*n;
            % Post Detection combining 
            yHat =  sign(sum(conj(h).*y,1)./sum(h.*conj(h),1)); 
            % Receiver - hard decision decoding
            ipHat = real(yHat)>0;
            % Counting the Errors
            nErr(jj,ii) = size(find([ip- ipHat]),2); 
            iter = iter+1;
        end
    end
end
% Simulated BER
simBer = nErr/N; 
EbN0Lin = 10.^(Eb_N0_dB/10);
% Theoretical BER
theoryBer_nRx1 = 0.5.*(1-1*(1+1./EbN0Lin).^(-0.5)); 
p = 1/2 - 1/2*(1+1./EbN0Lin).^(-1/2);
theoryBer_nRx2 = p.^2.*(1+2*(1-p)); 
%Ploting the Graph
close all
figure
semilogy(Eb_N0_dB,theoryBer_nRx1,'bp-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer(1,:),'mo-','LineWidth',2);
semilogy(Eb_N0_dB,theoryBer_nRx2,'rd-','LineWidth',2);
semilogy(Eb_N0_dB,simBer(2,:),'ks-','LineWidth',2);
axis([0 50 10^-6 0.5])
grid on
legend('M=1 (theoretical)', 'M=1 (simulation)', 'M=4 (theoretical)', 'M=4 (simulation)');
xlabel('SNR(Eb/No) in dB');
ylabel('Bit Error Rate');
title('BER for BPSK with Post Detection Combining');