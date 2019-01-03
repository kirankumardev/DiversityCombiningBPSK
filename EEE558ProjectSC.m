clear
% Number of Bits or Symbols
N = 10^6; 
% Transmitter part
% Generating 0,1 with equal probability
ip = rand(1,N)>0.5; 
% BPSK modulation 0 -> -1; 1 -> 0
s = 2*ip-1; 
% Array of Diversity factor M = 1 and M = 4
nRx = [1 4];
% Multiple SNR(dB) values
Eb_N0_dB =[0:50];
iter = 0;
while(iter<10^6)
    for jj = 1:length(nRx)
        for ii = 1:length(Eb_N0_dB)
            % White Gaussian Noise, 0dB variance
            n = 1/sqrt(2)*[randn(nRx(jj),N) + j*randn(nRx(jj),N)]; 
            % Channel
            h = 1/sqrt(2)*[randn(nRx(jj),N) + j*randn(nRx(jj),N)]; 
            % Channel and noise Noise addition
            sD = kron(ones(nRx(jj),1),s);
            y = h.*sD + 10^(-Eb_N0_dB(ii)/20)*n;
            % Finding the Power of the channel on all receiver
            hPower = h.*conj(h);
            % Finding the maximum power
            [hMaxVal ind] = max(hPower,[],1);
            hMaxValMat = kron(ones(nRx(jj),1),hMaxVal);
            % Selecting the chain with the maximum power
            ySel = y(hPower==hMaxValMat);
            hSel = h(hPower==hMaxValMat);
            % Selection Combining
            yHat = ySel./hSel;
            % Reshaping to get the proper matrix dimension 
            yHat = reshape(yHat,1,N); 
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
theoryBer_nRx2 = 0.5.*(1-2*(1+1./EbN0Lin).^(-0.5) + (1+2./EbN0Lin).^(-0.5));
% Plotting the graph
close all
figure
semilogy(Eb_N0_dB,theoryBer_nRx1,'bp-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer(1,:),'mo-','LineWidth',2);
semilogy(Eb_N0_dB,theoryBer_nRx2,'rd-','LineWidth',2);
semilogy(Eb_N0_dB,simBer(2,:),'ks-','LineWidth',2);
axis([0 50 10^-6 0.5])
grid on
legend('M=1 (theoretical)', 'M=1 (Simulation)', 'M=4 (theoretical)', 'M=4 (Simulation)');
xlabel('SNR(Eb/No) in dB');
ylabel('Bit Error Rate');
title('BER for BPSK with Selection Combining');