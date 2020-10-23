% Script to evaluate the LOS probability for 3GPP model.
% Looks quite similar to what is presented in paper, at least for UMi
% UMa a bit different, might be because of a dependence on Receiver height
% as well as 2D outdoor distance

close all

rx = Receiver(0,0,0);
tx = Transmitter(0,1,0);
model_umi = Model_3gpp();
model_uma = Model_3gpp("UMa");

N = 300;
range = 1:N;
psi = zeros(N,1);
psa = zeros(N,1);

for i = range
    [~, pLOS] = model_umi.rollLOS(rx,tx);
    psi(i)=pLOS;
    [~, pLOS] = model_uma.rollLOS(rx,tx);
    psa(i)=pLOS;
    tx=tx.movePosition(0,1,0);
end

plot(range,psi);
hold;
plot(range,psa);
grid on;
legend("UMI", "UMA");
title("LOS Probability in 3GPP Model");
xlabel("2D-Distance between Tx and Rx");