% Script to evaluate the LOS probability for 3GPP model.
% Looks quite similar to what is presented in paper, at least for UMi
% UMa a bit different, might be because of a dependence on Receiver height
% as well as 2D outdoor distance

close all

rx = Receiver(0,0,1.5);
tx = Transmitter(0,1,4);
model_umi = Model_3gpp("UMi");
model_uma = Model_3gpp("UMa");

steps = 5000;
ys = 1:steps;
pls_UMi_LOS = zeros(steps,1);
pls_UMa_LOS = zeros(steps,1);
isLOS = true;
for y = ys   
    tx.setPosition(0,y+9,10)
    pls_UMi_LOS(y) = model_umi.largeScalePathloss(rx,tx,isLOS);
    pls_UMa_LOS(y) = model_uma.largeScalePathloss(rx,tx,isLOS);
    pls_UMi_NLOS(y) = model_umi.largeScalePathloss(rx,tx,false);
    pls_UMa_NLOS(y) = model_uma.largeScalePathloss(rx,tx,false);
end
semilogx(ys+9,pls_UMi_LOS,"b");
hold on
semilogx(ys+9,pls_UMa_LOS,"r");
hold on
semilogx(ys+9,pls_UMi_NLOS,"--b");
hold on
semilogx(ys+9,pls_UMa_NLOS,"--r");
hold on

grid on
legend("UMi LOS", "UMa LOS", "UMi NLOS", "UMa NLOS", 'Location','northwest');
title("Pathloss in 3GPP Model (hUE = 1.5m, hBS = 4m)");
xlabel("2D-Distance between Tx and Rx");
