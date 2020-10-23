% Script to debug the smallScaleFading of 3gpp model.

close all

rx = Receiver(0,0,1.5);
tx = Transmitter(0,150,20);
model_umi = Model_3gpp("UMi");
model_uma = Model_3gpp("UMa");

tx.setPosition(0,50,10)

% channel = model_umi.applyModel(rx,tx);
it = 10000;
vals = zeros(it,1);

Ns_clusters = zeros(it,1);
DS = zeros(it,1);
ASD = zeros(it,1);

for i = 1:it
    results = model_umi.smallScaleParams(rx,tx,false);
    
    DS(i)= results('DS');
    ASD(i) = results('ASD');
    SF(i) = results('SF');
end
mean(SF)
% mean(ASD)
% mean(DS)
cdfplot(DS);