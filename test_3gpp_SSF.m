% Script to debug the smallScaleFading of 3gpp model.

close all

rx = Receiver(0,0,1.5);
tx = Transmitter(0,150,20);
model_umi = Model_3gpp("UMi");
model_uma = Model_3gpp("UMa");

tx.setPosition(0,50,10)

channel = model_umi.applyModel(rx,tx);
it = 1000;
vals = zeros(it,1);

Ns_clusters = zeros(it,1);
isLoss = zeros(it,1);
DS = zeros(it,1);
allDelays = cell(it);
allPowers = cell(it);

for i = 1:it
    results = model_umi.applyModel(rx,tx,false);
    
    parameters = results('parameters');
    N_clusters = results('N_clusters');
    clusterDelays = results('clusterDelays');
    clusterPowers = results('clusterPowers');
    isLOS = results('isLOS');
    Hust = results('channel');
    angles = results('angles');   
    
    DS(i) = parameters('DS');
    Ns_clusters(i) = N_clusters;
    isLoss(i) = isLOS;
    allDelays{i} = clusterDelays;
    allPowers{i} = clusterPowers;
end

rmsDelaySpreads = zeros(it,1);
delays = zeros(1,1);
for i=1:it
    rmsDelaySpreads(i) = RMS_delaySpread(allPowers{i}, allDelays{i});
    delays = cat(1,delays,allDelays{i});
end

cdfplot(rmsDelaySpreads*1e9)
title('RMS delay spread for UMi, 28 GHz, 3GPP Model, NLOS, d=150m')
xlabel('RMS delay spread [ns]')
ylabel('Cumulative Distribution Function (CDF)')
% corrcoef(vals)
