fc = 29;
muDS = -7.31;
sigDS = 0.6^2;

rDS = 2.06;

it = 10000;
N_clusters = 6;
M_rays = 10;

DS_log = normrnd(muDS, sigDS, it,1);
DS = 10.^DS_log;

delays_clusters = zeros(it,N_clusters);
delays_rays = zeros(it,N_clusters,M_rays);
powers_clusters = zeros(it,N_clusters);
powers_rays = zeros(it,N_clusters,M_rays);

rmss =  zeros(it,1);
xis = [.3 6 9 12 15];

for i = 1:it
    tau = -rDS * DS(i) * log(rand(N_clusters,1));
    delays_clusters(i,:) = sort(tau-min(tau))*1e9;
    
    delays_rays(i,:,:) = delays_clusters(i,:)'*ones(1,M_rays);
    temp = delays_rays(i,:,1) + [0 0 0 0 5 5 10 10 5 0]';
    delays_rays(i,:,:) = reshape(temp',1,N_clusters,M_rays);
    
    powers_clusters(i,:) = exp(-delays_clusters(i,:)/1e9 * (rDS -1) / rDS / DS(i)) .* 10.^(normrnd(0,xi,N_clusters,1)/10)';
    powers_clusters(i,:) = powers_clusters(i,:) / sum(powers_clusters(i,:));
    powers_rays(i,:,:) = powers_clusters(i,:)'*ones(1,M_rays)/M_rays;
    
    rmss(i) = RMS_delaySpread(reshape(powers_rays(i,:,:),N_clusters*M_rays,1), reshape(delays_rays(i,:,:),N_clusters*M_rays,1));
end

cdfplot(DS*1e9)
xlim([0 180])
