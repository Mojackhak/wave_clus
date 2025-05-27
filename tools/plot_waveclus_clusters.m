function plot_waveclus_clusters(spikes, cluster_class)
% plot_waveclus_clusters  – Quick clone of the Wave_Clus cluster panels
%
%   spikes        N × L   matrix of waveforms (µV or AD units)
%   cluster_class N × 2   [cluster_ID , spike_time]

%% ---- basic bookkeeping -----------------------------------------------
cluster_id = cluster_class(:,1);
good_clus   = setdiff( unique(cluster_id) , 0 );   % drop class 0 (noise)
nClus = numel(good_clus);

% arrange subplots in a square-ish grid
nCols = ceil(sqrt(nClus));
nRows = ceil(nClus / nCols);
colors = lines(max( nClus , 7 ));                  % default MATLAB colormap

figure('Color','w'); set(gcf,'Name','Wave_Clus-style plot');

%% ---- plot each cluster ------------------------------------------------
for ii = 1:nClus
    c = good_clus(ii);
    idx = find(cluster_id == c);
    wav = spikes(idx,:);         % all spikes in this cluster
    mu  = mean(wav,1);           % mean waveform
    
    % -------- waveform axes ---------
    subplot(nRows, nCols, ii); hold on
    plot(wav','Color',[ colors(ii,:) 0.15 ])      % raw spikes (faint)
    plot(mu,'k','LineWidth',1.5)                  % mean waveform
    title(sprintf('Cluster %d :  # %d',c,numel(idx)), ...
          'FontWeight','bold','FontSize',10)
    
    xlim([1 size(spikes,2)])
    ylim( [min(mu)-30 , max(mu)+30] )             % tweak if needed
    xlabel('sample'); ylabel('µV')
end

end 