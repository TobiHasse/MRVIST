function reach_length = channel_reach_length( riv, B, us_lim_meters, ...
    pix_per_chan, x_start, x_end )
% Purpose:  Create frequency distributions for the length of channel within
%           each reach and save figure 25 from chapter 4
% Inputs    riv             struct with river planforms
%           B               channel half width (meters)
%           usl_lim_meters  The extreme upstream limit of the model domain
%                           (in meters)
%           pix_per_chan    channel width in pixels
%           x_start x_end   start and end pixels of the reaches
%           ideally the pixel coordinate of the border between reaches
%           would be input as well.
%
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     2018-2019, edited October 2021

% public function dependencies
%           distributionPlot from the MATLAB file Exchange FEX
%           http://www.mathworks.com/matlabcentral/fileexchange/
%           37105-plot-spread-points-beeswarm-plot

%% find the channel length in the upstream and downstream reaches

% HARD CODE boundaries
% us = 1542.63;      % upstream in meters
% ms = 15276.63;      % 'mid'point separating reaches
% ds = 28086.63;      % dowstream end of downstream reach
ms_pixel = 1325;    % 'mid'point separating reaches in pixels

pixel_sz = 2 * B / pix_per_chan;

% determine meter coordinate for reach boundaries
% -1 because pixel index is inclusive
us = us_lim_meters + pixel_sz * ( x_start -1 ); % upstream start
ms = us_lim_meters + pixel_sz *   ms_pixel;     % between reaches
ds = us_lim_meters + pixel_sz *   x_end;        % downstream end

reach_length = zeros(size(2,numel(riv)));

for t = 1:6845  % numel(riv)
    X = riv(t).Xcl;
    Y = riv(t).Ycl;
    % channel length at each node ( algorithm from curvars_TRH )
    dS = sqrt( diff(X).^2 + diff(Y).^2);
    dS = [dS;0];                         % dS is one element longer than X
    idus = find( X>us & X<ms);           % find indixes within reach
    idds = find( X>ms & X<ds);
    reach_length(1,t) = sum( dS(idus));
    reach_length(2,t) = sum( dS(idds));
end


% DEBUG make some plots
hf = figure(20);            % planform with channel nodes marked by reach
    plot(X,Y)
    axis equal
    hold on
    plot(X(idus),Y(idus), 'marker','x')
    plot(X(idds),Y(idds), 'marker','o')
    hold off
    hf.Units = 'inches';
    hf.Position = [.2 .2 6.25 4.2];
    set( hf, 'color', 'w', 'PaperPosition', hf.Position );
    drawnow
    
figure(21)                  % time series
    plot(reach_length')
    legend('Upstream','Dowstream')
    ylabel('Reach length (m)')

figure(22)
    boxplot(reach_length','labels',{'Upstream','Downstream'})
    ylabel('Reach length (m)')

figure(23) % in case the formatted figure fails
    distributionPlot(reach_length')
%% Make pretty plot and save
hf = figure(24);
    clf
    hdv = distributionPlot(reach_length(:,4000:6845)',...
        'color',[.3 .3 .3],'showMM',6,'xnames',{'Upstream','Downstream'});
    box on
    set(hdv{2}(1), 'color', [1 1 1 ]);  % change the color of the x
    ylabel('Length of channel in reach (km)')

    hold on % add markers for mean
    hmv = plot(hdv{3}.XTick,mean(reach_length(:,4000:6845),2),...
        'xw','MarkerSize',12,'linewidth',2);
    hold off
    set(hdv{2}(:), 'color', [0 0 0 ]); % set quantiles to black                  

    ha = hf.Children;
    ha.YTickLabel = ha.YTick/1000;
    ha.XLim = [ 0.3 2.7 ];
    ha.Position = [.1 .09 .865 .875];
    hf.Units = 'inches';
    hf.Position = [4 4  5.8 6];
    set( hf, 'color', 'w', 'PaperPosition', hf.Position );
drawnow
outfile = 'reach length violins 4000 6845.png';
print('-painters', '-dpng', '-r600', outfile)
% plot2svg([files{i},' reach length violins 4000 6845.svg'])
end % function channel_reach_length
