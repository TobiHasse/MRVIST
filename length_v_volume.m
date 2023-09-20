function length_v_volume( file, reach_length, dt,...
    pt_bar_dist, vert_dist, pb_age_dist, vd_age_dist )
% Purpose:  Create frequency distributions for the length of channel within
%           each reach and save figure 25 from chapter 4
% Inputs    files           cell array of input files
%           reach_length    length of reach(es) in meters
%           *_dist          storage time distributions
%           *_age_*         age distributions
%
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     August, 2019, edited October 2021

% public function dependencies
%           distributionPlot from the MATLAB file Exchange FEX
%           http://www.mathworks.com/matlabcentral/fileexchange/
%           37105-plot-spread-points-beeswarm-plot
%
% private functions dependencies
% alpha_marker          set alpha value on figure markers ~Tobias Hasse

% Compare reach length to deposition & eroding sediment

% remove terraces ( likely redundant ) the original storage distribution 
% includes 'terrace' from prior to model start. removed from updated files
pt_bar_dist = fliplr(tril(fliplr(pt_bar_dist),-2)); %-2
pb_age_dist = fliplr(tril(fliplr(pb_age_dist),-2)); %-2 

% the age distributions include the amount deposited 
depos_vol.v = vd_age_dist(:,6845);
depos_vol.p = pb_age_dist(:,6845);
erode_vol.v = sum(vert_dist,2);
erode_vol.p = sum(pt_bar_dist,2);

px = [1:length(depos_vol.v)] * dt;

% time series plots
figure(10)
%     subplot(2,1,1)
plot(px,[.1*reach_length', depos_vol.v, erode_vol.v])
%         hold on
%         plot(both_reach_ln, erode_vol.v)
%         hold off
xlabel('reach length and overbank deposition & erosion')
legend('length','deposition','erosion')
figure(11)
%     subplot(2,1,2)
plot(px,[.1*reach_length', depos_vol.p, erode_vol.p])

%         plot(both_reach_ln, depos_vol.p)
%         hold on
%         plot(both_reach_ln, erode_vol.p)
%         hold off
xlabel('reach length and point bar deposition & erosion')
legend('length','deposition','erosion')


%% Bivariate comparison matrix
figure(12)
clf
subplot(4,4,1) % plot bivariate comparison
hp =    plot(reach_length,depos_vol.p,'marker','.','linestyle','none');
alpha_marker(hp,0.1) % make markers semitransparent and monochrome
xlabel('length');ylabel('pt depo');
subplot(4,4,5)
hp =    plot(reach_length,depos_vol.v,'marker','.','linestyle','none');
alpha_marker(hp,0.1)
xlabel('length');ylabel('vt depo');
subplot(4,4,9)
hp =    plot(reach_length,erode_vol.p,'marker','.','linestyle','none');
alpha_marker(hp,0.1)
xlabel('length');ylabel('pt erode');
subplot(4,4,13)
hp =    plot(reach_length,erode_vol.v,'marker','.','linestyle','none');
alpha_marker(hp,0.1)
xlabel('length');ylabel('vt erode');

subplot(4,4,6)
hp =    plot(depos_vol.p,depos_vol.v,'marker','.','linestyle','none');
alpha_marker(hp,0.1)
xlabel('pt depo');ylabel('vt depo');
alpha_marker(hp,0.1)
subplot(4,4,10)
hp =    plot(depos_vol.p,erode_vol.p,'marker','.','linestyle','none');
alpha_marker(hp,0.1)
xlabel('pt depo');ylabel('pt erode');
subplot(4,4,14)
hp =    plot(depos_vol.p,erode_vol.v,'marker','.','linestyle','none');
alpha_marker(hp,0.1)
xlabel('pt depo');ylabel('vt erode');

subplot(4,4,11)
hp =    plot(depos_vol.v,erode_vol.p,'marker','.','linestyle','none');
alpha_marker(hp,0.1)
xlabel('vt depo');ylabel('pt erode');
subplot(4,4,15)
hp =    plot(depos_vol.v,erode_vol.v,'marker','.','linestyle','none');
alpha_marker(hp,0.1)
xlabel('vt depo');ylabel('vt erode');

subplot(4,4,16)
hp =    plot(erode_vol.p,erode_vol.v,'marker','.','linestyle','none');
alpha_marker(hp,0.1)
xlabel('pt erode');ylabel('vt erode');
alpha_marker(hp,0.1)

%% DEBUG test crosscorr function and plot crosscorr outputs
% X1 = reach_length/mean(reach_length);
% X2 = erode_vol.p  /mean(erode_vol.p);
% 
% [XCF, lags, bounds] = crosscorr(X1, X2);
% figure(13)
% clf
% subplot(2,1,1)
% crosscorr(reach_length,erode_vol.p)
% subplot(2,1,2)
% h =plot(lags,XCF);
% set(h,'color','k','marker','o','markerfacecolor','k','linestyle','none')
% set(gca,'xtick',[-10 0 10],'xgrid','on','xminorgrid','on')
% ylim([0 1])
% xlim([-10 10])
%% Cross correlation plots Ch4 fig27 
names = {'both reaches ln','pt dep vol','vt dep vol',...
                           'pt erosion','vt erosion'};
fig_names = {'Reach Length','Point Bar','Flood','Point Bar','Flood'};
data = [reach_length',depos_vol.p,depos_vol.v,erode_vol.p,erode_vol.v];

hcc = figure(14);
clf
n=numel(names);
for xx = 1:n    % plot all the cross correlations in a lower triangle
    plt_num = n * (xx - 1);                 % subplot index
    for yy = 1:xx
        plt_num = plt_num + 1;
        hf = subplot( n,n, plt_num);
        [XCF, lags, bounds] = crosscorr( data(:,xx), data(:,yy), 10);
        h = plot( lags, XCF );
        hf.YLim = [0 1];
        hf.XLim = [-10 10];
        set( h, 'color', 'k', 'marker', 'o', 'markerfacecolor', 'k',...
            'markersize', 3, 'linestyle', 'none' )
        set( hf, 'xtick', [-10:2:10], 'xgrid', 'on', 'ygrid', 'on' )
        hf.YTick = [ 0:.2:1 ];
        hf.GridAlpha = .3;
        title([names{xx},' lagged by ',names{yy}])
        % add labels to diagonal subplots
        if isequal(xx,yy)
            text( -8,.1, fig_names{xx} )
        end
    end
end

% add legend
hf = subplot(n,n,5);
[XCF, lags, bounds] = crosscorr( data(:,4), data(:,4), 10 );
h = plot( lags, XCF );
hf.YLim = [0 1];
hf.XLim = [-10 10];
set( h, 'color', 'k', 'marker', 'o', 'markerfacecolor', 'k',...
    'markersize', 3, 'linestyle', 'none' )
set( hf, 'xtick', [-10:2:10], 'xgrid', 'on', 'ygrid', 'on' )
hf.GridAlpha = .3;

% remove space between subplots
hccp = hcc.Children;
for i = 1:numel(hccp)
    hccp(i).XTickLabel = [];
    hccp(i).YTickLabel = [];
    hccp(i).Title = [];
    hccp(i).Position([ 3 4 ]) = [ .16282  .17267 ];
end

% format legend & add back labels
hf.XTickLabel = { -10; '';'';'';'';0;'';'';'';''; 10}; % mostly blank
hf.YTick = [0:.2:1];
hf.YTickLabel = {0 ; '';'R^2';'';''; 1};
text( -9.2, 0.88, 'KEY', 'fontsize', 14 )
text( -4  , 0.25, 'Title' )
text( -10 , 0.1 , '(Row & Column)' )
xlabel('Lag time steps')
hf.Position([ 1 2 ]) = [ 0.72 0.6 ];              % reposition
text( -40,  0.7, 'Deposition', 'fontsize', 14)    % other labels
text( 0  , -1.3, 'Erosion'   , 'fontsize', 14)

% set gray background for deposition related subplots
set(hccp([ 4,5,8,9,11:15 ]), 'color', 1-.15*[1 1 1] )

hcc.Units = 'inches';
hcc.Position = [2 2 6.25 4.2];
set( hcc, 'color', 'w', 'PaperPosition', hcc.Position );
drawnow
fname = [file, ' length deposition erosion crosscorr.png' ];
print('-painters', '-dpng', '-r600', fname)
% savefig(14,[files{i},' length deposition erosion crosscorr.fig'])
% plot2svg([files{i},' length deposition erosion crosscorr.svg'])

%% violins
hfv = figure(33);
clf
hdv = distributionPlot( data( 4000:6845, 2:end ), 'color', [.3 .3 .3], ...
    'showMM', 6, 'distWidth', 1.2, 'xnames',...
    {'Pt Bar Dep','Flood Dep','Pt Bar Erosion','Flood Erosion'});
box on

hold on % add markers for mean
hmv = plot(hdv{3}.XTick,mean(data( 4000:6845, 2:end )),...
    'xw','MarkerSize',12,'linewidth',2);
hold off
set(hdv{2}(1), 'color', [1 1 1 ], 'linewidth', 2); % change color of the x
set(hdv{2}(:), 'color', [0 0 0 ]);                 % set quantiles to black

ylabel('Volume of sediment (full pixels)')

hfv.Children.Position = [.1 .09 .865 .875];
hfv.Units = 'inches';
hfv.Position = [2 6.5 5.8 6];
set( hfv, 'color', 'w', 'PaperPosition', hfv.Position );
drawnow
fname = [file,' both erosion deposition 4000 6845.png'];
print('-painters', '-dpng', '-r600', fname)
% plot2svg([files{i},' both erosion deposition 4000 6845.svg'])

end
