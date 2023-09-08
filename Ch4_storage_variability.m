function Ch4_storage_variability(pt_bar_dist,vert_dist,dt)
% Purpose:  Show the variability through time of certain statistical
%           measures of storage time including the time series of median 
%           storage time accross 2800 time steps beginning
%           with time step 4000 and plot this in 7 subplot panels as if it
%           were a monthly calendar with 7 weeks
%           Make a figure showing the frequency distribution of the above
%           time series and others
%           Show the autocorrelation of the above time series
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     2018-2019, edited October 2021

% public function dependencies
%           distributionPlot from the MATLAB file Exchange FEX
%           http://www.mathworks.com/matlabcentral/fileexchange/
%           37105-plot-spread-points-beeswarm-plot

% private functions dependencies
%          pdfcdfpctile.m
%          Ch4fig22_calendar_plot
%          Ch4fig23_violins_storage_time 
%          Ch4fig24_autocorr
%          Storage time distributions from the complete model _storage_ run

px = [1:length(vert_dist)-1]'*dt;
%% DEBUG
% load in data used for dissertation figure
% cd('C:\Users\User\Documents\MATLAB\MacMatlab\variability')
% 
% load('TRH 205k storage time variability.mat','ma');
% maTRA = ma;


%% extract median from storage distributions

% correct method, but results are 30 years younger than dissertation
% (dt = 30) because the last column in the storage time distributions are
% all zeros.  The second to last colum represents sediment stored for 30
% years (not 60) and therefore the last column should be excluded when 
% finding the median.  
% pdfcdfpctile.m finds the index of the value and then scales it by 
% dt to obtain the age at that percentile
[p, c, ma] = pdfcdfpctile( fliplr( vert_dist( : , 1:end-1 ) + ...
    pt_bar_dist( : , 1:end-1 ) ), dt, 0.50);

% make calendar
figno = 58;
Ch4fig22_calendar_plot( ma, dt, figno )

%% matches dissertation fig 22 and includes the error that makes ages
% 30 years (dt) older than they should be (notice end vs end -1 

% [p, c, ma] = pdfcdfpctile( fliplr( vert_dist( : , 1:end   ) + ...
%                                  pt_bar_dist( : , 1:end   ) ), dt, 0.50);
% %
% figno = 59;
% Ch4fig22_calendar_plot(ma, dt, figno)

%% DEBUG compare both methods
% plts =5;
% figure(7)
% subplot(plts,1,1)
% % plot(ma+30-maTRA(1:length(ma))) %           +30
% plot(ma-maTRA(1:length(ma))) %
% subplot(plts,1,2)
% plot(ma,maTRA(1:length(ma)))
% subplot(plts,1,3)
% plot(ma)
% subplot(plts,1,4)
% plot(maTRA(1:length(ma)))
% subplot(plts,1,5)
% plot(c(1:end,1:5))

%% Violin plots of storage time Ch4fig23
% get the PDF's CDF's, and 70th percentile values
[p, c, pa70] = pdfcdfpctile( fliplr( vert_dist( : , 1:end-1 ) + ...
    pt_bar_dist( : , 1:end-1 ) ), dt, 0.70);
% compute the residence times from the PDF distribution
avea =sum(bsxfun(@times,p,px'),2)./sum(p,2);

% make violin plots of storage times for the portion of the simulation
% analyzed in the dissertation, from 120 to 205 kyr (4000:6845) HARD CODED
figno = 60;

Ch4fig23_violins_storage_time( avea(4000:6845), pa70(4000:6845), ...
    ma(4000:6845), dt, figno )

% make autocorrelation plot
% slightly lower autocorrelation when using just a portion of the data
% compared with using the entire vectors
figno = 49;
Ch4fig24_autocorr( ma(4000:6845), pa70(4000:6845), avea(4000:6845),...
    dt, figno )

%% DEBUG
% load in data used for dissertation figure
% cd('C:\Users\User\Documents\MATLAB\MacMatlab\variability')
% 
% load('TRH 205k storage time variability.mat','ma');
% figno = 48;
% Ch4fig24_autocorr(ma, pa70, avea, dt, figno)


end % function Ch4_storage_variability


function Ch4fig24_autocorr(ma, pa70, avea, dt, figno)
% Purpose   This function will make the autocorrelation figure for an input
%           vector(s)
% Author:   Tobias Hasse tobiack@udel.edu

hf = figure(figno);
clf
    hf.Units = 'inches';    % set figure units to inches
    set(hf,'color','w',...
        'position',[.5 1 5.8 4]) % inches

[acf lags bounds] = autocorr(ma);
lags = lags * dt;
px = lags;
hp = stem(px,acf,'color','k','markerfacecolor','k');
ylim([ 0 1 ])
hp = hf.Children;   % Get plot handle from figure
hp.Position = [ .1 .12 .87 .82];  % .13 .11434 .775 .81066
hp.XTick = px(1:2:end);
hp.YTick = [0 bounds(1) 0.5 1];
hp.YGrid = 'on';
hp.GridAlpha = .5;

xlabel('Lag (years)')
hy = ylabel( 'Sample Autocorrelation');
hy.Position(1) = -35; % set position of ylabel

title('Autocorrelation of Median storage time analyzed')
set( hf, 'color', 'w', 'PaperPosition', hf.Position );
print('-painters', '-dpng', '-r600', 'median autocorr.png')

% subplot(3,1,1)
% autocorr(ma)
% set(gca,'xticklabel',[0:2:20]*dt)
% xlabel('Lag (years)')
% title('Autocorrelation of Median storage time')
% subplot(3,1,2)
% autocorr(avea)
% set(gca,'xticklabel',[0:2:20]*dt)
% xlabel('Lag (years)')
% title('Autocorrelation of Average storage time')
% subplot(3,1,3)
% autocorr(pa70)
% set(gca,'xticklabel',[0:2:20]*dt)
% xlabel('Lag (years)')
% title('Autocorrelation of 70th percentile storage time')

end % function Ch4fig24_autocorr


function Ch4fig23_violins_storage_time(aa,p7,md,dt,figno)
% Purpose   This function will make violin plots for 3 distributions and
%           overlay multiple plots to show the structure of each in the
%           most appropriate y-axis scaling and save a file in the current
%           folder. 
%           This code is excerpted from calendar_plot.m and formerly from
%           an_howardPC.m from the Windows 8.1 machine and
%           creates figure 23 from Chapter 4

%           An earlier version of this code used a customized version of
%           distributionPlot.m with the option showwMM 7 added to
%           distributionPlot to combine options 2 and 6.  The code here
%           hase been updated to add option 2 (plot the mean) after the
%           rest of the distribution is created by distributionPlot.  This
%           is accomplished by using hold on, plot(stuff) hold off and plot
%           formatting code is slightly adjusted
%           legacy code is commented out with the note: "% showMM 7"
% Author:   Tobias Hasse tobiack@udel.edu

% Inputs    ma              a vector
%           figno           figure number
%           dt              the time step in years

%% DEBUG select a subset of the input data
% aa = avea(4000:6845);
% p7 = pa70(4000:6845);
% md = ma  (4000:6845);
%%
hf = figure(figno); 
    clf
    hf.Units = 'inches';    % set figure units to inches
    set(hf,'color','w',...
        'position',[.5 1 6.25 6.25]) % inches
        %  'position',[245 1000 600 600]) % 752 564 % pixels
%
h1 = subplot(3,1,1);    % all 3
    hd1 = distributionPlot([aa,p7,md],...
        'color',[.4 .4 .4],'showMM',6,'histOpt',2,'distWidth',.9,...
        'xnames',{'Residence Time','70^t^h Percentile','Median'}); 
    % showMM option 7 = a combo of options 2 & 6 customized by Tobias Hasse
    hy = ylabel('Storage time (years)');
    hy.Position(1)=.25; % set position of ylabel
    box on
    %     h1.YLim =([0 21000]); % Hard coded option
    h1.YLim =([min([aa;p7;md])-dt max([aa;p7;md])]); % temporary y limits
    h1.XLim =[.4 3.8];          % sets padding space next to violins
    % figure subheading
    text( h1.XLim(1) + 0.03 * diff( h1.XLim ), ...
          h1.YLim(2) - 0.00 * diff( h1.YLim ), 'A)', 'fontsize', 16)

    hold on
    hm1 = plot(h1,hd1{3}.XTick,mean([aa,p7,md]),...
        'xw','MarkerSize',12,'linewidth',2);
    hold off

h2 = subplot(3,1,2); % two violins
    hd2 =    distributionPlot([p7,md],...
        'color',[.4 .4 .4],'showMM',6,'histOpt',2);
    h2.XTickLabel = [];
    box on
    %     h2.YLim =([250 3725]);
    h2.YLim =([min([p7;md])-dt max([p7;md])]);
    h2.XLim =[.4 2.7];
    text( h2.XLim(1) + 0.05 * diff( h2.XLim ), ...
          h2.YLim(2) - 0.0 * diff( h2.YLim ), 'B)', 'fontsize', 16)

    hold on
    hm2 = plot(h2,hd2{3}.XTick,mean([p7,md]),...
        'xw','MarkerSize',12,'linewidth',2);
    hold off

h3 = subplot(3,1,3); % just the median violin
    hd3 =    distributionPlot(md,...
        'color',[.4 .4 .4],'showMM',6,'histOpt',2);
    h3.XTickLabel = [];
    box on
    %     h3.YLim = ([300 780]);
    h3.TickLength  =  h3.TickLength*2;
    h3.YLim = ([min(md)-dt max(md)+dt]);
    h3.XLim =[.3 1.75];
    text( h3.XLim(1) + 0.1 * diff( h3.XLim ), ...
          h3.YLim(2) - 0.0 * diff( h3.YLim ), 'C)', 'fontsize', 16)

    hold on
    hm3 = plot(h3,hd3{3}.XTick,mean(md),...
        'xw','MarkerSize',12,'linewidth',2);
    hold off

% format the x marking the mean of each distribution 
% showMM 7
% ( required for customized distributionPlot and showMM 7)
% set(hd1{2}(1), 'color', [1 1 1 ], 'linewidth', 2);  
% set(hd2{2}(1), 'color', [1 1 1 ], 'linewidth', 2);  
% set(hd3{2}(1), 'color', [1 1 1 ], 'linewidth', 2);  
% format quantile lines of each distribution (standard distributionPlot)
set(hd1{2}(:), 'color', [0 0 0 ]);                  
set(hd2{2}(:), 'color', [0 0 0 ]);                  
set(hd3{2}(:), 'color', [0 0 0 ]);                  
% obtain and create yticks which label important quantiles from the median
% distribution violin
tmp = cell2mat( get( hd3{2}(:),'YData'));
yticks3 = round( unique( sort( [ reshape( tmp, 1,numel(tmp) ),...
    hm3.YData, max( md ), min( md ) ] ) ) , 0 );
% showMM 7
% yticks3 = round(unique(sort([hd3{2}(:).YData,max(md),min(md)])),0);
posi3 =[.59 .54 .275 .345];
set(h3,'position',posi3,...
    'ytick',     yticks3)
% the following ytick vector is for the dissertation which includes two
% errors: 1) off by 30 years 2) includes the 'terrace', that is, the
% erosion of the floodplain which has never previously been visited by the
% channel
% set(h3,'ytick',     [330 450 480 497 540 750])

% obtain and create yticks with label inportant quantiles from the 70th
% percentile violin but exclude the prior tick marks from yticks3
tmp =cell2mat(get(hd2{2}(:),'YData'));
yticks2 = setdiff( round( unique( sort( [ reshape( tmp, 1,numel(tmp) ),...
    hm2.YData, h2.YLim(2), min(p7)] ) ),0), yticks3);
% showMM 7
% yticks2 = setdiff( round( unique( sort( [ hd2{2}(:).YData,...
%     h2.YLim(2), min(p7) ] ) ), 0), yticks3);

% posi2 = [.31 .28 .575 .625];  % old positioning
posi2 = [.357 .28 .528 .625];   % new positioning
set(h2,'position',posi2,...
    'ytick',yticks2)
% the following ytick vector is for the dissertation which includes two
% errors: 1) off by 30 years 2) includes the 'terrace'
%     set(gca,'position',posi2,'ytick',[570 840 930 965 1050 3600])

% obtain uniqe yticks for the residence time
tmp =cell2mat(get(hd1{2}(:),'YData'));
yticks1 = setdiff( round( unique( sort( [ reshape( tmp, 1,numel(tmp) ),...
    hm1.YData, h1.YLim(2), min(aa)] ) ), 3, 'significant'), yticks2);
% showMM 7
% yticks1 = setdiff( round( unique( sort( [ hd1{2}(:).YData,...
%     h1.YLim(2), min(aa) ] ) ), 3, 'significant'), yticks2);
yticks1 = setdiff( yticks1 , yticks3 );
posi = [.13 .11 .775 .815];         % new position
%     posi = [.13 .31 .775 .315];   % old position
set(h1,'position',posi,...
    'ytick',    yticks1)
% the following ytick vector is for the dissertation which includes two
% errors: 1) off by 30 years 2) includes the 'terrace'
% set(h1,      'ytick',     [539 1210 1830 2490 2940 20300],...
%         'yticklabel',     [539 1210 1830 2490 2940 20300])

[log10(h1.YLim(2)) log10(h2.YLim(2)) log10(h3.YLim(2)) ];
tmp        = 10^(  floor( log10( h1.YLim(2)))-1);
h1.YLim(2) = ceil( h1.YLim(2)/tmp + 4 * h1.YLim(2) / 100 / tmp)*tmp;
tmp        = 10^(  floor( log10( h2.YLim(2)))-1);
h2.YLim(2) = ceil( h2.YLim(2)/tmp + 4 * h2.YLim(2) / 100 / tmp)*tmp;
tmp        = 10^(  floor( log10( h3.YLim(2)))-1);
h3.YLim(2) = ceil( h3.YLim(2)/tmp + 4 * h3.YLim(2) / 100 / tmp)*tmp;
    
% set( hf, 'color', 'w', 'PaperPosition', hf.Position );
% drawnow
% outfile = 'violin storage times both reaches.png';
% print('-painters', '-dpng', '-r600', outfile)
    
%%
end % function violins_storage_time


function Ch4fig22_calendar_plot( ma, dt, figno )
% Purpose   This function will make the calendar plot of median storage
%           time and save it to file in the current folder
%           This code is excerpted from calendar_plot.m and formerly from
%           an_howardPC.m from the Windows 8.1 machine and
%           creates figure 22 from Chapter 4
% Author:   Tobias Hasse tobiack@udel.edu

% Inputs    ma              a vector
%           figno           figure number
%           dt              the time step in years

px=[1:length(ma)]'*dt;      % x coordinate for plots

hf = figure(figno);         % get figure handle

hf.Units = 'inches';        % set figure units and positoin
hf.Position = [1 1 5.8 8.31];

nplots = 7;                 % number of 'weeks' in the calendar (sub)plots

for i = 1:nplots
    i_s = 4000 + 400 * (i-1);   % start index
    i_e = 4000 + 400 *  i;      % end   index
    h = subplot( nplots, 1, i );
    plot( px( i_s:i_e ), ma( i_s:i_e ), 'color', 'k');
    posi = h.Position;      % adjust subplot positions
    posi = posi + [ -.04 0 .11 0];
    posi(4)=.1212;
    h.Position = posi;
    h.YLim = [ min( ma( 4000:6845 )), max( ma( 4000:6845 ) ) ];
    yl = h.YLim;
    %         ylim(yl)
    %         yl = [ min( ma( 4000:6845 )), max( ma( 4000:6845 ) ) ];
    %         ylim(yl)
    % annotate lower left corner of each subplot with starting age
    text( h.XLim(1), h.YLim(1) + 0.1 * diff( h.YLim ),...
        strcat('\leftarrow', num2str( h.XLim( 1 ) / 1000 ), ' kyr'))
    
    if isequal(i,1)
        title(strcat('Calander plot of median storage time variability',...
            ' (all sediment) 120 - 204 kyr' ) )
    end
    % label of average median storage time
    lab_mu = round( mean( ma( 4000:6845 )), 0);
    % y-axis label on middle subplot because there is a lot of text
    if isequal(i,ceil(nplots/2))
        ylabel(strcat('Median age of eroding sediment (years) (\mu = ',...
            num2str(lab_mu),')'))
    end
    
    if i<nplots % format all but the last subplot
        box on
        set(h,'xtick',[])
    else        % x-axis labels assuming kyr HARD CODED
        rng = ( h.XTick(end) - h.XTick(1)) / 1000;
        xLabs = 0:rng/( numel(h.XTick)-1):rng;
        set(h,'xticklabel', strcat('+',num2str([ xLabs ]'),' kyr'))
    end
    % format y-axis on each subplot
    set(h,'ytick',[yl(1),mean(ma(4000:6845)), yl(2)],...
        'ygrid','on','GridLineStyle','--','GridAlpha',1,...
        'color',1-.15*[1 1 1] * rem(i,2))
    try             % try to place the letter mu as a tick label
        set(h,'yticklabel',[ num2str(yl(1));'\mu';num2str(yl(2))])
    catch           % just use numbers
        set(h,'yticklabel',[yl(1),lab_mu,yl(2)])
    end
end
% set PaperPosition so that *.png file preserves figure size
set( hf, 'color', 'w', 'PaperPosition', hf.Position );
drawnow
outfile ='calendar plot all sediment both reaches median storage time.png';
print('-painters', '-dpng', '-r600', outfile)

%plot2svg('calendar plot all sediment both reaches median storage time.svg')
end % function Ch4fig22_calendar_plot
