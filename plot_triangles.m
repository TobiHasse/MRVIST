function plot_triangles(fname, time_step_years, ...
    pb_age_dist, vd_age_dist, pt_bar_dist, vert_dist, pb_ftt, vd_ftt)
% Purpose:  Take the storage time distributions from all model steps and 
%           create figures showing the distributions as in the 
%           figures 28-39 of Chapter 4
%           code excerpted from plot_figs_new.m
% Inputs    fname               input file name
%           times_step_years    the depositional model time step
%           flood (vertical, v) deposit and
%           point bar (pt_bar, pb) deposit distributions for
%           age     of stored material, 
%           dist    storage time of eroding material and
%           ftt     life expectancy or Forward Transit Time
% Outputs   figure .png files with the data file name prepended
%           figure .svg files " " " "
%           images of the indexed data (using imwrite) to import into the
%                  svg files for higher qualit yimages
% dependancies  cmapsPUB        customized colormaps 
%               set_ax_full     adjust axis to maximize 5.8 inch wide paper
%               set_ax_zoom     format zoomed in figure
%
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     August, 2019, edited November 2021

%% get the colormap data
cdata = cmapsPUB('topo',0);
% select the zoomed in area of the plot
ax_zoom = [.5 750.5 3999.5 4749.5 ]; % powerpoint, png & svg output
ax_zm = [ 1 751 4000 4750 ];         % imwrite output

%% AGE DISTRIBUTION figures point bar then flood deposits
plt = fliplr(vd_age_dist);               % copy to temporary variable
mxp = max(plt(:));                       % needed to index data for imwrite
%
hva = figure(11);
    clf
    imagesc(plt)
    h = colorbar('location','east','direction','normal'); 
    colormap(flipud(cdata(2:end,:)))        % colors scaled to the data
    title( sprintf( 'Age dist of overbank material. Total amount %.2e ',...
        sum( plt(:) ) ) )
    ylabel(h,'Volume of sediment remaining (full pixels)','fontsize',12)
    xlabel('Age of floodplain sediment (kyr)')

    drawnow
    set_ax_full(hva,'upper',time_step_years) % format figure
    outfile = [fname,' Age _ flood.png'];
    % print to file: pretty good output, reliable
    print('-painters', '-dpng', '-r600', outfile)
    % save as svg for editing and annotation, poor raster image quality
    saveas(hva, [outfile(1:end-4),'.svg'], 'svg')
    % imwrite the data to import high resolution images into a vector 
    % graphics editor.  
    % 254 to make 255 indexes to match length of colormap
    % + 1 so that 0 is white, greater than 0 is on the color map 
    % XREF indicates this file should be an external reference to the svg
    imwrite( ceil(254*plt./mxp) + 1, ...
        flipud(cdata(2:end,:)),[outfile(1:end-4),' XREF.png'])

    set_ax_zoom(hva, time_step_years, ax_zoom) % format figure
    outfile = [fname,' zoom Age _ flood.png'];
    print('-painters', '-dpng', '-r600', outfile)
    saveas(hva, [outfile(1:end-4),'.svg'], 'svg')
    
    imwrite(ceil(254*plt( ax_zm(3):ax_zm(4) , ax_zm(1):ax_zm(2) )./mxp)...
        + 1, flipud(cdata(2:end,:)),[outfile(1:end-4),' XREF .png'])

    close(hva)

%% corollary figure for flood deposits
plt = fliplr(pb_age_dist);
% mxp = max(plt(:)); % comment out to scale colors to vd_age_dist above
plt( plt > mxp ) = mxp; % 
 
hpa = figure(12);
    clf
    imagesc(fliplr(pb_age_dist))
    h=colorbar('location','east','direction','normal');
    % find how long the color bar should be so the colorful parts match
    c_extend = round( 255 * max( pb_age_dist(:) ) / ...
        max( vd_age_dist(:) ) - 255 );      % approximately 412
    cbrn=repmat(cdata(2,:),c_extend,1);     % match colors with flood fig
    colormap(flipud([cbrn;cdata(2:end,:)]))
    %     colormap(flipud(cdata(2:end,:)))  % colors scaled to this data
    title( sprintf( 'Age dist of pt bar material. Total amount %.2e ',...
        sum( pb_age_dist(:) ) ) )
    ylabel(h,'Volume of sediment remaining (full pixels)','fontsize',12)
    xlabel('Age of point bar sediment (kyr)')
    
    drawnow               % needed to create figure.Children for formatting
    set_ax_full(hpa,'upper',time_step_years) % format figure
    outfile = [fname,' Age _ point bar.png'];
    print('-painters', '-dpng', '-r600', outfile)
    saveas(hpa, [outfile(1:end-4),'.svg'], 'svg')
    imwrite( ceil(254*plt./mxp) + 1, ...
        flipud(cdata(2:end,:)),[outfile(1:end-4),' XREF.png'])
        
    set_ax_zoom(hpa, time_step_years, ax_zoom) % format figure
    outfile = [fname,' zoom Age _ point bar.png'];
    print('-painters', '-dpng', '-r600', outfile)
    saveas(hpa, [outfile(1:end-4),'.svg'], 'svg')
    imwrite(ceil(254*plt( ax_zm(3):ax_zm(4) , ax_zm(1):ax_zm(2) )./mxp)...
        + 1, flipud(cdata(2:end,:)),[outfile(1:end-4),' XREF .png'])

    close(hpa)
%% STORAGE TIME DISTRIBUTIONS point bar then flood deposits
plt = fliplr(pt_bar_dist);
mxp = max(plt(:));

hps = figure(21);
    clf
    imagesc(plt)
    title('pt bar storage time')
    colormap(flipud(cdata(2:end,:)))
    h = colorbar('location','east','direction','normal'); 
    ylabel(h,'Fraction of eroding sediment older','fontsize',12)
    xlabel('Storage time of eroded point bar sediment (kyr)')

    drawnow
    set_ax_full(hps,'upper',time_step_years) % format figure
    outfile = [fname,' Storage _ point bar.png'];
    print('-painters', '-dpng', '-r600', outfile)
    saveas(hps, [outfile(1:end-4),'.svg'], 'svg')
    imwrite( ceil(254*plt./mxp) + 1, ...
        flipud(cdata(2:end,:)),[outfile(1:end-4),' XREF.png'])

    set_ax_zoom(hps, time_step_years, ax_zoom) % format figure
    outfile = [fname,' zoom Storage _ point bar.png'];
    print('-painters', '-dpng', '-r600', outfile)
    saveas(hps, [outfile(1:end-4),'.svg'], 'svg')
    imwrite(ceil(254*plt( ax_zm(3):ax_zm(4) , ax_zm(1):ax_zm(2) )./mxp)...
        + 1, flipud(cdata(2:end,:)),[outfile(1:end-4),' XREF .png'])

    close(hps)
%
plt = fliplr(vert_dist);
mxp = max(plt(:));

hvs = figure(22);
    imagesc(plt)
    title('overbank storage time')
    colormap(flipud(cdata(2:end,:)))
    h= colorbar('location','east','direction','normal'); 
    ylabel(h,'Fraction of eroding sediment older','fontsize',12)
    xlabel('Storage time of eroded flood sediment (kyr)')

    drawnow
    set_ax_full(hvs,'upper',time_step_years) % format figure
    outfile = [fname,' Storage _ flood.png'];
    print('-painters', '-dpng', '-r600', outfile)
    saveas(hvs, [outfile(1:end-4),'.svg'], 'svg')
    imwrite( ceil(254*plt./mxp) + 1, ...
        flipud(cdata(2:end,:)),[outfile(1:end-4),' XREF.png'])

    set_ax_zoom(hvs, time_step_years, ax_zoom) % format figure
    outfile = [fname,' zoom Storage _ flood.png'];
    print('-painters', '-dpng', '-r600', outfile)
    saveas(hvs, [outfile(1:end-4),'.svg'], 'svg')
    imwrite(ceil(254*plt( ax_zm(3):ax_zm(4) , ax_zm(1):ax_zm(2) )./mxp)...
        + 1, flipud(cdata(2:end,:)),[outfile(1:end-4),' XREF .png'])
    
    close(hvs)
%%  FORWARD TRANSIT TIME __ LIFE EXPECTANCY point bar then flood
plt = pb_ftt;
mxp = max(plt(:));

hpl = figure(61);
    clf
    imagesc(pb_ftt)
    title('Point Bar Life Expectancy')
    colormap(flipud(cdata(2:end,:)))
    h= colorbar('location','east','direction','normal'); 
    ylabel(h,'Fraction of sediment cohort remaining','fontsize',12)
    xlabel('Life expectancy of sediment cohorts (ky)')

    drawnow
    set_ax_full(hpl,'lower',time_step_years) % format figure
    outfile = [fname,' Life expectancy _ point bar.png'];
    print('-painters', '-dpng', '-r600', outfile)
    saveas(hpl, [outfile(1:end-4),'.svg'], 'svg')
    imwrite( ceil(254*plt./mxp) + 1, ...
        flipud(cdata(2:end,:)),[outfile(1:end-4),' XREF.png'])
    
    set_ax_zoom(hpl, time_step_years, ax_zoom) % format figure
    outfile = [fname,' zoom Life expectancy _ point bar.png'];
    print('-painters', '-dpng', '-r600', outfile)
    saveas(hpl, [outfile(1:end-4),'.svg'], 'svg')
    imwrite(ceil(254*plt( ax_zm(3):ax_zm(4) , ax_zm(1):ax_zm(2) )./mxp)...
        + 1, flipud(cdata(2:end,:)),[outfile(1:end-4),' XREF .png'])

    close(hpl)
%
plt = vd_ftt;
mxp = max(plt(:));

hvl = figure(62);
    clf
    imagesc(vd_ftt)
    title('Flood Life Expectancy')
    colormap(flipud(cdata(2:end,:)))
    h= colorbar('location','east','direction','normal'); 
    ylabel(h,'Fraction of sediment cohort remaining','fontsize',12)
    xlabel('Life expectancy of sediment cohorts (ky)')

    drawnow
    set_ax_full(hvl,'lower',time_step_years) % format figure
    outfile = [fname,' Life expectancy _ flood.png'];
    print('-painters', '-dpng', '-r600', outfile)
    saveas(hvl, [outfile(1:end-4),'.svg'], 'svg')
    imwrite( ceil(254*plt./mxp) + 1, ...
        flipud(cdata(2:end,:)),[outfile(1:end-4),' XREF.png'])
 
    set_ax_zoom(hvl, time_step_years, ax_zoom) % format figure
    outfile = [fname,' zoom Life expectancy _ flood.png'];
    print('-painters', '-dpng', '-r600', outfile)
    saveas(hvl, [[outfile(1:end-4),'.svg'],'.svg'], 'svg')
    imwrite(ceil(254*plt( ax_zm(3):ax_zm(4) , ax_zm(1):ax_zm(2) )./mxp)...
        + 1, flipud(cdata(2:end,:)),[outfile(1:end-4),' XREF .png'])

    
end                                     % function plot_tirangles

function set_ax_full(hf, cBarPosition, dt)
% Purpose:  set the axis properties for publication of the triangle
%           distribution plots for the 5.8 inch wide dissertation format
% Inputs    hf              the figure handle
%           cBarPosition    'upper' or 'lower' the color bar position
%           dt              the depositional time step in years
%
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     August, 2019, edited November 2021

% adjust caxis to make white background make caxis minimum less than data
colordata = colormap;               
ca = caxis;                     
caxis( [ -ca(2) / ( size( colordata, 1 ) -1 ) + 10^-10 ca(2) ] ) 
  

ylabel('Model time step (kyr)')      % add the y label, same for all plots
hf.Units = 'inches';                    % set the units and position
hf.Position = [.5 .5 5.8 5.8];
set( hf, 'color', 'w', 'PaperPosition', hf.Position );

ha = hf.Children(2);                    % get the axis handle
ha.Title = [];                          % delete the figure title

ha.XTickLabel = ha.XTick * dt / 1000;   % tick labels in thousands of years
ha.YTickLabel = ha.YTick * dt / 1000;
ha.LineWidth = 1.5;

p1 = .07;                               % lower left coordinate of axes
p2 = .925;                              % axes size (square for square fig)
ha.Position = [p1 p1 p2 p2];            % maximize plot area
ha.XRuler.TickLabelGapOffset = -3;      % tighten up axis labels
ha.YRuler.TickLabelGapOffset = -2;
ha.YLabel.VerticalAlignment = 'baseline';
ha.XLabel.VerticalAlignment = 'cap';


hc = hf.Children(1);                    % colorbar handle and positioning 
switch cBarPosition  % normal east position [ 0.9088 0.1048 0.0479 0.8519 ]
    case 'upper'  
        hc.Position = [0.9288    0.2048    0.0479    0.7519]; 
    case 'lower'
        hc.Position = [0.9288    0.1048    0.0479    0.7519]; 
    otherwise 
        disp('please choose colorbar position upper or lower')
end
end                                     % function set_ax_full

function set_ax_zoom(hf, dt, ax_lims)
% Purpose:  set the axis properties for publication of the triangle
%           distribution plots for the 5.8 inch wide dissertation format
% Inputs    hf          the figure handle
%           dt          the depositional time step in years
%           ax_lims     limits of the zoomed in figure
%
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     August, 2019, edited November 2021

% adjust caxis to make white background make caxis minimum less than data
colordata = colormap;               
ca = caxis;                     
caxis( [ -ca(2) / ( size( colordata, 1 ) -1 ) + 10^-10 ca(2) ] ) 
  

ylabel('Model time step (kyr)')      % add the y label, same for all plots
hf.Units = 'inches';                    % set the units and position
hf.Position = [.5 .5 5.8 6.4];
set( hf, 'color', 'w', 'PaperPosition', hf.Position );

ha = hf.Children(2);                    % get the axis handle
ha.Title = [];                          % delete the figure title
ha.XLim = ax_lims(1:2);
ha.YLim = ax_lims(3:4);

ha.XTickLabel = ha.XTick * dt / 1000;   % tick labels in thousands of years
ha.YTickLabel = ha.YTick * dt / 1000;
ha.LineWidth = 1.5;

% p1 = .07;                             % lower left coordinate of axes
% p2 = .855;                            % axes size (square for square fig)
% ha.Position = [p1 p1 p2 p2];          % maximize plot area
ha.XRuler.TickLabelGapOffset = -3;      % tighten up axis labels
ha.YRuler.TickLabelGapOffset = -1;
ha.YLabel.VerticalAlignment = 'baseline';
ha.XLabel.VerticalAlignment = 'cap';


hc = hf.Children(1);                    % colorbar handle and positioning 
hc.Visible = 'off';                     % turn off old colorbar
h = colorbar('location','southoutside','direction','reverse'); 
ylabel(h,hc.Label.String,'fontsize',11) % use old label for new colorbar
% Position colorbar and figure
% 'normalized'
% h.Position = [.07 .068 .925 .028]; % 0.0826    0.1103    0.8998    0.0408
% ha.Position = [.07 .156 .925 .837];% 0.0700    0.2288    0.9250    0.7662
% % ha.Position = [.07 .18 .925 .789]; % 5.8x6.8 inch
% 'pixels'
h.Units = 'pixels';
h.Position = [40 42 515 20];
ha.Units = 'pixels';
ha.Position = [40 97 515 515];

%%
end                                     % function set_ax_zoom

