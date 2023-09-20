function alpha_marker(hPlot,alph)
% Purpose   set alpha value for marker face color on figure
% Author    Tobias Hasse tobiack@udel.edu
% Date      Summer 2019
drawnow;  % drawnow is needed or the handles dont exist for next lines
view(0,90) % for plot3 non isometric view
hMarkers = hPlot.MarkerHandle;
hMarkers.EdgeColorData = uint8(255*[0;0;0;alph]);
hMarkers.FaceColorData = uint8(255*[1;1;1;  0]);
end
