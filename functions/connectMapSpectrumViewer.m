function connectMapSpectrumViewer(action, varargin)
%CONNECTMAPSPECTRUMVIEWER Add click and arrow-key voxel spectrum viewing.
%
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

switch lower(action)
    case 'setupfigure'
        figMap = varargin{1};
        spectrumFigureSize = varargin{2};
        set(figMap, 'Tag', 'SlowMRSIMapViewer', 'KeyPressFcn', @moveSelectedSpectrumVoxel);
        setappdata(figMap, 'SpectrumFigureSize', spectrumFigureSize);
        setappdata(figMap, 'SpectrumViewerData', loadViewerDataFromBase());
        setappdata(figMap, 'SpectrumViewerState', localEmptyViewerState());
    case 'enableclick'
        img = varargin{1};
        mapName = varargin{2};
        nz = varargin{3};
        img.UserData = struct('mapName', mapName, 'nz', nz);
        img.ButtonDownFcn = @showClickedSpectra;
        img.HitTest = 'on';
        img.PickableParts = 'all';
    otherwise
        error('Unknown map spectrum viewer action "%s".', action);
end
end

function showClickedSpectra(src, ~)
ax = ancestor(src, 'axes');
figMap = ancestor(src, 'figure');
clickPoint = ax.CurrentPoint;
n2 = round(clickPoint(1,1));
n1 = round(clickPoint(1,2));
clickInfo = src.UserData;

showVoxelSpectra(figMap, ax, n1, n2, clickInfo.nz, clickInfo.mapName);
end

function moveSelectedSpectrumVoxel(src, event)
figMap = findMapFigure(src);
if isempty(figMap) || ~isgraphics(figMap)
    return
end

selected = getappdata(figMap, 'SelectedSpectrumVoxel');
if isempty(selected) || ~isfield(selected, 'ax') || ~isgraphics(selected.ax)
    return
end

switch event.Key
    case 'leftarrow'
        selected.n2 = selected.n2 - 1;
    case 'rightarrow'
        selected.n2 = selected.n2 + 1;
    case 'uparrow'
        selected.n1 = selected.n1 - 1;
    case 'downarrow'
        selected.n1 = selected.n1 + 1;
    otherwise
        return
end

viewerData = getappdata(figMap, 'SpectrumViewerData');
dataSize = size(viewerData.ful);
selected.n1 = max(1, min(dataSize(1), selected.n1));
selected.n2 = max(1, min(dataSize(2), selected.n2));

showVoxelSpectra(figMap, selected.ax, selected.n1, selected.n2, selected.nz, selected.mapName);
end

function figMap = findMapFigure(src)
if isgraphics(src, 'figure') && strcmp(src.Tag, 'SlowMRSIMapViewer')
    figMap = src;
    return
end

figMap = findobj(0, 'Type', 'figure', 'Tag', 'SlowMRSIMapViewer');
if ~isempty(figMap)
    figMap = figMap(1);
end
end

function showVoxelSpectra(figMap, ax, n1, n2, nz, mapName)
viewerData = getappdata(figMap, 'SpectrumViewerData');
if isempty(viewerData) || ~isfield(viewerData, 'ful')
    viewerData = loadViewerDataFromBase();
    setappdata(figMap, 'SpectrumViewerData', viewerData);
end

dataSize = size(viewerData.ful);
if n1 < 1 || n2 < 1 || nz < 1 || ...
        n1 > dataSize(1) || n2 > dataSize(2) || nz > dataSize(3)
    return
end

updateSelectionMarker(figMap, ax, n1, n2);
setappdata(figMap, 'SelectedSpectrumVoxel', ...
    struct('n1', n1, 'n2', n2, 'nz', nz, 'mapName', mapName, 'ax', ax));

sigFul = squeeze(viewerData.ful(n1,n2,nz,:));
sigPar = squeeze(viewerData.par(n1,n2,nz,:));
sigDif = squeeze(viewerData.dif(n1,n2,nz,:));

sigFul = apodizeGW_v2(sigFul, viewerData.EPSIDur, viewerData.linewidthHz.ful, 0);
sigPar = apodizeGW_v2(sigPar, viewerData.EPSIDur, viewerData.linewidthHz.par, 0);
sigDif = apodizeGW_v2(sigDif, viewerData.EPSIDur, viewerData.linewidthHz.dif, 0);

viewerState = getappdata(figMap, 'SpectrumViewerState');
if ~isValidViewerState(viewerState)
    viewerState = createSpectrumFigure(figMap, viewerData.X_ppm, sigFul, sigPar, sigDif, mapName, n1, n2, nz);
    setappdata(figMap, 'SpectrumViewerState', viewerState);
else
    updateSpectrumFigure(viewerState, sigFul, sigPar, sigDif, mapName, n1, n2, nz);
end

drawnow limitrate
end

function viewerData = loadViewerDataFromBase()
viewerData = struct();
viewerData.X_ppm = evalin('base', 'X_ppm');
viewerData.EPSIDur = evalin('base', 'EPSIDur');
viewerData.linewidthHz = evalin('base', 'outputLinewidthHz');
viewerData.ful = evalin('base', 'xfData_ful_Filt');
viewerData.par = evalin('base', 'xfData_par_Filt');
viewerData.dif = evalin('base', 'xfData_dif_Filt');
end

function updateSelectionMarker(figMap, ax, n1, n2)
marker = getappdata(figMap, 'SpectrumSelectionMarker');
if isValidMarker(marker) && isequal(marker.ax, ax)
    set(marker.circle, 'XData', n2, 'YData', n1);
    set(marker.cross, 'XData', n2, 'YData', n1);
    return
end

delete(findobj(figMap, 'Tag', 'SelectedSpectrumVoxel'));
hold(ax, 'on')
circle = plot(ax, n2, n1, 'r+', 'MarkerSize', 10, 'LineWidth', 1.4, ...
    'Tag', 'SelectedSpectrumVoxel');
cross = plot(ax, n2, n1, 'w.', 'MarkerSize', 8, 'LineWidth', 1.2, ...
    'Tag', 'SelectedSpectrumVoxel');
set([circle cross], 'HitTest', 'off', 'PickableParts', 'none');
hold(ax, 'off')
setappdata(figMap, 'SpectrumSelectionMarker', struct('ax', ax, 'circle', circle, 'cross', cross));
end

function tf = isValidMarker(marker)
tf = isstruct(marker) && ...
    isScalarGraphics(marker.ax) && ...
    isScalarGraphics(marker.circle) && ...
    isScalarGraphics(marker.cross);
end

function state = createSpectrumFigure(figMap, xPpm, sigFul, sigPar, sigDif, mapName, n1, n2, nz)
figSpec = figure(6);
clf(figSpec);
set(figSpec, 'Name', 'SLOW-MRSI voxel spectra', 'NumberTitle', 'off', ...
    'Tag', 'SlowMRSISpectrumViewer', 'KeyPressFcn', @moveSelectedSpectrumVoxel, ...
    'Color', 'w');

axFull = subplot(1,2,1, 'Parent', figSpec);
Xlim = [min(xPpm), max(xPpm)];
fulLine = plot(axFull, xPpm, 10000*real(sigFul), 'LineWidth', 1.2); hold(axFull, 'on')
parLine = plot(axFull, xPpm, 10000*real(sigPar), 'LineWidth', 1.2); hold(axFull, 'off')
xlim(Xlim)
grid(axFull, 'on')
set(axFull, 'XDir', 'reverse')
legend(axFull, {'SLOW-ful', 'SLOW-par'}, 'Location', 'northeast')
title(axFull, sprintf('%s map selection', mapName), 'Interpreter', 'none')

axDiff = subplot(1,2,2, 'Parent', figSpec);
difLine = plot(axDiff, xPpm, 10000*real(sigDif), 'm', 'LineWidth', 1.2);
set(axDiff, 'XDir', 'reverse')
xlim(Xlim)
grid(axDiff, 'on')
xlabel(axDiff, 'ppm')
ylabel(axDiff, 'a.u.')
title(axDiff, sprintf('voxel: nx=%d, ny=%d, nz=%d', n1, n2, nz), 'Interpreter', 'none')
legend(axDiff, {'SLOW-dif'}, 'Location', 'northwest')

spectrumFigureSize = getappdata(figMap, 'SpectrumFigureSize');
if isempty(spectrumFigureSize)
    spectrumFigureSize = [12, 5];
end
setFigureSize(figSpec, spectrumFigureSize);

state = struct( ...
    'fig', figSpec, ...
    'axFull', axFull, ...
    'axDiff', axDiff, ...
    'fulLine', fulLine, ...
    'parLine', parLine, ...
    'difLine', difLine);
localTightYLimits(axFull, [10000*real(sigFul(:)); 10000*real(sigPar(:))]);
localTightYLimits(axDiff, 10000*real(sigDif(:)));
end

function updateSpectrumFigure(state, sigFul, sigPar, sigDif, mapName, n1, n2, nz)
set(state.fulLine, 'YData', 10000*real(sigFul));
set(state.parLine, 'YData', 10000*real(sigPar));
set(state.difLine, 'YData', 10000*real(sigDif));
title(state.axFull, sprintf('%s map selection', mapName), 'Interpreter', 'none')
title(state.axDiff, sprintf('voxel: nx=%d, ny=%d, nz=%d', n1, n2, nz), 'Interpreter', 'none')
localTightYLimits(state.axFull, [10000*real(sigFul(:)); 10000*real(sigPar(:))]);
localTightYLimits(state.axDiff, 10000*real(sigDif(:)));
end

function state = localEmptyViewerState()
state = struct( ...
    'fig', gobjects(0), ...
    'axFull', gobjects(0), ...
    'axDiff', gobjects(0), ...
    'fulLine', gobjects(0), ...
    'parLine', gobjects(0), ...
    'difLine', gobjects(0));
end

function tf = isValidViewerState(state)
tf = isstruct(state) && ...
    isScalarGraphics(state.fig) && ...
    isScalarGraphics(state.axFull) && ...
    isScalarGraphics(state.axDiff) && ...
    isScalarGraphics(state.fulLine) && ...
    isScalarGraphics(state.parLine) && ...
    isScalarGraphics(state.difLine);
end

function tf = isScalarGraphics(graphicsHandle)
tf = isscalar(graphicsHandle) && isgraphics(graphicsHandle);
end

function localTightYLimits(ax, values)
values = values(isfinite(values));
if isempty(values)
    return
end

yMin = min(values);
yMax = max(values);
if yMin == yMax
    margin = max(1, abs(yMin)*0.05);
else
    margin = 0.05*(yMax - yMin);
end
ylim(ax, [yMin - margin, yMax + margin])
end
