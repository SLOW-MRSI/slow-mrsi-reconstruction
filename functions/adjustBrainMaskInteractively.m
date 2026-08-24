function [brainMask, params, mask_main, mask_wat] = adjustBrainMaskInteractively( ...
    waterMapSmooth, initialThresholdFraction, initialSkullMarginVoxels, previewMode)
%ADJUSTBRAINMASKINTERACTIVELY Preview and tune water-mask skull stripping.

params = struct();
params.thresholdFraction = initialThresholdFraction;
params.skullMarginVoxels = initialSkullMarginVoxels;
params.minObjectVoxels = 100;
params.closeRadius = 2;

while true
    [brainMask, mask_main, mask_wat, params.thresholdValue, params.skullMarginVoxels] = buildBrainMaskFromWaterMap( ...
        waterMapSmooth, params.thresholdFraction, params.skullMarginVoxels, ...
        params.minObjectVoxels, params.closeRadius);

    plotBrainMaskPreview(brainMask, previewMode, waterMapSmooth, params);

    if ~strcmpi(previewMode, 'on')
        return
    end

    [action, params] = waitForBrainMaskFigureAction(params);
    if strcmp(action, 'accept')
        [brainMask, mask_main, mask_wat, params.thresholdValue, params.skullMarginVoxels] = buildBrainMaskFromWaterMap( ...
            waterMapSmooth, params.thresholdFraction, params.skullMarginVoxels, ...
            params.minObjectVoxels, params.closeRadius);
        plotBrainMaskPreview(brainMask, previewMode, waterMapSmooth, params);
        return
    end
    if strcmp(action, 'cancel')
        return
    end
end
end

function [action, params] = waitForBrainMaskFigureAction(params)
fig = figure(4);
set(fig, 'Name', 'Brain mask adjustment', 'NumberTitle', 'off', ...
    'ToolBar', 'figure', 'MenuBar', 'figure');

controlPanel = uipanel(fig, 'Units', 'normalized', 'Position', [0.02 0.012 0.96 0.075], ...
    'BorderType', 'none');
uiFontSize = 9;
buttonFontSize = 9;

uicontrol(controlPanel, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.00 0.22 0.15 0.44], 'String', 'Threshold fraction', ...
    'HorizontalAlignment', 'left', 'FontSize', uiFontSize);
thresholdEdit = uicontrol(controlPanel, 'Style', 'edit', 'Units', 'normalized', ...
    'Position', [0.16 0.20 0.10 0.50], ...
    'String', num2str(params.thresholdFraction), 'FontSize', uiFontSize);

uicontrol(controlPanel, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.31 0.22 0.11 0.44], 'String', 'Skull margin', ...
    'HorizontalAlignment', 'left', 'FontSize', uiFontSize);
skullMarginEdit = uicontrol(controlPanel, 'Style', 'edit', 'Units', 'normalized', ...
    'Position', [0.43 0.20 0.08 0.50], ...
    'String', num2str(params.skullMarginVoxels), 'FontSize', uiFontSize);

uicontrol(controlPanel, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.62 0.16 0.10 0.58], 'String', 'Update', ...
    'Callback', @updateMask, 'FontSize', buttonFontSize);
uicontrol(controlPanel, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.75 0.16 0.10 0.58], 'String', 'Accept', ...
    'Callback', @acceptMask, 'FontSize', buttonFontSize);
uicontrol(controlPanel, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.88 0.16 0.10 0.58], 'String', 'Cancel', ...
    'Callback', @cancelMask, 'FontSize', buttonFontSize);

setappdata(fig, 'BrainMaskAction', 'cancel');
setappdata(fig, 'BrainMaskParams', params);
drawnow
uiwait(fig)

if isgraphics(fig)
    action = getappdata(fig, 'BrainMaskAction');
    params = getappdata(fig, 'BrainMaskParams');
else
    action = 'cancel';
end

    function updateMask(~, ~)
        updatedParams = readParamsFromControls(params);
        setappdata(fig, 'BrainMaskParams', updatedParams);
        setappdata(fig, 'BrainMaskAction', 'update');
        uiresume(fig)
    end

    function acceptMask(~, ~)
        updatedParams = readParamsFromControls(params);
        setappdata(fig, 'BrainMaskParams', updatedParams);
        setappdata(fig, 'BrainMaskAction', 'accept');
        if isgraphics(controlPanel)
            set(controlPanel, 'Visible', 'off')
        end
        drawnow
        uiresume(fig)
    end

    function cancelMask(~, ~)
        setappdata(fig, 'BrainMaskAction', 'cancel');
        uiresume(fig)
    end

    function updatedParams = readParamsFromControls(currentParams)
        updatedParams = currentParams;
        thresholdFraction = str2double(get(thresholdEdit, 'String'));
        skullMarginVoxels = str2double(get(skullMarginEdit, 'String'));
        if isfinite(thresholdFraction)
            updatedParams.thresholdFraction = max(0, thresholdFraction);
        end
        if isfinite(skullMarginVoxels)
            updatedParams.skullMarginVoxels = max(0, round(skullMarginVoxels));
        end
    end
end
