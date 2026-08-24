function updateTsEPSIProgress(progress, overallFraction, stageName, blockFraction, detailText)
%UPDATETSEPSIPROGRESS Update the SLOW-MRSI progress window and handle cancel.
% SLOW-MRSI reconstruction package
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

if nargin < 5
    detailText = '';
end

if ~isstruct(progress) || ~isfield(progress, 'enabled') || ~progress.enabled
    return
end

if isempty(progress.h) || ~isgraphics(progress.h)
    return
end

if getappdata(progress.h, 'canceling')
    error('MATLAB:TsEPSIUserCanceled', 'SLOW-MRSI reconstruction canceled by user.');
end

progressHandles = [progress.stageText, progress.detailText, progress.blockText, ...
    progress.overallText, progress.barPatch];
if any(~isgraphics(progressHandles))
    return
end

overallFraction = min(max(overallFraction, 0), 1);
blockFraction = min(max(blockFraction, 0), 1);

set(progress.stageText, 'String', stageName);
set(progress.detailText, 'String', detailText);
set(progress.blockText, 'String', sprintf('Block: %5.1f%%', 100*blockFraction));
set(progress.overallText, 'String', sprintf('Overall: %5.1f%%', 100*overallFraction));
set(progress.barPatch, 'XData', [0 overallFraction overallFraction 0], 'YData', [0 0 1 1]);

if overallFraction >= 1
    set(progress.h, 'CloseRequestFcn', 'delete(gcbf)');
    if isfield(progress, 'actionButton') && isgraphics(progress.actionButton)
        set(progress.actionButton, 'String', 'Close', 'Callback', 'delete(gcbf)');
    end
end

drawnow
end
