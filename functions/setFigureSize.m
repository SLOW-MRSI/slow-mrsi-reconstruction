function setFigureSize(figHandle, figSizeInches)
%SETFIGURESIZE Set only figure width and height while preserving location.
%
% Author/Maintainer: Dr. Guodong Weng, University of Bern
% Date: 2026-07-04

figHandle.Units = 'inches';
figPosition = figHandle.Position;
figPosition(3:4) = figSizeInches;
figHandle.Position = figPosition;
end
