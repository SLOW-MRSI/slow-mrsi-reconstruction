function [Y] = linspaceGW(a,b,N)
%LINSPACEGW Linspace variant used by the SLOW-MRSI workflow.
% SLOW-MRSI package maintenance: Dr. Guodong Weng, University of Bern, 2026-07-04

Y = linspace(a, b - (b-a)/N,  N);

end
