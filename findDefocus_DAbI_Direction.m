function [z] = findDefocus_DAbI_Direction(Nk_up, Nk_down, CTF_vir, mycoord, imStack, overlap_mask, angle_rad, z,row1, row2, col1,col2)
% findDefocus_DAbI_Direction determines the sign of the defocus value (z)
% by comparing the fringe density by adding known defocus aberration on 
% both side. The direction where the fringes are denser should be the
% direction of defocus.

% This method is only used for large defocus distance.

% Inputs:
% - Nk_up / Nk_down: Estimated fringe density for z +/- ∆z cases.
% - CTF_vir:         Virtual, known defocus aberration.
% - mycoord:         Illumination wavevectors (pixel shift).
% - imStack:         Two raw images (2-LED acquisition).
% - overlap_mask:    Binary mask for overlapping fringe region S.
% - angle_rad:       Rotation angle to align fringes horizontally.
% - z:               Estimated magnitude of defocus.
% - row1, row2, col1, col2: Indices to crop the aligned region of interest.

% Output:
% - z: Final signed defocus (positive or negative), based on fringe asymmetry.

% To be released


end