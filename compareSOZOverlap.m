%%========================================================================
% Script:     compareSOZOverlap.m
% Purpose:    Compare a functional heatmap volume against a resection mask:
%             - Resample resection mask into heatmap space
%             - Threshold heatmap at 90th percentile to predict seizure onset
%             - Compute true/false overlaps and volumes
%             - Report quantitative metrics of prediction accuracy
%
% Author:     Anya Trubelja
% Created:    2025-07-02
% Last Edit:  2025-07-26
%
% Requirements:
%   • MATLAB
%   • SPM12 on MATLAB path (spm12-master)
%   • Image Processing Toolbox (for prctile)
%
% Inputs (configure at top):
%   heatmap_nii - Path to float32 heatmap NIfTI
%   mask_nii    - Path to binary resection mask NIfTI
%
% Outputs:
%   • Console report of voxel counts, volumes (mm^3), and percentages:
%       - Resection mask size
%       - Predicted-onset volume
%       - True overlap metrics
%       - False positive metrics
%
% Workflow:
%   0) Initialize SPM environment
%   1) Load heatmap and mask volumes
%   2) Resample mask into heatmap voxel space
%   3) Threshold heatmap at 90th percentile
%   4) Compute voxel volumes from affine
%   5) Calculate counts, volumes, and percentages
%   6) Print summary report
%%========================================================================

%% User inputs — adjust to your file locations
heatmap_nii = '/Users/neurallabpostgrads/Desktop/Resection_masks/Subject_020/auto_heatmap_smooth_020.nii';  % Smoothed heatmap volume
mask_nii    = '/Users/neurallabpostgrads/Desktop/Resection_masks/Subject_020/resection_mask_020.nii';        % Binary resection mask

%% 0. Prepare SPM environment for NIfTI I/O
spmPath = 'Toolboxes/spm12-master';       % Path to SPM12 folder
addpath(spmPath);                         % Add SPM to MATLAB path
spm('defaults','fmri');                   % Set fMRI defaults
spm_jobman('initcfg');                    % Initialize job manager

%% 1. Load heatmap and resection mask volumes
Vh = spm_vol(heatmap_nii);                % Header for heatmap
H  = spm_read_vols(Vh);                   % Heatmap data [X×Y×Z]
Vm = spm_vol(mask_nii);                   % Header for mask
M0 = spm_read_vols(Vm) > 0;               % Binary mask [X'×Y'×Z'] logical

%% 2. Resample mask into heatmap voxel space
% Build grid of heatmap voxel indices in homogeneous coords
[xh,yh,zh] = ndgrid(1:Vh.dim(1), 1:Vh.dim(2), 1:Vh.dim(3));
nVox       = numel(xh);
homH       = [xh(:)'; yh(:)'; zh(:)'; ones(1,nVox)];

% Transform heatmap voxels to world RAS coordinates
RAS_world  = Vh.mat * homH;                % 4×N matrix

% Transform world RAS to mask voxel coordinates
voxM       = Vm.mat \ RAS_world;        % 4×N matrix
xm = round(voxM(1,:));                    % Mask-space i indices
ym = round(voxM(2,:));                    % Mask-space j indices
zm = round(voxM(3,:));                    % Mask-space k indices

% Identify valid mappings within mask bounds
valid = xm>=1 & xm<=Vm.dim(1) & ...
        ym>=1 & ym<=Vm.dim(2) & ...
        zm>=1 & zm<=Vm.dim(3);

% Initialize resampled mask M in heatmap space
M = false(Vh.dim);
% Map valid heatmap voxels to mask voxels
idxH = sub2ind(Vh.dim, xh(valid), yh(valid), zh(valid));
idxM = sub2ind(Vm.dim, xm(valid), ym(valid), zm(valid));
M(idxH) = M0(idxM);

%% 3. Threshold heatmap to create predicted-onset binary mask
% Extract nonzero heatmap values
vals = H(H > 0);
% Compute 90th percentile threshold (top 10% of nonzero voxels)
threshold90pct = prctile(vals, 90);
% Generate predicted mask: voxels >= threshold
predictedMask = H >= threshold90pct;

%% 4. Compute voxel volume from affine matrix
% Voxel size in mm along each axis
voxelSizes_mm   = sqrt(sum(Vh.mat(1:3,1:3).^2, 1));  % [dx dy dz]
voxelVolume_mm3 = prod(voxelSizes_mm);               % mm³ per voxel

%% 5. Calculate counts, volumes, and percentages
nResVox     = nnz(M);                               % Resection voxels count
volRes_mm3  = nResVox * voxelVolume_mm3;             % Resection volume (mm³)

nPredVox    = nnz(predictedMask);                   % Predicted voxels count
volPred_mm3 = nPredVox * voxelVolume_mm3;            % Predicted volume (mm³)

nTrueVox    = nnz(M & predictedMask);               % True positives count
volTrue_mm3 = nTrueVox * voxelVolume_mm3;            % True positive volume (mm³)

nFalseVox   = nnz(~M & predictedMask);               % False positives count
volFalse_mm3= nFalseVox * voxelVolume_mm3;           % False positive volume (mm³)

pctTrue  = 100 * nTrueVox  / nResVox;               % % of resection covered
pctFalse = 100 * nFalseVox / nPredVox;              % % of preds outside resection

%% 6. Display results summary
fprintf('\n=== 90%%‐percentile cutoff = %.4g — Evaluation Metrics ===\n\n', threshold90pct);

fprintf('Resection mask:          %6d voxels (%.1f mm³)\n', nResVox, volRes_mm3);
fprintf('Predicted-onset mask:    %6d voxels (%.1f mm³)\n\n', nPredVox, volPred_mm3);

fprintf('True overlap (in mask & predicted):\n');
fprintf('  %6d voxels (%.1f mm³), %.1f%% of resection volume\n\n', nTrueVox, volTrue_mm3, pctTrue);

fprintf('False positives (predicted outside mask):\n');
fprintf('  %6d voxels (%.1f mm³), %.1f%% of predicted volume\n\n', nFalseVox, volFalse_mm3, pctFalse);

