%%========================================================================
% Script:     makeSOZNifti.m
% Purpose:    Load a surface‐based heatmap (.fig), convert vertex intensities 
%             to voxel space, write unsmoothed and smoothed NIfTI volumes, 
%             and display central orthogonal slices.
%
% Author:     Anya Trubelja
% Created:    2025‑06‑04
% Last Edit:  2025‑07‑26
%
% Inputs:
%   figFile      - Path to saved heatmap .fig file
%   t1File       - Path to template anatomical NIfTI (.nii)
%   x_offset_mm  - Manual R/L translation (mm)
%   y_offset_mm  - Manual A/P translation (mm)
%   z_offset_mm  - Manual S/I translation (mm)
%
% Outputs:
%   heatmap_nosmooth.nii   - Float32 NIfTI of raw projected heatmap
%   heatmap_smooth.nii     - Float32 NIfTI after 3D Gaussian smoothing
%   Figures showing sagittal, coronal, axial slices of the smoothed map
%
% Workflow:
%   1. Extract vertices & color data from .fig
%   2. Map to voxel indices using SPM’s affine
%   3. Write and verify unsmoothed volume
%   4. Smooth in MATLAB; write and verify smoothed volume
%   5. Display central slices
%
% Requirements:
%   • MATLAB with SPM12 on path
%   • Image Processing Toolbox for imgaussfilt3
%
%%========================================================================

%% USER SETTINGS — adjust these three parameters to your data
figFile    = 'SOZHeatmap020.fig';              % Path to your saved heatmap .fig
t1File     = '/Users/neurallabpostgrads/Desktop/Resection_masks/Subject_020/preop_aligned2_020.nii';  % Path to template anatomical NIfTI
% The offsets below let you manually translate the overlay in mm:
x_offset_mm = 0;   % Right (+) / Left (–) translation
y_offset_mm = 0;   % Anterior (+) / Posterior (–) translation
z_offset_mm = 0;   % Superior (+) / Inferior (–) translation

%% 0. Extract surface geometry and data from heatmap .fig
% Open the .fig invisibly and locate the patch object containing the surface
fig = openfig(figFile, 'invisible');
h   = findobj(fig, 'Type', 'patch');

% Retrieve vertex coordinates and associated color data
vertex_pial    = h.Vertices;           % [nVerts × 3] matrix of (X,Y,Z) positions
vertex_colors  = h.FaceVertexCData;    % [nVerts × 3] RGB values per vertex
cmap           = colormap(fig);        % Figure’s colormap
nC             = size(cmap,1);         % Number of colors

% Convert each vertex color to a normalized density value [0,1]
vertex_density = zeros(size(vertex_colors,1), 1);
for ii = 1:size(vertex_colors,1)
    % Find nearest colormap entry (minimize squared RGB distance)
    [~, idx] = min( sum( (cmap - vertex_colors(ii,:)).^2, 2 ) );
    % Scale index to [0,1]
    vertex_density(ii) = (idx - 1) / (nC - 1);
end

% Clean up
close(fig);

%% 1. Configure SPM environment
% Add SPM’s toolbox folder to MATLAB path and initialize defaults
spmPath = 'Toolboxes/spm12-master';  
addpath(spmPath);
spm('defaults', 'fmri');
spm_jobman('initcfg');

%% 2. Load anatomical template and enforce float32 output
V      = spm_vol(t1File);    % Header info for template NIfTI
V.dt   = [16 0];             % Force data type to float32
dim    = V.dim;              % Volume dimensions [X Y Z]
affine = V.mat;              % 4×4 RAS→voxel affine transform
outVol = zeros(dim);         % Preallocate output volume

%% 2.5. Apply manual translation to surface in RAS space
fprintf('Applying manual offset (R/L: %+d mm, A/P: %+d mm, S/I: %+d mm)\n', ...
        x_offset_mm, y_offset_mm, z_offset_mm);

% Translate RAS coordinates directly
vertex_pial(:,1) = vertex_pial(:,1) + x_offset_mm;
vertex_pial(:,2) = vertex_pial(:,2) + y_offset_mm;
vertex_pial(:,3) = vertex_pial(:,3) + z_offset_mm;

%% 3. Project surface vertices into voxel space
nVert   = size(vertex_pial, 1);
% Convert to homogeneous RAS, then to voxel indices
RAS_hom = [vertex_pial'; ones(1, nVert)];
VOX_hom = affine \ RAS_hom;

% Round to nearest voxel center
i = round(VOX_hom(1, :));
j = round(VOX_hom(2, :));
k = round(VOX_hom(3, :));

% Identify which indices fall within the volume bounds
valid = i>=1 & i<=dim(1) & j>=1 & j<=dim(2) & k>=1 & k<=dim(3);
inds  = sub2ind(dim, i(valid), j(valid), k(valid));

% Map densities into the output volume
outVol(inds) = vertex_density(valid);

fprintf('Mapped nonzero voxels: %d of %d (%.2f%%)\n', ...
        nnz(outVol), numel(outVol), 100*nnz(outVol)/numel(outVol) );

%% 4. Write unsmoothed NIfTI and verify integrity
V.fname = 'heatmap_nosmooth.nii';     % Output filename
spm_write_vol(V, outVol);             % Write volume as float32

% Reload and count nonzero voxels to confirm write
W  = spm_vol(V.fname);
Y2 = spm_read_vols(W);

fprintf('Verification (unsmoothed) → nonzeros: %d / %d (%.2f%%)\n\n', ...
        nnz(Y2), numel(Y2), 100*nnz(Y2)/numel(Y2) );

%% 5. Smooth heatmap in MATLAB and write smoothed NIfTI
% Compute voxel dimensions from affine
voxelSizes = sqrt(sum(affine(1:3,1:3).^2, 1));  
% Define Gaussian smoothing (FWHM) in mm
fwhm_mm   = [5 5 5];
sigma_mm  = fwhm_mm / (2*sqrt(2*log(2)));
sigma_vox = sigma_mm ./ voxelSizes;    % Convert to voxel units

% Apply 3D Gaussian filter (replicate padding at edges)
smoothedVol = imgaussfilt3(outVol, sigma_vox, 'Padding', 'replicate');

% Prepare header for smoothed volume
Vs       = V;
Vs.fname = 'heatmap_smooth.nii';
Vs.dt    = [16 0];             % Ensure float32 again
spm_write_vol(Vs, smoothedVol);

% Reload and verify
U  = spm_vol(Vs.fname);
Z2 = spm_read_vols(U);

fprintf('Verification (smoothed) → nonzeros: %d / %d (%.2f%%)\n\n', ...
        nnz(Z2), numel(Z2), 100*nnz(Z2)/numel(Z2) );

%% 6. Visualize central orthogonal slices of smoothed heatmap
midX = round(dim(1)/2);
midY = round(dim(2)/2);
midZ = round(dim(3)/2);

figure;
% Sagittal view (X)
subplot(1,3,1);
imagesc(squeeze(Z2(midX,:,:))'); 
axis equal tight off;
title('Sagittal');   
colormap(jet); 
colorbar;

% Coronal view (Y)
subplot(1,3,2);
imagesc(squeeze(Z2(:,midY,:))'); 
axis equal tight off;
title('Coronal');    
colormap(jet); 
colorbar;

% Axial view (Z)
subplot(1,3,3);
imagesc(Z2(:,:,midZ));        
axis equal tight off;
title('Axial');      
colormap(jet); 
colorbar;
