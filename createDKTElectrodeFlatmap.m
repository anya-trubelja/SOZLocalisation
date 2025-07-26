function createDKTElectrodeFlatmap(subject_id, hemi, regionNames, electrodeXLSX)
%%========================================================================
% Function:   createDKTElectrodeFlatmap
% Purpose:    Generate a flatmap patch for user-specified DKT regions,
%             integrate electrode-contact vertices, grow geodesic buffers,
%             visualize, flatten, and open both original and flattened patch
%             in Freeview.
%
% Author:     Anya Trubelja
% Created:    2025-04-29
% Last Edit:  2025-07-26
%
% Requirements:
%   • MATLAB
%   • Freesurfer installation (FREESURFER_HOME set)
%   • computeElectrodeContactTable.m on MATLAB path
%   • freesurferMatlabLibrary-master on MATLAB path
%
% Usage:
%   createDKTElectrodeFlatmap('sub20','rh', {...}, '/path/to/MRIcentered.xlsx')
%
% Inputs:
%   subject_id    - Freesurfer subject folder name (string)
%   hemi          - Hemisphere: 'lh' or 'rh'
%   regionNames   - Cell array of DKT region name strings
%   electrodeXLSX - Path to MRI-centered electrode locations spreadsheet
%
% Workflow:
%   0) Initialize Freesurfer env and MATLAB paths
%   1) Define subject, surface, and label directories
%   2) Load and validate pial and inflated surfaces
%   3) Mask ROI vertices from DKT labels
%   4) Read and project electrodes via helper function
%   5) Grow a 10 mm geodesic buffer around contacts
%   6) Compute regionFaces and border edges
%   7) Visualize inflated surface + ROI + buffers + electrodes
%   8) Visualize patch + electrodes inside ROI
%   9) Write DKT patch file (.patch)
%  10) Flatten patch and open in Freeview
%%========================================================================

%% 0. FreeSurfer environment
fsHome = getenv('FREESURFER_HOME');
if isempty(fsHome)
  fsHome = '/Applications/freesurfer/';
  setenv('FREESURFER_HOME', fsHome);
end
[~, sdir] = system(sprintf('bash -lc "source %s/SetUpFreeSurfer.sh; echo $SUBJECTS_DIR"', fsHome));
sdir = strtrim(sdir);
if ~isempty(sdir)
  setenv('SUBJECTS_DIR', sdir);
end
subjects_dir = getenv('SUBJECTS_DIR');
if isempty(subjects_dir)
  subjects_dir = '/Applications/freesurfer/subjects';
  setenv('SUBJECTS_DIR', subjects_dir);
  warning('Using default SUBJECTS_DIR: %s', subjects_dir);
end
fstoolbox = 'Toolboxes/freesurferMatlabLibrary-master';
addpath(genpath(fstoolbox)); 

%% 1. Define directories
subj_dir  = fullfile(subjects_dir, subject_id);
surf_dir  = fullfile(subj_dir, 'surf');
label_dir = fullfile(subj_dir, 'label', 'roi_labels');

%% 2. Read + validate pial + inflated surfaces
pial_surf_file     = fullfile(surf_dir, sprintf('%s.pial',    hemi));
inflated_surf_file = fullfile(surf_dir, sprintf('%s.inflated',hemi));

[vertex_pial_raw, faces_pial_raw] = freesurfer_read_surf(pial_surf_file);
if min(faces_pial_raw(:)) == 0
    faces_pial = faces_pial_raw + 1;
else
    faces_pial = faces_pial_raw;
end
if max(faces_pial(:)) > size(vertex_pial_raw,1) || min(faces_pial(:)) < 1
    error('faces_pial indices outside [1 .. %d].', size(vertex_pial_raw,1));
end

[vertex_inflated_raw, faces_inflated_raw] = freesurfer_read_surf(inflated_surf_file);
if min(faces_inflated_raw(:)) == 0
    faces_inflated = faces_inflated_raw + 1;
else
    faces_inflated = faces_inflated_raw;
end
if max(faces_inflated(:)) > size(vertex_inflated_raw,1) || min(faces_inflated(:)) < 1
    error('faces_inflated indices outside [1 .. %d].', size(vertex_inflated_raw,1));
end

vertex_pial     = vertex_pial_raw;
vertex_inflated = vertex_inflated_raw;

fprintf('… Loaded and validated pial+inflated surfaces.\n');

%% 3. Mask ROI vertices from DKT
nV_pial = size(vertex_pial,1);
validV  = false(nV_pial, 1);
for i = 1:numel(regionNames)
  lblFile = fullfile(label_dir, sprintf('%s.%s.label', hemi, regionNames{i}));
  if ~exist(lblFile, 'file')
    error('Missing label: %s', lblFile);
  end
  fid = fopen(lblFile, 'r'); fgetl(fid);
  C    = textscan(fid, '%d %*f %*f %*f'); fclose(fid);
  verts = C{1} + 1;
  validV(verts) = true;
end
validV_ROI = validV;  % copy of the pure ROI mask
fprintf('… Built ROI mask from DKT labels.\n');

%% 4. Read + project electrodes via helper
T = computeElectrodeContactTable(subj_dir, hemi, electrodeXLSX);
if isempty(T)
    warning('No grey-matter contacts (isgrey==1). Exiting.');
    return;
end
elecV = T.PialVertexIdx;
validV(elecV) = true;

% Find electrodes truly inside the DKT ROI
roiVerts   = find(validV_ROI);
isInROI    = ismember(elecV, roiVerts);
elecInROI  = elecV(isInROI);

fprintf('… Computed electrode_contact_matrix (%d × 9).\n', height(T));

%% 5. Grow a 10 mm geodesic buffer around electrodes (on inflated mesh)
ed1_all = sort(faces_inflated(:,[1 2]), 2);
ed2_all = sort(faces_inflated(:,[2 3]), 2);
ed3_all = sort(faces_inflated(:,[3 1]), 2);
allEdges = [ed1_all; ed2_all; ed3_all];
[uniqueE_all, ~, ic_all] = unique(allEdges,'rows');
us = uniqueE_all(:,1);
vs = uniqueE_all(:,2);
w  = sqrt(sum((vertex_inflated(us,:) - vertex_inflated(vs,:)).^2, 2));
G = graph(us, vs, w);

uElecV   = unique(elecV);
Dall     = distances(G, uElecV);
dmin_all = min(Dall, [], 2);
bufV     = find(dmin_all <= 10);
validV(bufV) = true;

fprintf('… Grown 10 mm geodesic buffer around electrodes (%d vertices in buffer).\n', numel(bufV));

%% 6. Compute regionFaces & border edges using only validV_ROI 
patch_idx_ROI = find(validV_ROI);
maskFaces     = all(ismember(faces_inflated, patch_idx_ROI), 2);
regionFaces   = faces_inflated(maskFaces, :);

ed1_pf      = sort(regionFaces(:,[1 2]), 2);
ed2_pf      = sort(regionFaces(:,[2 3]), 2);
ed3_pf      = sort(regionFaces(:,[3 1]), 2);
edgesAll_pf = [ed1_pf; ed2_pf; ed3_pf];
[uniqueE_pf, ~, ic_pf] = unique(edgesAll_pf, 'rows');
counts_pf   = accumarray(ic_pf, 1);
borderE     = uniqueE_pf(counts_pf == 1, :);

fprintf('… Computed regionFaces (%d faces) and border edges (%d edges).\n', size(regionFaces,1), size(borderE,1));

%% 7. Figure 1: Inflated + ROI Patch + Buffers + Borders + Electrodes
figure('Name','Inflated Surface + ROI Patch + Buffers + Electrodes','Color','w');
hold on;

% Entire inflated mesh
patch( ...
    'Vertices',     vertex_inflated, ...
    'Faces',        faces_inflated, ...
    'FaceColor',    [0.8 0.8 0.8], ...
    'EdgeColor',    'none', ...
    'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.4, ...
    'DiffuseStrength', 0.6 ...
);

% ROI patch (semi-transparent blue)
patch( ...
    'Vertices',     vertex_inflated, ...
    'Faces',        regionFaces, ...
    'FaceColor',    [0.3 0.6 0.9], ...
    'EdgeColor',    'none', ...
    'FaceAlpha',    0.3, ...
    'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.2, ...
    'DiffuseStrength', 0.7 ...
);

% Border edges (dark blue lines)
for i = 1:size(borderE,1)
    v1 = borderE(i,1); v2 = borderE(i,2);
    xline = [vertex_inflated(v1,1), vertex_inflated(v2,1)];
    yline = [vertex_inflated(v1,2), vertex_inflated(v2,2)];
    zline = [vertex_inflated(v1,3), vertex_inflated(v2,3)];
    plot3(xline, yline, zline, 'Color', [0 0 0.5], 'LineWidth', 1.0);
end

% Buffer (green dots)
scatter3( ...
    vertex_inflated(bufV,1), ...
    vertex_inflated(bufV,2), ...
    vertex_inflated(bufV,3), ...
    36, 'g', 'filled' ...
);

% All electrodes (red dots)
scatter3( ...
    vertex_inflated(elecV,1), ...
    vertex_inflated(elecV,2), ...
    vertex_inflated(elecV,3), ...
    80, 'r', 'filled' ...
);

axis equal off;
camlight headlight;
view(3);
title(sprintf('Subject %s, %s Inflated + ROI Patch + Buffers + Electrodes', subject_id, hemi), 'FontSize', 14);
hold off;

%% 8. Figure 2: Patch + Electrodes Inside ROI Only
figure('Name','Patch + Electrodes (in ROI)','Color','w');
hold on;

patch( ...
    'Vertices',     vertex_inflated, ...
    'Faces',        regionFaces, ...
    'FaceColor',    [0.3 0.6 0.9], ...
    'EdgeColor',    'none', ...
    'FaceAlpha',    0.5, ...
    'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.3, ...
    'DiffuseStrength', 0.7 ...
);

scatter3( ...
    vertex_inflated(elecInROI,1), ...
    vertex_inflated(elecInROI,2), ...
    vertex_inflated(elecInROI,3), ...
    100, 'r', 'filled' ...
);

axis equal off;
camlight headlight;
view(3);
title(sprintf('Subject %s: %s ROI Patch + Electrodes (in ROI)', subject_id, hemi), 'FontSize', 14);
hold off;

%% 9. Write patch to file (using patch_idx_ROI and inflated coordinates)
patch_idx = patch_idx_ROI;           % the vertices of the DKT‐only patch
V3d = vertex_inflated;               % use inflated coordinates for V3d

patch_file = fullfile(surf_dir, sprintf('%s.dkt.electrode.patch', hemi));
flat_file  = fullfile(surf_dir, sprintf('%s.dkt.electrode.flat.patch', hemi));
fid = fopen(patch_file,'wb','b');
fwrite(fid, int32(-1), 'int32');
fwrite(fid, int32(numel(patch_idx)), 'int32');

for ii = 1:numel(patch_idx)
  vid = patch_idx(ii);
  isB = any(borderE(:) == vid);             % whether this vertex is on the border
  raw = int32(vid) * (isB * -2 + 1);         % if border, multiply by -1; else keep positive
  fwrite(fid, raw, 'int32');
  fwrite(fid, single(V3d(vid,:)), 'float32');
end
fclose(fid);

%% 10. Flatten and view in Freeview
system(sprintf('"%s/bin/mris_flatten" %s %s', fsHome, patch_file, flat_file));
infl = fullfile(surf_dir, sprintf('%s.inflated', hemi));
system(sprintf('"%s/bin/freeview" -f %s:patch=%s &', fsHome, infl, patch_file));
system(sprintf('"%s/bin/freeview" -f %s:patch=%s &', fsHome, infl, flat_file));

end

%%========================================================================
%% Helper: plot_mesh_simple =============================================
function plot_mesh_simple(vertex, faces)
% Quick patch renderer for zero-based faces
%   vertex: [N×3] coords; faces: [M×3] zero-based

    % If transposed, fix orientation
    if size(faces,1)==3 && size(faces,2)~=3
        faces = faces';
    end
    faces = faces + 1;  % to 1-based

    patch('Vertices',vertex,'Faces',faces, ...
          'FaceColor',[0.8 0.8 0.8],'EdgeColor','none', ...
          'FaceLighting','gouraud','AmbientStrength',0.4,'DiffuseStrength',0.6);
    axis equal off; camlight headlight; material dull;
end
%%========================================================================
