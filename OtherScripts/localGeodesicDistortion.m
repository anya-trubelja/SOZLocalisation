%%========================================================================
% Script:     localGeodesicDistortion.m
% Purpose:    Quantify and visualize local distance distortion between 3D
%             geodesic and 2D flatmap distances across the cortical surface.
%
% Author:     Anya Trubelja
% Created:    2024-12-12
% Last Edit:  2025-07-26
%
% Requirements:
%   • MATLAB
%   • Freesurfer installation (FREESURFER_HOME & SUBJECTS_DIR)
%   • Freesurfer MATLAB I/O (read_surf, read_patch)
%   • Graph toolbox for geodesic computations
%   • bfs_neighborhood helper to obtain k-ring vertices
%
% Inputs (set in user section):
%   subject          - Subject folder name (e.g. 'subject020')
%   hemi             - Hemisphere code: 'lh' or 'rh'
%   flat_patch_file  - Path to .flat.patch file for flatmap
%
% Outputs:
%   • local_distance_distortion: [%nVerts×1] signed % error between
%     2D flatmap and 3D geodesic distances in each vertex’s local 2-ring region
%   • Scatter heatmap of distortion on flatmap
%%========================================================================

%% User inputs — define your subject and hemisphere here
subject         = 'subject020';           % Subject folder name within SUBJECTS_DIR
hemi            = 'rh';                   % Hemisphere: 'lh' or 'rh'
flat_patch_file = '/Applications/freesurfer/subjects/subject020/surf/rh.dkt.electrode.flat.patch';

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

%% 1. Load inflated surface mesh
if ~exist(surf_path,'file')
    error('Inflated surface not found: %s', surf_path);
end
[vertices, faces] = read_surf(surf_path);   % Zero-based faces
% Convert to MATLAB 1-based indexing and transpose for consistency
faces_inflated  = (faces + 1)';  % 3×nFaces
vertex_inflated = vertices';     % 3×nVerts

%% 2. Load flatmap patch for vertex XY coordinates
flatmap = read_patch(flat_patch_file);
nVerts  = size(vertex_inflated, 2);
nFaces  = size(faces_inflated, 2);

%% 3. Optional: preview decimated mesh (commented)
%[decel_faces, decel_verts] = reducepatch(faces, vertices, 0.1);
%figure;
%patch('Faces', decel_faces, 'Vertices', decel_verts, ...
%      'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',1);
%axis equal off; camlight headlight; lighting gouraud;

%% 4. Build adjacency list of 1-ring neighbors for each vertex
adj_list = cell(nVerts, 1);
for fIdx = 1:nFaces
    f = faces_inflated(:, fIdx);
    % Add each pair of vertices in the face as neighbors
    adj_list{f(1)} = unique([adj_list{f(1)}, f(2), f(3)]);
    adj_list{f(2)} = unique([adj_list{f(2)}, f(3), f(1)]);
    adj_list{f(3)} = unique([adj_list{f(3)}, f(1), f(2)]);
end

%% 5. Prepare flatmap XY lookup for each vertex
flatmap_coords = nan(2, nVerts);
idxs = flatmap.vno + 1;  % Convert zero-based vertex numbers
flatmap_coords(1, idxs) = flatmap.x;
flatmap_coords(2, idxs) = flatmap.y;

% Initialize distortion array
local_distance_distortion = nan(nVerts, 1);

%% 6. Compute local distance distortion for each vertex
for v = 1:nVerts
    % a) Get 2-ring neighborhood via BFS
    neighborhood = bfs_neighborhood(v, adj_list, 2);
    if numel(neighborhood) < 10, continue; end  % Skip small regions

    % b) Build weighted graph for 3D geodesic distances
    edges3D = [];
    weights3D = [];
    for vi = neighborhood
        nbrs = intersect(adj_list{vi}, neighborhood);
        for vj = nbrs
            if vi < vj  % avoid duplicates
                edges3D(end+1, :) = [find(neighborhood==vi), find(neighborhood==vj)];
                weights3D(end+1) = norm(vertex_inflated(:,vi) - vertex_inflated(:,vj));
            end
        end
    end
    if isempty(edges3D), continue; end
    G = graph(edges3D(:,1), edges3D(:,2), weights3D, numel(neighborhood));

    % Compute geodesic distances from center
    center_idx = find(neighborhood==v);
    dists3D = distances(G, center_idx);

    % c) Compute 2D Euclidean distances on flatmap
    flat_local = flatmap_coords(:, neighborhood);
    if any(isnan(flat_local), 'all'), continue; end
    centerXY = flat_local(:, center_idx);
    diffs2D  = flat_local - centerXY;
    dists2D  = sqrt(sum(diffs2D.^2, 1));

    % d) Compute signed % error (ignore self-distance)
    mask = dists3D > eps;
    signed_err = (dists2D(mask) - dists3D(mask)) ./ dists3D(mask);
    local_distance_distortion(v) = 100 * mean(signed_err);
end

%% 7. Plot distortion heatmap on flatmap
idxs = flatmap.vno + 1;
dist = local_distance_distortion(idxs);
valid = ~isnan(dist);
x = flatmap.x(valid);
y = flatmap.y(valid);
d = dist(valid);

% Clip extremes at 99th percentile
thr = prctile(abs(d), 99);
d(d>thr) = thr; d(d<-thr) = -thr;

figure;
scatter(x, y, 12, d, 'filled');
colormap('jet'); caxis([-thr thr]);
colorbar('Label','Signed Distortion (%)');
title('Local Geodesic Distance Distortion');
xlabel('Flatmap X'); ylabel('Flatmap Y');
axis equal off;

