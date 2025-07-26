%%========================================================================
% Script:     ATMultilateration6.m
% Purpose:    Perform multilateration-based seizure onset zone (SOZ) localization
%             using electrode onset times, generate heatmap on pial surface,
%             and visualize results on flatmap and 3D pial mesh.
%
% Author:     Anya Trubelja
% Created:    2025-06-10
% Last Edit:  2025-07-26
%
% Adapted from: Thomas Hartigan's THMultilateration5.m
% Available at: https://tcdud.sharepoint.com/:u:/r/sites/TCD365-NeuralEngineering/Shared%20Documents/General/Neurosurgical%20Planning/Thomas%20Hartigan/Matlab%20Code/Thomas%27s%20matlab%20code/THMultilateration5.m?csf=1&web=1&e=L7dRF1
% 
% Requirements:
%   • MATLAB
%   • Freesurfer MATLAB tools (freesurferMatlabLibrary-master)
%   • Wavelet toolbox (toolbox_wavelets)
%   • Graph toolbox (toolbox_graph)
%   • computeElectrodeContactTable.m
%   • plot_mesh.m (helper for 3D surface plotting)
%
% Inputs (modify paths as needed):
%   file_path         - Path to Excel with electrode VertexIdx and OnsetTime
%   pial_file         - Freesurfer surface file labeled rh.pial.T1 or lh.pial.T1
%   flat_patch_file   - Path to .flat.patch file for electrode flatmap
%
% Workflow:
%   0) Initialize FreeSurfer and add toolboxes
%   1) Import electrode data and pial surface
%   2) Map electrode contacts to flatmap coordinates
%   3) Parse onset times and compute time differences
%   4) Compute pairwise propagation distances
%   5) Trilateration across all triplet combinations, buffer by 5 mm
%   6) Plot electrodes and origin points on flatmap
%   7) Project origin points back to pial mesh vertices
%   8) Compute Gaussian-blurred vertex density map
%   9) Plot heatmap on pial surface with lighting and labels
%  10) Compute 2D/3D clustering coherency metrics
%%========================================================================

%% User inputs — adjust to your file locations
file_path = '/Users/neurallabpostgrads/Desktop/Electrode_locations/Electrodes_020/OnsetTimes.xlsx';
pial_file = '/Applications/freesurfer/subjects/anya_thesis_subjects/subject020/surf/rh.pial.T1';
flat_patch_file = '/Applications/freesurfer/subjects/anya_thesis_subjects/subject020/surf/rh.dkt.electrode.flat.patch';

%% 0. FreeSurfer environment setup
% Ensure FREESURFER_HOME and SUBJECTS_DIR are set, then add MATLAB toolboxes
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
% Add required toolboxes for surface I/O, wavelets, and graph operations
addpath(genpath('Toolboxes/freesurferMatlabLibrary-master'));
addpath(genpath('Toolboxes/toolbox_wavelets'));
addpath(genpath('Toolboxes/toolbox_graph'));

%% 1. Import electrode contact times and pial mesh
data      = readtable(file_path, 'FileType', 'spreadsheet');
% Load pial vertex and face arrays from Freesurfer
[vertex_pial, faces_pial] = freesurfer_read_surf(pial_file);
% Load flatmap patch with electrode connections
flatmap = read_patch(flat_patch_file);

%% 2).Extract and map electrode contact indices to flatmap XY
% Excel provides VertexIdx referencing full pial mesh
VertexIdx = data.VertexIdx;
% Flatmap .ind is zero-based; convert to MATLAB 1-based
patchVtx = flatmap.ind + 1;
% Find matching flatmap indices for each contact
[found, loc] = ismember(VertexIdx, patchVtx);
if any(~found)
    error('Some VertexIdx values not in flat.patch; verify hemisphere and file.');
end
% Build an N×2 matrix of electrode XY coordinates on flatmap
electrode_coords = [flatmap.x(loc), flatmap.y(loc)];

%% 3. Parse raw onset times (HH:MM:SS:cc) and compute relative times
rawStr = string(data.OnsetTime);
parts  = split(rawStr, ':');
assert(size(parts,2)==4, 'OnsetTime must have 4 fields HH:MM:SS:cc');
H  = double(parts(:,1)); Mi = double(parts(:,2));
S  = double(parts(:,3)); CC = double(parts(:,4));
% Absolute times in seconds
t_sec = H*3600 + Mi*60 + S + CC/100;
% Relative to reference electrode (first contact)
time_points = t_sec - t_sec(1);
assert(numel(time_points) == numel(VertexIdx), 'Mismatch in counts');

disp('Relative onset times (s):'); disp(time_points);

%% 4. Compute pairwise propagation distances (mm) given velocity
N = size(electrode_coords,1);
speed = 300; % propagation speed mm/s
% Initialize distance matrix
distances = Inf(N);
for i = 1:N
    for j = setdiff(1:N,i)
        dt = abs(time_points(i) - time_points(j)); % seconds
        distances(i,j) = dt * speed;
    end
end

%% 5. Trilaterate origins for all 3-contact combinations within 5 mm
combs  = nchoosek(1:N,3);
origin_points = [];
for idx = 1:size(combs,1)
    % Select triplet
    e = combs(idx,:);
    p = electrode_coords(e,:);
    r12 = distances(e(1), e(2));
    r23 = distances(e(2), e(3));
    r31 = distances(e(3), e(1));
    % Trilateration on 2D flatmap
    origin = trilaterate(p(1,:),r12, p(2,:),r23, p(3,:),r31);
    % Check validity and geodesic proximity (<=5 mm)
    if ~any(isnan(origin)) && isWithin5mmOfAnyVertex(origin,flatmap)
        origin_points(end+1,:) = origin; %#ok<SAGROW>
        fprintf('Combo %3d → origin [%.2f, %.2f]\n', idx, origin);
    end
end

%% 6. Plot electrodes and origin points on flatmap
figure; plot(flatmap.x, flatmap.y,'.k','DisplayName','Flatmap vertices'); hold on;
scatter(electrode_coords(:,1), electrode_coords(:,2), 50,'b','filled','DisplayName','Electrodes');
scatter(origin_points(:,1), origin_points(:,2),100,'r','filled','DisplayName','Origins');
xlabel('X (flatmap)'); ylabel('Y (flatmap)'); title('Electrodes and Origins on Flatmap');
legend; grid on; hold off;

%% 7. Project origin points back to nearest pial vertices (3D)
nOrig = size(origin_points,1);
assert(nOrig>0,'No origin points found; check trilateration.');
vno   = zeros(nOrig,1);
distv = zeros(nOrig,1);
for i = 1:nOrig
    diffs = vertex_pial - [origin_points(i,:), zeros(1)]; % assume Z=0 on flat
    d2    = sqrt(sum(diffs(:,1:2).^2,2));
    [distv(i), idx] = min(d2);
    vno(i) = flatmap.ind(idx)+1;
end

%% 8. Compute Gaussian-blurred density on pial mesh
sigma = 10; % mm standard deviation
V    = vertex_pial;
D    = pdist2(V, V(vno,:));           % distance matrix
W    = exp(-D.^2/(2*sigma^2));        % Gaussian weights
vertex_density = sum(W,2);             % density per vertex
% Normalize to [0,1]
vertex_density = (vertex_density - min(vertex_density)) / (max(vertex_density)-min(vertex_density));

%% 9. Plot heatmap on pial surface
figure; h=patch('Vertices',V,'Faces',faces_pial, ...
    'FaceVertexCData',vertex_density,'FaceColor','interp', ...
    'EdgeColor','none','FaceAlpha',0.8);
view(3); axis equal tight off; colormap(jet); colorbar;
camlight headlight; lighting gouraud; material shiny;
% Annotate axes
lims = axis; mid = mean(lims);
text(lims(2),mid(3),mid(5),'Lateral','FontWeight','bold');
text(lims(1),mid(3),mid(5),'Medial','FontWeight','bold');
text(mid(1),lims(4),mid(5),'Anterior','FontWeight','bold');
text(mid(1),lims(3),mid(5),'Posterior','FontWeight','bold');

%% 10. Compute clustering coherency (2D and 3D)
cent2d = mean(origin_points,1);
cent3d = mean(V(vno,:),1);
avg2d  = mean(vecnorm(origin_points-cent2d,2,2));
avg3d  = mean(vecnorm(V(vno,:)-cent3d,2,2));
cov2d  = cov(origin_points); [V2d,E2d] = eig(cov2d);
cov3d  = cov(V(vno,:)); [V3d,E3d] = eig(cov3d);
disp('2D centroid, avg dist, eigenvals/vects:'); disp(cent2d); disp(avg2d); disp(diag(E2d)); disp(V2d);
disp('3D centroid, avg dist, eigenvals/vects:'); disp(cent3d); disp(avg3d); disp(diag(E3d)); disp(V3d);


