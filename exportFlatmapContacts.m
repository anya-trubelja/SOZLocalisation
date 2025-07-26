%%========================================================================
% Script:     exportFlatmapContacts.m
% Purpose:    Read a Freesurfer flat‐patch file, map electrode contacts to  
%             pial‐surface vertices within that patch, visualize the patch  
%             border and contacts on the inflated surface, and export the  
%             filtered contact list to Excel.
%
% Author:     Anya Trubelja
% Created:    2025‑05‑17
% Last Edit:  2025‑07‑26
%
% Requirements:
%   • MATLAB
%   • Freesurfer installation with MATLAB interface (FREESURFER_HOME)
%   • computeElectrodeContactTable.m in MATLAB path
%   • freesurferMatlabLibrary‑master in MATLAB path
%
% Inputs (user‑configurable at top of script):
%   subject_id      - Subject identifier string (e.g. 'sub20')
%   hemi            - Hemisphere specifier: 'lh' or 'rh'
%   patchFile       - Full path to .patch file for flatmap region
%   electrodeXLSX   - Path to electrode location Excel file (MRI‐centered)
%   outputXLSX      - Desired output path for filtered contacts Excel file
%
% Outputs:
%   • Figure showing inflated surface with:
%       - Patch region (green)
%       - Patch border edges (blue)
%       - All electrode contacts (red)
%   • Excel spreadsheet of contacts within the patch, with columns:
%       VertexIdx, Electrode, Contact, x, y, z, isgrey
%
% Workflow Overview:
%   0. Set up Freesurfer environment and MATLAB paths
%   1. Read patch file: extract vertex indices defining region
%   2. Compute mapping of electrode contacts → pial vertices
%   3. Filter contacts to those inside the patch
%   4. Read inflated surface, build border edges around patch
%   5. Plot surface, patch, border, and contacts for QC
%   6. Write filtered contact table to Excel
%
%%========================================================================

%% USER SETTINGS — adjust these parameters to your data
subject_id    = 'sub20';  % Subject identifier
hemi          = 'rh';     % Hemisphere: 'lh' or 'rh'
patchFile     = '/Applications/freesurfer/subjects/sub20/surf/rh.dkt.electrode.flat.patch';  % Flatmap .patch file path
electrodeXLSX = '/Users/neurallabpostgrads/Desktop/Electrode_locations/Electrodes_020/MRIcenteredHP.xlsx';  % Input electrode locations
outputXLSX    = '/Users/neurallabpostgrads/Desktop/Flatmap_Contacts_HP.xlsx';  % Output for filtered contacts

%% 0. Freesurfer environment setup
% Ensure FREESURFER_HOME is set, source setup script, and add MATLAB library
fsHome = getenv('FREESURFER_HOME');
if isempty(fsHome)
    fsHome = '/Applications/freesurfer';
    setenv('FREESURFER_HOME', fsHome);
end
system(sprintf('bash -lc "source %s/SetUpFreeSurfer.sh"', fsHome));

subjects_dir = getenv('SUBJECTS_DIR');
if isempty(subjects_dir)
    subjects_dir = '/Applications/freesurfer/subjects';
    setenv('SUBJECTS_DIR', subjects_dir);
end
subj_dir = fullfile(subjects_dir, subject_id);
surf_dir = fullfile(subj_dir, 'surf');

% Add Freesurfer MATLAB toolbox for surface IO
fstoolbox = 'Toolboxes/freesurferMatlabLibrary-master';
addpath(genpath(fstoolbox));

%% 1. Read patch file
% Open binary patch: first int32 is magic, second is vertex count
fid = fopen(patchFile, 'rb', 'ieee-be');
if fid < 0
    error('Cannot open patch file: %s', patchFile);
end

fread(fid, 1, 'int32');        % Skip magic number
nV = fread(fid, 1, 'int32');   % Number of vertices in patch
patchIdx = zeros(nV, 1);

% Read each entry: index (int32) then x,y,z (single) which we skip
for i = 1:nV
    idx = fread(fid, 1, 'int32');
    patchIdx(i) = abs(idx);     % Use absolute to correct sign
    fread(fid, 3, 'single');    % Discard coordinates
end
fclose(fid);

% Convert 0-based indices if present to MATLAB 1-based
if any(patchIdx == 0)
    patchIdx = patchIdx + 1;
end

fprintf('Patch vertices: %d (indices range %d – %d)\n', nV, min(patchIdx), max(patchIdx));

%% 2. Map electrodes to pial vertices
% computeElectrodeContactTable returns table with contact→vertex mapping
T = computeElectrodeContactTable(subj_dir, hemi, electrodeXLSX);
if isempty(T)
    warning('No grey-matter contacts found. Exiting.');
    return;
end
fprintf('Total grey-matter contacts: %d (vert %d – %d)\n', ...
        height(T), min(T.PialVertexIdx), max(T.PialVertexIdx));

%% 3. Filter contacts inside patch
inMask = ismember(T.PialVertexIdx, patchIdx);  % Logical mask of contacts within ROI
fprintf('Contacts inside patch: %d\n', sum(inMask));
T2 = T(inMask, :);  % Subset table

%% 4. Load inflated surface and determine patch border
% Read inflated mesh vertices (Vi) and faces (Fi_raw)
surfFile = fullfile(surf_dir, sprintf('%s.inflated', hemi));
[Vi, Fi_raw] = freesurfer_read_surf(surfFile);
Fi = Fi_raw + (min(Fi_raw(:)) == 0);  % Ensure faces are 1-based indices

% Sanity check: patch indices must be within vertex count
assert(max(patchIdx) <= size(Vi,1), 'Patch index exceeds mesh vertex count.');

% Build all edges from faces (each triangle gives 3 edges)
e12 = [Fi(:,1), Fi(:,2)];
e23 = [Fi(:,2), Fi(:,3)];
e31 = [Fi(:,3), Fi(:,1)];
allE = sort([e12; e23; e31], 2);  % Sort each row for consistency

% Identify edges with exactly one endpoint in patch (border)
in1 = ismember(allE(:,1), patchIdx);
in2 = ismember(allE(:,2), patchIdx);
borderEdges = allE(xor(in1, in2), :);
fprintf('Border edges count: %d\n', size(borderEdges,1));

%% 5. Visualization: patch, border, contacts
figure('Name', 'Patch-on-Inflated Debug', 'Color', 'w');
hold on;

% Draw semi-transparent inflated surface mesh
patch('Vertices', Vi, 'Faces', Fi, ...
      'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Overlay patch region (green markers)
scatter3(Vi(patchIdx,1), Vi(patchIdx,2), Vi(patchIdx,3), ...
         15, 'g', 'filled', 'MarkerFaceAlpha', 0.4);

% Plot border edges (blue lines)
for i = 1:size(borderEdges,1)
    v = borderEdges(i,:);
    plot3(Vi(v,1), Vi(v,2), Vi(v,3), 'b-', 'LineWidth', 1.5);
end

% Overlay all electrode contacts (red markers)
CA = Vi(T.PialVertexIdx, :);
scatter3(CA(:,1), CA(:,2), CA(:,3), 60, 'r', 'filled');

axis equal off;
camlight headlight;
view(3);
title(sprintf('%s %s patch on inflated', subject_id, hemi), 'FontSize', 14);
hold off;

%% 6. Export filtered contact table to Excel
if ~isempty(T2)
    % Extract vertex indices and contact labels
    vertexIdx = T2.PialVertexIdx;
    E = T2.Shaft;    % Electrode shaft labels
    C = T2.Contact;  % Contact numbers

    % Retrieve XYZ coordinates for each contact vertex
    coords = Vi(vertexIdx, :);
    x = coords(:,1);
    y = coords(:,2);
    z = coords(:,3);

    % Flag all as grey-matter contacts
    isgrey = ones(height(T2), 1);

    % Build output table with named variables
    Out = table(vertexIdx, E, C, x, y, z, isgrey, ...
        'VariableNames', {'VertexIdx', 'Electrode', 'Contact', 'x', 'y', 'z', 'isgrey'});

    % Write to spreadsheet
    writetable(Out, outputXLSX, 'FileType', 'spreadsheet');
    fprintf('Wrote %d contacts → %s\n', height(Out), outputXLSX);
else
    warning('No contacts to write.');
end
