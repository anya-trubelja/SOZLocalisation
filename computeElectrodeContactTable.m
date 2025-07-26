function T = computeElectrodeContactTable(subjectDir, hemi, electrodeXLSX)
%%========================================================================
% Function:   computeElectrodeContactTable
% Purpose:    Load Freesurfer pial surface and electrode location Excel file,
%             filter for grey-matter contacts, project each contact to the
%             nearest pial vertex, and return a detailed table of results.
%
% Author:     Anya Trubelja
% Created:    2025-05-14
% Last Edit:  2025-07-26
%
% Requirements:
%   • MATLAB
%   • Freesurfer MATLAB toolbox (freesurfer_read_surf)
%   • Electrode locations XLSX with columns:
%       Shaft (string), Contact (numeric), X, Y, Z (double), isgrey (0/1)
%
% Usage:
%   T = computeElectrodeContactTable(subjectDir, hemi, electrodeXLSX);
%
% Inputs:
%   subjectDir    - Full path to Freesurfer subject directory
%   hemi          - Hemisphere code: 'lh' or 'rh'
%   electrodeXLSX - Path to MRI-centered electrode locations Excel file
%
% Output:
%   T : table with columns [nGrey×9]:
%       Shaft, Contact, X_before, Y_before, Z_before,
%       PialVertexIdx, X_proj, Y_proj, Z_proj
%%========================================================================

    %% 1) Load pial surface mesh
    pialSurfFile = fullfile(subjectDir, 'surf', sprintf('%s.pial', hemi));
    if ~exist(pialSurfFile, 'file')
        error('computeElectrodeContactTable: missing pial surface: %s', pialSurfFile);
    end
    % Read vertices and faces from Freesurfer format
    [vertex_pial_raw, faces_pial_raw] = freesurfer_read_surf(pialSurfFile);
    % Ensure face indices are 1-based
    if min(faces_pial_raw(:)) == 0
        faces_pial = faces_pial_raw + 1;
    else
        faces_pial = faces_pial_raw;
    end
    % Validate face indices range
    nVert = size(vertex_pial_raw, 1);
    if any(faces_pial(:) < 1) || any(faces_pial(:) > nVert)
        error('Face indices out of bounds [1..%d]', nVert);
    end
    vertex_pial = vertex_pial_raw;  % Alias for clarity
    fprintf('Loaded pial surface: %d vertices, %d faces.\n', nVert, size(faces_pial,1));

    %% 2) Read electrode Excel file
    % Read raw cell array; expect at least header + 6 columns
    C = readcell(electrodeXLSX);
    [nRows, nCols] = size(C);
    if nRows < 2 || nCols < 6
        error('Excel file must have ≥1 header + data rows and ≥6 columns');
    end
    % Parse columns: shaft names, contact numbers, XYZ coords, grey flag
    shaftAll     = string(C(2:end,1));
    contactAll   = cell2mat(C(2:end,2));
    x_beforeAll  = cell2mat(C(2:end,3));
    y_beforeAll  = cell2mat(C(2:end,4));
    z_beforeAll  = cell2mat(C(2:end,5));
    isgreyAll    = cell2mat(C(2:end,6));

    %% 3) Filter for grey-matter contacts
    idxGrey = find(isgreyAll == 1);
    nGrey = numel(idxGrey);
    if nGrey == 0
        warning('No grey-matter contacts found. Returning empty table.');
        T = table(string.empty, [], [], [], [], [], [], [], [], ...
            'VariableNames', {'Shaft','Contact','X_before','Y_before','Z_before', ...
                              'PialVertexIdx','X_proj','Y_proj','Z_proj'});
        return;
    end
    % Prepare arrays for filtered contacts
    shafts   = shaftAll(idxGrey);
    contacts = contactAll(idxGrey);
    Xb       = x_beforeAll(idxGrey);
    Yb       = y_beforeAll(idxGrey);
    Zb       = z_beforeAll(idxGrey);
    Pidx     = zeros(nGrey,1);
    Xp       = zeros(nGrey,1);
    Yp       = zeros(nGrey,1);
    Zp       = zeros(nGrey,1);

    %% 4) Project each contact to nearest pial vertex
    for k = 1:nGrey
        pt = [Xb(k), Yb(k), Zb(k)];
        % Compute Euclidean distances to all vertices
        diffs = vertex_pial - pt;
        dists = sqrt(sum(diffs.^2,2));
        [~, vid] = min(dists);
        % Record nearest vertex index and its coordinates
        Pidx(k) = vid;
        Xp(k)   = vertex_pial(vid,1);
        Yp(k)   = vertex_pial(vid,2);
        Zp(k)   = vertex_pial(vid,3);
    end
    fprintf('Projected %d contacts onto pial surface.\n', nGrey);

    %% 5) Assemble and return result table
    T = table(shafts, contacts, Xb, Yb, Zb, Pidx, Xp, Yp, Zp, ...
        'VariableNames', {'Shaft','Contact','X_before','Y_before','Z_before', ...
                          'PialVertexIdx','X_proj','Y_proj','Z_proj'});
end
