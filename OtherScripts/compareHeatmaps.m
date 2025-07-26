%%========================================================================
% Script:     compareHeatmaps.m
% Purpose:    Load two surface‐based heatmap .fig files, extract scalar intensity
%             maps, compute their difference and correlation statistics, and
%             visualize the original maps, their difference, and a scatter plot.
%
% Author:     Anya Trubelja
% Created:    2025-07-12
% Last Edit:  2025-07-26
%
% Requirements:
%   • MATLAB
%   • Two .fig files each containing a Patch object with FaceVertexCData
%
% Inputs (set at top of script):
%   fig1_file - Filename of first heatmap .fig (automated)
%   fig2_file - Filename of second heatmap .fig (manual)
%
% Outputs:
%   • Console output of mean difference, max absolute difference, and Pearson R
%   • Figure window showing:
%       1) Map 1
%       2) Map 2
%       3) Difference map (Map2 − Map1)
%       4) Scatter plot of Map1 vs. Map2 intensities colored by difference
%%========================================================================

%% User parameters — adjust paths to your .fig files
fig1_file = 'SOZHeatmap023.fig';        % Automated SOZ heatmap
fig2_file = 'ManualSOZHeatmap023.fig';  % Manual SOZ heatmap

%% 1. Load figures invisibly
% Open without displaying to extract patch data
hF1 = openfig(fig1_file, 'invisible');
hF2 = openfig(fig2_file, 'invisible');

%% 2. Locate the Patch objects
p1 = findobj(hF1, 'Type', 'Patch');  % Should contain surface geometry + CData
p2 = findobj(hF2, 'Type', 'Patch');
assert(~isempty(p1) && ~isempty(p2), 'Patch objects not found in one or both FIG files.');

%% 3. Extract geometry and raw color data
% Vertices (V), Faces (F), raw CData for each patch
V1 = p1.Vertices;  F1 = p1.Faces;  C1_raw = p1.FaceVertexCData;
V2 = p2.Vertices;  F2 = p2.Faces;  C2_raw = p2.FaceVertexCData;

%% 4. Convert raw CData to scalar intensity maps
% If CData is Nx3 RGB, convert via luminance; else if Nx1, use directly
if size(C1_raw,2)==3 && size(C1_raw,1)==size(V1,1)
    C1 = 0.299*C1_raw(:,1) + 0.587*C1_raw(:,2) + 0.114*C1_raw(:,3);
elseif numel(C1_raw)==size(V1,1)
    C1 = double(C1_raw(:));
else
    error('Cannot interpret CData of map1 (%dx%d).', size(C1_raw));
end

if size(C2_raw,2)==3 && size(C2_raw,1)==size(V2,1)
    C2 = 0.299*C2_raw(:,1) + 0.587*C2_raw(:,2) + 0.114*C2_raw(:,3);
elseif numel(C2_raw)==size(V2,1)
    C2 = double(C2_raw(:));
else
    error('Cannot interpret CData of map2 (%dx%d).', size(C2_raw));
end

%% 5. Compute difference map and summary statistics
dMap  = C2 - C1;                   % Difference per vertex
meanD = mean(dMap);                % Mean difference
maxD  = max(abs(dMap));            % Maximum absolute difference
R     = corr(C1, C2);              % Pearson correlation

fprintf('Mean Δ (map2−map1): %.4f\n', meanD);
fprintf('Max |Δ|         : %.4f\n', maxD);
fprintf('Pearson R       : %.4f\n', R);

%% 6. Visualize maps and their comparison
figure('Name','Heat‐Map Comparison','Units','normalized','Position',[.1 .1 .8 .6]);

% 6a) Map 1
subplot(2,3,1);
patch('Vertices',V1,'Faces',F1, ...
      'FaceVertexCData',C1, 'FaceColor','interp','EdgeColor','none', ...
      'CDataMapping','scaled');
view(3); axis tight off; daspect([1 1 1]);
title('Map 1: Automated'); colormap(gca,'jet'); colorbar;

% 6b) Map 2
subplot(2,3,2);
patch('Vertices',V2,'Faces',F2, ...
      'FaceVertexCData',C2, 'FaceColor','interp','EdgeColor','none', ...
      'CDataMapping','scaled');
view(3); axis tight off; daspect([1 1 1]);
title('Map 2: Manual');   colormap(gca,'jet'); colorbar;

% 6c) Difference map (Map2 − Map1)
subplot(2,3,3);
patch('Vertices',V1,'Faces',F1, ...
      'FaceVertexCData',dMap, 'FaceColor','flat','EdgeColor','none', ...
      'CDataMapping','scaled');
view(3); axis tight off; daspect([1 1 1]);
title('Difference (2−1)');
lim = max(abs(dMap)); caxis([-lim lim]);
colormap(gca,'parula'); colorbar;

% 6d) Scatter plot of intensities colored by difference
subplot(2,3,[4 5 6]);
scatter(C1, C2, 20, dMap, 'filled');
xlabel('Map 1 intensity'); ylabel('Map 2 intensity');
title(sprintf('Scatter R=%.2f', R));
axis equal tight; grid on; colormap(gca,'jet'); colorbar;

%% 7. Close figure handles to free resources
close(hF1);
close(hF2);

