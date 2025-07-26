%%========================================================================
% Script:     compare_mean_distortion_from_figs.m
% Purpose:    Load two flatmap figures (.fig) containing signed-percent
%             distortion values, extract distortion data, compute summary
%             statistics on absolute distortion (mean & std), and display
%             results in both text and a bar graph with error bars.
%
% Author:     Anya Trubelja
% Created:    2025-12-25
% Last Edit:  2025-07-26
%
% Requirements:
%   • MATLAB with Statistics and Machine Learning Toolbox
%   • Two .fig files each showing signed-percent distortion heatmap
%
% Usage:
%   1) Set figFile1 and figFile2 under “%% 1. Specify .fig filenames.”
%   2) Run this script in MATLAB:
%         >> compare_mean_distortion_from_figs
%   3) Review printed mean ± std values in Command Window and bar graph.
%%========================================================================

%% 1. Specify your two .fig filenames here
figFile1 = 'electrode_distort_020.fig';       % Generated flatmap distortion
figFile2 = 'Geodesdic_distort_sub_020.fig';   % Manual flatmap distortion

%% 2. Helper function: extract signed-percent distortion vector from .fig
function D = getDistortionVectorFromFig(fname)
    % Open figure invisibly to avoid display
    hFig = openfig(fname, 'invisible');
    % Find the first axes object
    ax = findobj(hFig, 'Type', 'axes');
    if isempty(ax)
        close(hFig);
        error('No axes found in "%s".', fname);
    end
    ax = ax(1);

    % 2.1) Try image objects (imagesc, imshow)
    hImg = findobj(ax, 'Type', 'Image');
    if ~isempty(hImg)
        C = get(hImg(1), 'CData');
        D = C(:);
        close(hFig);
        return;
    end

    % 2.2) Try surface objects (surf, trisurf)
    hSurf = findobj(ax, 'Type', 'Surface');
    if ~isempty(hSurf)
        cdata = get(hSurf(1), 'CData');
        if ~isempty(cdata)
            D = cdata(:);
        else
            zdata = get(hSurf(1), 'ZData');
            D = zdata(:);
        end
        close(hFig);
        return;
    end

    % 2.3) Try patch objects (flatmap irregular grid)
    hPatch = findobj(ax, 'Type', 'Patch');
    if ~isempty(hPatch)
        fv = get(hPatch(1), 'FaceVertexCData');
        if ~isempty(fv)
            D = fv(:);
        else
            cd = get(hPatch(1), 'CData');
            if ~isempty(cd)
                D = cd(:);
            else
                close(hFig);
                error('Patch found but no FaceVertexCData or CData in "%s".', fname);
            end
        end
        close(hFig);
        return;
    end

    % 2.4) Try scatter objects
    hSc = findobj(ax, 'Type', 'Scatter');
    if ~isempty(hSc)
        D = get(hSc(1), 'CData');
        D = D(:);
        close(hFig);
        return;
    end

    % 2.5) Try contour objects
    hCon = findobj(ax, 'Type', 'Contour');
    if ~isempty(hCon)
        z = get(hCon(1), 'ZData');
        if ~isempty(z)
            D = z(:);
        else
            cd = get(hCon(1), 'CData');
            D = cd(:);
        end
        close(hFig);
        return;
    end

    % If none matched, throw an error
    close(hFig);
    error('No Image/Surface/Patch/Scatter/Contour found in "%s".', fname);
end

%% 3. Extract the two signed-distortion vectors
try
    D1_signed = getDistortionVectorFromFig(figFile1);
catch ME
    error('Failed to extract from %s: %s', figFile1, ME.message);
end
try
    D2_signed = getDistortionVectorFromFig(figFile2);
catch ME
    error('Failed to extract from %s: %s', figFile2, ME.message);
end
% Ensure numeric double precision
D1_signed = double(D1_signed);
D2_signed = double(D2_signed);

%% 4. Compute absolute distortion values
A1 = abs(D1_signed);
A2 = abs(D2_signed);

%% 5. Summary statistics: mean(|D|), std(|D|), sample size
n1     = numel(A1);
meanA1 = mean(A1);
stdA1  = std(A1);

n2     = numel(A2);
meanA2 = mean(A2);
stdA2  = std(A2);

% Display numerical results in Command Window
fprintf('\n=== Absolute Distortion Summary ===\n');
fprintf('Generated Flatmap: N = %6d   mean(|D|) = %8.4f   std(|D|) = %8.4f\n', n1, meanA1, stdA1);
fprintf('Manual Flatmap   : N = %6d   mean(|D|) = %8.4f   std(|D|) = %8.4f\n', n2, meanA2, stdA2);

%% 6. Bar graph with error bars ± STD
figure('Name','Mean Absolute Distortion','Color','w');
means = [meanA1, meanA2];
stds  = [stdA1,  stdA2];
b = bar(1:2, means, 0.5, 'FaceColor', [0.2 0.6 0.8]);
hold on;
errorbar(1:2, means, stds, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
hold off;

% Configure axes and labels
xlim([0.5 2.5]);
set(gca, 'XTick', [1 2], 'XTickLabel', {'Generated','Manual'});
ylabel('Mean |Distortion| (signed %)');
title('Comparison of Mean Absolute Distortion ± STD');
grid on;

% Annotate numeric mean±std on each bar
y_offset = 0.02 * max(means + stds);
for i = 1:2
    text(i, means(i) + stds(i) + y_offset, ...
         sprintf('%.2f ± %.2f', means(i), stds(i)), ...
         'HorizontalAlignment','center', 'FontSize',10);
end

