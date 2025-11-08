% PC-MRI-GE loading, masking, and displacement calculation
% Written by Mahsa Karamzadeh
% Last updated: 2025-09-27
%
% This script loads pcMRI data from a .mat file (with variable "processed_data"),
% computes timing of the cardiac phases, and prepares data for further analysis.


close all; clc; clear;
%% ----------------------------- Setup & Params -----------------------------

% Palette (for consistent plots)
Yellow = [0.99 0.99 0.00];
Blue   = [0.00 0.45 0.74];
Red    = [0.85 0.33 0.10];
Green  = [0.47 0.67 0.19];

% Figure layout placeholders 
width = 0.33; height = 0.41; fig1x = 0; fig1y = 0.05; padx = 0.0; pady = 0.1; 

% ---- User-adjustable parameters ----
VENC = 2;         % velocity encoding [cm/s] (keep units explicit)
heartRateBpm = 72;     % beats per minute used for timing

% Voxel size [mm] 
vsz_x = 0.6; vsz_y = 0.6; vsz_z = 3.0;

%% ----------------------------- Select Data -------------------------------

[matFile, matPath] = uigetfile({'*.mat','MAT-files (*.mat)'}, 'Select pcMRI data');
if isequal(matFile,0)
    error('No file selected.');
end
filePath = fullfile(matPath, matFile);

% Default save directory: a "pcMRI_Result" folder next to the data file
defaultSaveDir = fullfile(matPath, 'pcMRI_Result');
if ~exist(defaultSaveDir, 'dir')
    mkdir(defaultSaveDir);
end
saveDir = defaultSaveDir;  

%% ----------------------------- Subject ID --------------------------------
% Try to infer a subject/patient ID from parent folder names.
% Your data format may differ; tweak as needed.
[parent1, ~, ~] = fileparts(matPath);       % parent folder
[parent2, lastFolder, ~] = fileparts(parent1);  % grandparent name
if ~isempty(lastFolder)
    Patient_ID = strrep(lastFolder, '_', ' ');
else
    [~, stem, ~] = fileparts(filePath);
    Patient_ID = stem;
end

%% ----------------------------- Load Data ---------------------------------
S = load(filePath);
if ~isfield(S, 'processed_data')
    error('Variable "processed_data" not found in %s', filePath);
end

% Convert velocity to cm/s (assuming file stores mm/s; change if needed)
vel = double(S.processed_data) / 10;   % mm/s → cm/s

% Dimensions: accept R x C x T or R x C x Z x T
nd = ndims(vel);
if nd == 3
    [sz_r, sz_c, NoPh] = size(vel); 
    nSlices = 1;
elseif nd == 4
    [sz_r, sz_c, nSlices, NoPh] = size(vel); 
else
    error('processed_data must be 3D or 4D (R x C x [Z] x T). Got ndims=%d.', nd);
end

%% ----------------------------- Timing ------------------------------------
RR_sec  = 60 / heartRateBpm;              % seconds per cardiac cycle
del_t   = RR_sec / NoPh;                  % seconds between consecutive phases
t_phases = 1:NoPh;                        % phase indices (1..NoPh)
t_sec    = del_t * t_phases;              % timestamps [s], start-of-bin

%% ----------------------------- Sanity Print ------------------------------
fprintf('Loaded: %s\n', filePath);
fprintf('Patient ID: %s\n', Patient_ID);
fprintf('Size: %dx%dx%dx%d (R x C x Z x T)\n', sz_r, sz_c, nSlices, NoPh);
fprintf('Velocity range (cm/s): [%.2f, %.2f]\n', min(vel,[],'all'), max(vel,[],'all'));
fprintf('RR = %.3f s   Δt = %.3f s   phases = %d   VENC = %.2f cm/s\n', RR_sec, del_t, NoPh, VENC);
fprintf('Save dir: %s\n', saveDir);

%% ----------- Correcting inhomogeneity in the velocity field -------------

% Assumes vel is 3-D: [Ny × Nx × Nt]
[Ny, Nx, Nt] = size(vel);

% --- 1) Compute per-voxel temporal mean (background bias) ---
mean_vel = mean(vel, 3);           % Ny×Nx

% --- 2) Subtract bias from each time-frame ---
vel_corr = vel - mean_vel;         % Ny×Nx×Nt

% --- 3) Quick check: mean after correction (should be near zero) ---
mean_vel_corr = mean(vel_corr, 3); 

% --- Parameters for display ---
phase    = 1;          % which time-frame to show
limits   = [ -2   2;   % color limits for first figure
            -0.3  0.3];% color limits for second figure
fontSize = 14;         % font size for titles and labels

% --- Plot before/after correction for the chosen phase ---
for i = 1:size(limits,1)
    clim = limits(i,:);

    figure('Name', sprintf('Phase %d  (caxis = [%g %g])', phase, clim(1), clim(2)), ...
           'Position', [150 150 900 400]);

    % BEFORE correction
    subplot(1,2,1);
    imagesc(vel(:,:,phase));
    axis equal tight off;
    caxis(clim);
    title(sprintf('Phase %d  BEFORE', phase), 'FontSize', fontSize);
    cb = colorbar;
    cb.Label.String   = 'Velocity [cm/s]';
    cb.Label.FontSize = fontSize;
    set(gca,'FontSize',fontSize);

    % AFTER correction
    subplot(1,2,2);
    imagesc(vel_corr(:,:,phase));
    axis equal tight off;
    caxis(clim);
    title(sprintf('Phase %d  AFTER', phase), 'FontSize', fontSize);
    cb = colorbar;
    cb.Label.String   = 'Velocity [cm/s]';
    cb.Label.FontSize = fontSize;
    set(gca,'FontSize',fontSize);

    colormap(parula);
end


%% ---------------- View all phases & choose a phase for masking ----------

% Ask user if they want to see the whole cardiac-cycle movie
yn = questdlg('Would you like to see the phase movie', 'Movie?', {'Yes', 'No'});


% Create figure for movie
fig1 = figure(1);
set(fig1,'Units','normalized','Position',[0.25 0.25 0.5 0.5]);

if strcmpi(yn,'Yes')
    ax = axes('Parent',fig1,'Position',[0.05 0.05 0.9 0.9]);  % leave room for title
    for k = 1:NoPh
        imshow(vel(:,:,k),[-VENC VENC],'Parent',ax);
        title(sprintf('Phase = %d',k),'FontWeight','bold');
        axis off
        pause(0.3);               % pause to create movie effect
    end
end

% Prompt user for phase number (default = 17)
if strcmpi(yn,'Yes') || strcmpi(yn,'No')
    answer = inputdlg({'Enter phase # to plot phase image'}, ...
                      'Phase Input', [1 40], {'17'});
    ik = round(str2double(answer{1}));
else
    ik = 17;
end
ik = max(1, min(NoPh, ik));       % keep phase within valid range

%% ---------------- Show selected phase & choose ROI ----------------------

fig2 = figure(2);
set(fig2,'Units','normalized','Position',[0.25 0.25 0.5 0.5]);

ax = axes('Parent',fig2,'Position',[0 0 1 1]);
imshow(vel(:,:,ik),[-VENC VENC],'Parent',ax);
title(sprintf('Phase = %d',ik),'FontWeight','bold');
axis off

% Let user adjust contrast interactively
imcontrast

% Let user click two points: top-left and bottom-right of ROI
title('Click 2 points: top-left and bottom-right of ROI');

roi = drawpoint;  pos = roi.Position;  tly = round(pos(2));  tlx = round(pos(1));
roi = drawpoint;  pos = roi.Position;  bry = round(pos(2));  brx = round(pos(1));

close(fig2);

%%
% Main loop
while true
fig_nums = [1,2, 3, 4, 5, 6,7, 8, 9];
for fn = fig_nums
    if isgraphics(fn)
        try, close(fn); end
    end
end
clf;
    clear roi_poly pos_poly mask_poly mask vel_mask vel_mask2 displacement_circ
    clear cin rin roi_vel_avg_circle roi_disp_circle_temp

    % s is three layers of the velocity scaled to be [0,1] to make it RGB later
    s(:,:,1) = (vel(:,:,ik)+VENC)/(2*VENC);
    s(:,:,2) = (vel(:,:,ik)+VENC)/(2*VENC);
    s(:,:,3) = (vel(:,:,ik)+VENC)/(2*VENC);

    fig2 = figure(2);
    fig2.Units = "normalized";
    fig2.Position = [fig1x fig1y+height+pady-0.05 width height];
    imagesc(s(tly:bry, tlx:brx, :));
    title(['Phase = ', num2str(ik)]);

    roi_poly = drawpolygon;
    pos_poly = roi_poly.Position;
    xadjusted_poly = pos_poly(:,1) + tlx; % because we cropped the image 
    yadjusted_poly = pos_poly(:,2) + tly;
    mask_poly = round(pos_poly,0) + [tlx - 1, tly - 1];
    mask = poly2mask(mask_poly(:,1), mask_poly(:,2), sz_r, sz_c);

    % To force the final displacement to be zero
    corrected_vel = vel;
    mean_vel_all = mean(corrected_vel, 3);
    corrected_vel = corrected_vel - mean_vel_all;

       for i = 1:sz_r
        for j = 1:sz_c
            if mask(i, j) == 1
                s(i, j, :) = Green;
            end
        end
    end

    mmcircrad = 1.5;
    numit = 50;  % number of points on the circle
    theta = 0:2*pi/numit:2*pi;
    pixelsize = 0.6;
    cirpix = mmcircrad / pixelsize;

    % Finding x and y cordinates of ones in mask
    % row is y and column is x
    [rones, cones] = find(mask);  %coordinates are based on actual image 

    for C = 1:nnz(mask)
        vel_mask(C,:) = corrected_vel(rones(C), cones(C), :);

        % Defining 50 points on the edge of the circle
        for ii = 1:numit
            Xcirc(ii,C) = cirpix*cos(theta(ii)) + cones(C);
            Ycirc(ii,C) = cirpix*sin(theta(ii)) + rones(C);
        end
    end

    % Checking if all those 50 points on each circles are inside the mask
    d = 0; in = zeros(numit, nnz(mask));
    for b = 1:nnz(mask)
        in(:,b) = inpolygon(Xcirc(:,b), Ycirc(:,b), xadjusted_poly, yadjusted_poly);

        if nnz(in(:,b)) > numit / 1.5  % If only numit/1.5 is inside the polygone that is accepteble
            d = d + 1;
            cin(d) = cones(b);
            rin(d) = rones(b);  % rones is row which is y
            xcircle_small(:,d) = Xcirc(:,b);
            ycircle_small(:,d) = Ycirc(:,b);
            vel_mask2(d,:) = vel_mask(b,:);
        end
    end

    figure(6); plot(xadjusted_poly,yadjusted_poly); hold on;
    plot(cin, rin, '.'); axis equal; set(gca, 'YDir','reverse'); hold off;

    % For each small circle, check which pixels are inside those small circles
    for v = 1:size(cin,2)
        in_circle = inpolygon(cones, rones, xcircle_small(:,v), ycircle_small(:,v));
        [in_circle_row,~] = find(in_circle);
        for i = 1:NoPh
            roi_vel_avg_circle(i,v) = mean(vel_mask(in_circle_row,i));
        end
        mean_vel_circle(v) = mean(roi_vel_avg_circle(:,v));  % Temporal avg of velocity [cm/s]
        roi_disp_circle_temp(:,v) = 10000 * cumtrapz(del_t, roi_vel_avg_circle(:,v));
        displacement_circ(v) = max(roi_disp_circle_temp(:,v)) - min(roi_disp_circle_temp(:,v));
    end

    [~,v_row_max] = max(displacement_circ);
    [~,v_row_min] = min(displacement_circ);
    fig4 = figure(4); imagesc(s(tly:bry, tlx:brx, :)); hold on;
    fig4.Units = "normalized";
    fig4.Position = [fig1x + width + padx, fig1y, width, height];
    plot(xcircle_small(:,v_row_max)-tlx, ycircle_small(:,v_row_max)-tly, 'r');
    plot(xcircle_small(:,v_row_min)-tlx, ycircle_small(:,v_row_min)-tly, 'b');
    axis equal;

    dismat = zeros(size(mask));
    for ppp = 1:size(cin,2)
        dismat(rin(ppp), cin(ppp)) = displacement_circ(ppp);
    end
    fig7=figure(7); imagesc(dismat); caxis([0, 900]); colorbar; axis equal;
    fig7.Units = "normalized"; fig7.Position = [fig1x fig1y width height];

    for q = 1:NoPh
        roi_corrected_vel_avg(q) = sum(mask .* corrected_vel(:,:,q), 'all') / nnz(mask);
        roi_vel_max(q) = max((mask .* vel(:,:,q)), [], 'all');
        roi_vel_min(q) = min((mask .* vel(:,:,q)), [], 'all');
    end
    roi_disp_max = 10000 * max(cumtrapz(del_t, roi_corrected_vel_avg));
    roi_disp_min = 10000 * min(cumtrapz(del_t, roi_corrected_vel_avg));

    figure(5);
    for t = 1:min(30, NoPh)
        subplot(5,6,t);
        imagesc(corrected_vel(tly:bry, tlx:brx, t), [-0.2, 0.2]);
        colorbar;
        title(['phase=' num2str(t)]);
    end

    fig3 = figure(3);
    fig3.Units = "normalized";
    fig3.Position = [fig1x+2*width+2*padx fig1y width height];

    subplot(2,1,1);
    plot(t_phases, roi_corrected_vel_avg, 'g', 'LineWidth', 3);
    yline(0, 'k--');
    xlabel('t [phases]'); ylabel('V_{avg}(t) [cm/s]');
    title("Patient ID = " + Patient_ID);

    subplot(2,1,2);
    plot(t_phases, 10000 * cumtrapz(del_t, roi_corrected_vel_avg), '.-', 'Color', Blue, 'LineWidth', 3);
    yline(0, 'k--');
    xlabel('t [phases]'); ylabel('Displacement(t) [micron]');

    annotation('textbox', [0.15, 0.1, 0.1, 0.1], 'String', ...
        {['Displacement= ' num2str(roi_disp_max - roi_disp_min, '%.0f')]}, ...
        'FontSize', 25, 'FitBoxToText', 'on', 'LineStyle', 'none', 'FontWeight', 'bold');


% Ask user whether to save the data

save_response = questdlg('Do you want to save the data for this region?', ...
                         'Save Data?', {'Yes', 'No'});
if strcmp(save_response, 'Yes')
    % Ask which region this is
    drawnow;
pause(0.1);  % allow GUI to update
    region_list = {'Cerebellum', 'Pons', 'Medulla', 'Spinal Cord', 'Other'};
region_idx = menu('Select the brain region:', region_list);

if region_idx ~= 0  % if user made a choice
    region_name = region_list{region_idx};
    region_tag = strrep(lower(region_name), ' ', '_');

    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    save_prefix = fullfile(savedir, [Patient_ID, '_', region_tag, '_', timestamp]);

    % Save figures
    saveas(fig3, [save_prefix '_plot_velocity_displacement.png']);
    saveas(fig4, [save_prefix '_max_min_displacement.png']);
    saveas(fig7, [save_prefix '_displacement_map.png']);

    % Save Excel data
    param_table = table(roi_disp_max - roi_disp_min, ...
                        max(roi_corrected_vel_avg), ...
                        min(roi_corrected_vel_avg), ...
                        'VariableNames', {'Displacement', 'Vel_Max', 'Vel_Min'});

    excel_file = fullfile(savedir, 'results_summary.xlsx');
    sheet_name = regexprep([Patient_ID, '_', region_tag], '[^\w]', '_');

    % Save Excel data to unified table
excel_file = fullfile(savedir, 'results_summary.xlsx');

% New data row
new_entry = table({Patient_ID}, {region_name}, ...
                  roi_disp_max - roi_disp_min, ...
                  max(roi_corrected_vel_avg), ...
                  min(roi_corrected_vel_avg), ...
                  'VariableNames', {'Patient_ID', 'Region', 'Displacement', 'Vel_Max', 'Vel_Min'});

% Load existing data if file is valid
if isfile(excel_file)
    try
        existing_data = readtable(excel_file);

        % Check headers
        required_vars = {'Patient_ID', 'Region', 'Displacement', 'Vel_Max', 'Vel_Min'};
        if ~all(ismember(required_vars, existing_data.Properties.VariableNames))
            warning('⚠️ Existing file has wrong format. Overwriting.');
            existing_data = new_entry;
        else
            % Overwrite if same Patient_ID + Region, else append
            match_idx = strcmp(existing_data.Patient_ID, Patient_ID) & ...
                        strcmp(existing_data.Region, region_name);

            if any(match_idx)
                existing_data(match_idx, :) = new_entry;
            else
                existing_data = [existing_data; new_entry];
            end
        end
    catch
        warning('⚠️ Could not read existing Excel. Replacing with new.');
        existing_data = new_entry;
    end
else
    existing_data = new_entry;
end

% Write final table
writetable(existing_data, excel_file);
disp(['✅ Saved/updated: ' region_name ' for ' Patient_ID]);

else
    disp('⚠️ Region selection canceled. Nothing was saved.');
end
else
disp('User did not want to save')
end
    % Ask user if they want to terminate or redo
    fig2 = figure(2);
    title("Click on the bottom right for re-do and top left for termination",FontSize=11)
    roi = drawpoint;
    Pos = roi.Position;
    ryi = round(Pos(2)) + tly - 1;
    rxi = round(Pos(1)) + tlx - 1;

    % To terminate click the pixel in the top left corner
    if rxi == tlx && ryi == tly
        disp('Terminated by user.');
        break;
    % To re-do segmentation
    elseif rxi == brx && ryi == bry
        disp('Redoing segmentation...');
        continue;
    end
end
