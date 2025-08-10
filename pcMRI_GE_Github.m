     %PC-MRI-GE loading and masking      08.22.2023

close all;  clear; clc;

orange_yellow = [224 163 46]/256;
orange = [210 120 44]/256;
purple = [94 60 108]/256;
Yellow = [255 255 0]/256;
Blue = [0 0.45 0.74];
Red = [0.85 0.33 0.10];
Green = [0.47 0.67 0.19];
width = 0.33;    height = 0.41; fig1x = 0; fig1y = 0.05; padx = 0.0; pady = 0.1;

[file,path] = uigetfile('C:\Users\mohseniankoochaksa.s\CCRC Dropbox\Mahsa Karamzadeh\Processed PC MRI\processed_chiari_data');
VENC = str2num(file(end-4));
file_path = convertCharsToStrings(path) + convertCharsToStrings(file);
savedir = ('C:\Users\mohseniankoochaksa.s\CCRC Dropbox\Mahsa Karamzadeh\Harvard_PCMRI\processed_chiari_data\Results');


[~, lastFolder] = fileparts(fileparts(path));
Patient_ID = lastFolder;
Patient_ID = strrep(Patient_ID, '_', ' ');

imageData = load(file_path()); 
% divide by 10 to convert from mm/s to cm/s
vel = imageData.processed_data/10;
sz_r = size(vel,1);        
sz_c = size(vel,2);
%pixel size - need to check
vsz_x = .6;
vsz_y = .6;
vsz_z = 3;
psz = [vsz_x, vsz_y vsz_z];
%need to check RR
RR=60/72;
NoPh = size(vel,3);
del_t = RR/NoPh;                    % time between two consecutive phases [s]
t_phases = 1:NoPh;                  % timing of the cardiac cycle [phases]
t = del_t*t_phases;                 % timing of the cardiac cycle [s]

%finding the min and max values and their indexes
[minValue,In_min] = min(vel,[],'all');
[maxValue,In_max] = max(vel,[],'all');
Data_size = size(vel);
[min_row, min_column, min_page]= ind2sub(Data_size,In_min);
[max_row, max_column, max_page]= ind2sub(Data_size,In_max);

mean_vel_each_pixel = mean(vel,3);       % temporal avg of velocity [cm/s]
figure(20)
imagesc(mean_vel_each_pixel);
colorbar;
% mean_vel_all = mean(vel,3);
% corrected_vel =vel(:,:,:)-mean_vel_all(:,:);


yn = questdlg('Would you like to see the phase movie', 'Movie?', {'Yes', 'No'});

fig4 = figure(4);
set(fig4,'Units','normalized','Position',[fig1x fig1y+height+pady width height])

% read the VENC from file name
    win = VENC;

if strcmp(yn, 'Yes')
    fig4 = figure(4);

    for k = 1: NoPh
        imshow(vel(:,:,k),[-win win]);
%         imshow(vel(:,:,k), [-win win]);
%         set(fig4,'Units','normalized','Position',[fig1x fig1y+height+pady width height])
        title(['Phase = ', num2str(k)], 'FontWeight','bold')
        axis off;
        pause(0.3)
    end
end

if strcmp(yn, 'No') || strcmp(yn, 'Yes')
    prompt = {'Enter phase # to plot phase image'};
    dlgtitle = 'Phase Input';
    dims = [1 40];
    definput = {'17','hsv'};
    UserInput = inputdlg(prompt, dlgtitle, dims, definput);
    ik = round(str2double(cell2mat(UserInput)));
elseif strcmp(yn, 'Cancel')
    ik = 17;
end

fig4 = figure(4);
imshow(vel(:,:,ik), [-win win]);


imcontrast
title('Click 2 pts: top left and bottom right');
roi = drawpoint;
Pos = roi.Position;
yi = Pos(2);            xi = Pos(1);
tlx = round(xi);        tly = round(yi);

roi = drawpoint;
Pos = roi.Position;
yi = Pos(2);            xi = Pos(1);
brx = round(xi);        bry = round(yi);

clear xi yi
close Figure 4

% Initialize a flag to control whether to draw the polygon
drawPolygon = true;
while drawPolygon==true

% s is three layers of the velocity scaled to be [0,1] to make it RGB later

s(:,:,1) = (vel(:,:,ik)+VENC)/(2*VENC);
s(:,:,2) = (vel(:,:,ik)+VENC)/(2*VENC);
s(:,:,3) = (vel(:,:,ik)+VENC)/(2*VENC);

fig4 = figure(4);
fig4.Units = "normalized";
fig4.Position = [fig1x fig1y+height+pady-0.05 width height];
imagesc(s(tly:bry, tlx:brx, :));
title(['Phase = ', num2str(ik)]);

mask = zeros(sz_r, sz_c);               % mask to be used in all calculation
mask_unwrapped = zeros(sz_r, sz_c);     % mask to know how many unwrapped pixels
vel_unwrapped = vel;                    % unwrapped velocity set equal to orig. vel. initially

for q = 1:NoPh
  roi_vel_avg(q) = sum(mask.*vel_unwrapped(:,:,q),'all')/...
  nnz(mask.*vel_unwrapped(:,:,q));           % avg vel [cm/s]
end



%draw polygone on the cropped image and make a mask
roi_poly = drawpolygon;
pos_poly = roi_poly.Position;
xadjusted_poly = pos_poly(:,1) + tlx; % to have the same coordinate as mask
yadjusted_poly = pos_poly(:,2) + tly;

mask_poly = round(pos_poly,0) + [tlx - 1, tly - 1];

mask = poly2mask(mask_poly(:,1), mask_poly(:,2),sz_r, sz_c);





%calculate the mean-vel of mask to subtract it from the whole to correct
%the offset. Offset of different location is different. That's why we
%select areas seperately to first subtract the offset and then correcting
%the velocity
for q = 1:NoPh
  roi_vel_avg(q) = sum(mask.*vel_unwrapped(:,:,q),'all')/...
  nnz(mask.*vel_unwrapped(:,:,q));           % avg vel [cm/s]
end
mean_vel = mean(roi_vel_avg);       % temporal avg of velocity [cm/s]


%subtracting the mean_vel from the velocity and then correct the velocity.
%This one is for correcting the error of machin
corrected_vel = vel(:,:,:)-mean_vel;

%To force the final displacement to be zero 
mean_vel_all = mean(corrected_vel,3);
corrected_vel =corrected_vel(:,:,:)-mean_vel_all(:,:);

check_vel = mean_vel_all > 0.05*VENC;
[row_check,col_check] = find(check_vel);

w(:,:,1) = (vel(:,:,17)+VENC)/(2*VENC);
w(:,:,2) = (vel(:,:,17)+VENC)/(2*VENC);
w(:,:,3) = (vel(:,:,17)+VENC)/(2*VENC);
for iii=1:size(row_check)
w(row_check(iii), col_check(iii), :) = Red;
end

figure(91);
imshow(w)

% %exclude pixels with mean>0.05*venc from the mask %remove red voxels
% for i = 1:length(row_check)
%     mask(row_check(i), col_check(i)) = 0;
% end

fig1 = figure(1);
fig1.Units = "normalized";
fig1.Position = [fig1x fig1y width height];
imagesc(s(tly:bry, tlx:brx, :))     % (grey scale cropped image)
title(['Phase = ', num2str(ik),',   TL = [' num2str(tly),', ',...
     num2str(tlx), '],   BR = [', num2str(bry), ', ',num2str(brx),']'])

for i = 1:sz_r
    for j = 1:sz_c
        % Check if the condition i & j > 0 is met
        if mask(i, j)==1
            s(i, j, :) = Green;
        end
    end
end


%copy lines from rest of the code and putting inside the loop
%defining circle for finding min and max displacement
mmcircarea = 30;
%mmcircrad = sqrt(mmcircarea/pi); % length of radius in mm
mmcircrad =2;
numit = 50; % number of points on the circle
theta = 0:2*pi/numit:2*pi;
voxelsize = 1; %1 voxel = 1mm  ***check these lines again after finding the pixel size
cirpix = mmcircrad/voxelsize;

%finding x and y cordinates of ones in mask
%row is y and column is x
[rones,cones] = find(mask); %coordinates are based on actual image 




 for C = 1:nnz(mask)
     vel_mask(C,:) = corrected_vel(rones(C), cones(C),:) ;

        %defining 50 points on the edge of the circle
    for ii=1:50
   Xcirc(ii,C) = cirpix*cos(theta(ii)) + cones(C);
   Ycirc(ii,C) = cirpix*sin(theta(ii)) + rones(C); 

 end
 end
 d=0;
in=zeros(numit,nnz(mask)); % checking if all those 50 points on each circles are inside the mask
 for b = 1:nnz(mask)
      in(:,b) = inpolygon(Xcirc(:,b),Ycirc(:,b),xadjusted_poly, yadjusted_poly);
   

      if nnz(in(:,b)) > numit/1.5 %if only numit/2 is inside the polygone that is accepteble
          d=d+1;
        cin(d) = cones(b);  %rones is row which is y
        rin(d) = rones(b);
        xcircle_small(:,d) = Xcirc(:,b);
        ycircle_small(:,d) = Ycirc(:,b);
        vel_mask2(d,:) = vel_mask(b,:);
      end
 end
 figure(6)

 plot(xadjusted_poly,yadjusted_poly)
 axis equal;
 set(gca, 'YDir','reverse')

 hold on
 plot(cin,rin,'.')
 axis equal;
 set(gca, 'YDir','reverse')
 hold off

%mask_circle = zeros(sz_r, sz_c); % mask to be used for circle calculation (min and max vel)

for v = 1:size(cin,2)
    %for each small circle, check which pixels are inside those small circles
    in_circle = inpolygon(cones, rones, xcircle_small(:,v),ycircle_small(:,v));
    [in_circle_row,~] = find(in_circle);
    %mask_circle(rones(in_circle_row), cones(in_circle_row)) = 1;


    for i=1:NoPh
        roi_vel_avg_circle(i,v) = mean(vel_mask(in_circle_row,i));
    end


mean_vel_circle(v) = mean(roi_vel_avg_circle(:,v));       % temporal avg of velocity [cm/s]

roi_disp_circle_temp(:,v) = 10000*(cumtrapz(del_t, roi_vel_avg_circle(:,v)));
displacement_circ(v) = max(roi_disp_circle_temp(:,v)) - min(roi_disp_circle_temp(:,v));

end


[max_disp,v_row_max] = max(displacement_circ);
[min_disp,v_row_min] = min(displacement_circ);
figure(8);
axis equal
imagesc(s(tly:bry, tlx:brx, :));
hold on
plot((xcircle_small(:,v_row_max)-tlx),(ycircle_small(:,v_row_max)-tly),'r')
plot((xcircle_small(:,v_row_min)-tlx),(ycircle_small(:,v_row_min)-tly),'b')
axis equal

dismat=mask*0;
for ppp = 1:size(cin,2)
qx=cin(ppp);
qy=rin(ppp);
qdisplacement=displacement_circ(ppp);
dismat(qy,qx)=qdisplacement; %must flip x and y for images. row, column for images
end



figure(9)
% it is for making CSF pixels black
% alpha= dismat<=750;
% imagesc(dismat,'AlphaData',alpha);
% set(gca,'color',[0 0 0]) % or any color you want
imagesc(dismat);

caxis([0, 900]);
colorbar;
axis equal

figure(10)

% Define the subplot layout
rows = 5;
cols = 6;
clims = [-1 1];
strainmat=mask*0;
for t = 1:30
    subplot(rows, cols, t); 
    for v = 1:size(cin,2)-1

        strain(v) = (roi_disp_circle_temp(t,v)-roi_disp_circle_temp(t,v+1))*0.01;
        qxx=cin(v);
        qyy=rin(v);

        qstrain=strain(v);
        strainmat(qyy,qxx)=qstrain; %must flip x and y for images. row, column for images

    end
imagesc(strainmat(min(rin):max(rin),min(cin):max(cin)),clims);
colorbar;
% axis equal
title(['phase=', num2str(t)], 'FontSize', 10);

end




    % CALCULATING PARAMETERS OF INTEREST: VELOCITY, Displacement, Q(T), AND SYS/DIAS VOLUMES
    if (i+tlx ~= tlx) && (j+tly ~= tly) && (i+tlx ~= brx) && (j+tly ~= bry)   % to exclude top left and bottom down pixel from calculation
        for q = 1:NoPh
            %roi_vel_avg(q) = sum(mask.*vel_unwrapped(:,:,q),'all')/...
                %nnz(mask.*vel_unwrapped(:,:,q));           % avg vel [cm/s]
            roi_corrected_vel_avg(q) = sum(mask.*corrected_vel(:,:,q),'all')/...
                nnz(mask.*vel_unwrapped(:,:,q));           % avg vel [cm/s]
            roi_vel_max(q) = max((mask.*vel_unwrapped(:,:,q)),[],'all'); % max vel [cm/s]
            roi_vel_min(q) = min((mask.*vel_unwrapped(:,:,q)),[],'all'); % max vel [cm/s]
        end
            %mean_vel = mean(roi_vel_avg);       % temporal avg of velocity [cm/s]
            roi_disp_min = 10000*min(cumtrapz(del_t, (roi_corrected_vel_avg)));
            roi_disp_max = 10000*max(cumtrapz(del_t, (roi_corrected_vel_avg)));
            vel_max = max(roi_corrected_vel_avg)
            vel_min = min(roi_corrected_vel_avg)

%figure for seeing the velocity of tissue
% Create a single figure for all velocities
figure(5);

% Define the subplot layout
rows = 5;
cols = 6;

for t = 1:30
    subplot(rows, cols, t); 
    clims = [-0.2 0.2];
    imagesc(corrected_vel(tly:bry, tlx:brx, t), clims);
    colorbar;
    title(['phase=', num2str(t)], 'FontSize', 10);

end


        %Area = round(nnz(mask)*psz(1)*psz(2)/100,2);   % Area [cm^2]
        Area = 1.;
        Q = Area*(roi_vel_avg - mean_vel);             % adjusted Q(t) [ml/s]
        Q_avg = round(mean(Q),3);                      % Q(t)_avg [ml/s]

        [~,index_dia] = find(Q >0);
        t_diastolic = del_t*nnz(index_dia);
        vol_dia = round(t_diastolic*mean(Q(index_dia)),3);      % approximate volume only

        [~,index_sys] = find(Q <0);
        t_systolic = del_t*nnz(index_sys);
        vol_sys = round(t_systolic*mean(Q(index_sys)),3);       % approximate volume only
    end

    % PLOTTING VELOCITY AND Q(T), AND UPDATING THE PLOT FOR EACH ITERATION
    fig3 = figure(3);
    fig3.Units = "normalized";
    fig3.Position = [fig1x+2*width+2*padx fig1y width height];
    subplot(2,1,1)
    % plot(t_phases,roi_vel_avg, '.-', 'Color', Blue, 'MarkerSize',10)
    % hold on
    % hold on
        plot(t_phases,roi_corrected_vel_avg, 'g', 'MarkerSize',14,'LineWidth',3)
    yline(0, 'k--')
ax = gca;
ax.FontSize = 16; 
ax.FontWeight="bold";

    hold off
    xlabel('t [phases]','FontWeight','bold');
    ylabel('V_{avg}(t) [cm/s]','FontWeight','bold');
    % text(0.05,0.2, ['mean-vel = ',num2str(mean_vel),' cm/s'],...
    %      'FontWeight','bold', 'Units','normalized')
    title("Patient ID = " + Patient_ID);

    subplot(2,1,2)
    plot(t_phases,Q, '.-', 'Color', Red, 'MarkerSize',30)
    plot(t_phases,10000*cumtrapz(del_t, roi_corrected_vel_avg), '.-', 'Color', Blue, 'MarkerSize',30,'LineWidth',3)
    ax = gca;
ax.FontSize = 16; 
ax.FontWeight="bold";

    hold on
    yline(0, 'k--')
    hold off
    xlabel('t [phases]','FontWeight','bold');
    ylabel('Displacement(t) [micron]','FontWeight','bold');
        % text(5, 2, ['max-disp = ', char(10), 'min-disp = ']);
        % text(11, 2, {num2str(roi_disp_max), num2str(roi_disp_min)});
%str = {['Max-disp = ' num2str(roi_disp_max)], ['Min-disp = ' num2str(roi_disp_min)], ['Displacement= ' num2str((roi_disp_max)-(roi_disp_min))]};
%text([2 10 20], [20 20 20], str);
str = { ['Displacement= ' num2str((roi_disp_max)-(roi_disp_min), '%.0f')]};
% text([2 10 20], [20 20 20], str);
% Set the position for the text
textHandle = annotation('textbox', [0.15, 0.1, 0.1, 0.1], 'String', str,'FontSize',25, 'FitBoxToText', 'on', 'LineStyle', 'none',FontWeight='bold');
yesgreen = 1;
count = 0;
countT = 1000;
drawPolygon=false;

while count <= countT
    count = count + 1;

    % COLORING THE PIXEL (BLUE = JUST SELECTED, GREEN = PREVIOUS SELECTION)
    fig1 = figure(1);
    fig1.Units = "normalized";
    fig1.Position = [fig1x fig1y width height];
    imagesc(s(tly:bry, tlx:brx, :))     % (grey scale cropped image)
    title(['Phase = ', num2str(ik),',   TL = [' num2str(tly),', ',...
        num2str(tlx), '],   BR = [', num2str(bry), ', ',num2str(brx),']'])
    subtitle(['Terminate = Top left pixel, Clear = Bottom right pixel'])
    
    roi = drawpoint;                    % (pixel selected)
    Pos = roi.Position;
    % in drawpoint it return the x and y cordinate but in an image x and y
    % are reversed
    yi = Pos(2);
    xi = Pos(1);
    if count > 1
        if yesgreen > 0.5
            s(ryi,rxi,:) = Green;   % (pixel after 1st round should be green except current)
        end
    end
    rxi = round(xi) + tlx -1
    ryi = round(yi) + tly -1

    s(ryi,rxi,:) = Blue;            % (current pixel appears in blue)
    imagesc(s(tly:bry, tlx:brx, :))

    fig4 = figure(4);
    title(['Phase = ', num2str(ik)]);


    % CHECKING THE LOGIC FOR CLICKING AND UNCLICKING THE PIXEL
    velo = (vel(ryi,rxi,ik)+VENC)/(2*VENC);

    if count == 1       % (just started)
        s(ryi,rxi,:) = Blue;
        mask(ryi,rxi) = 1;
    else                % (secound or more iteration)
        [a, b] = find(mask);
        IndexM = [a, b];
        logicM = (IndexM(:,1:2) == [ryi,rxi]);
        logicfind = logicM(:,1).*logicM(:,2);
        if isempty(find(logicfind == 1)) % (there is no match, keep it blue)
            s(ryi,rxi,:) = Blue;
            yesgreen = 1;
            mask(ryi,rxi) = 1;
        else                             % (there is a match, ungreen it = original)
            s(ryi,rxi,1) = velo;
            s(ryi,rxi,2) = velo;
            s(ryi,rxi,3) = velo;

            mask(ryi,rxi) = 0;
            yesgreen = 0;
        end
    end

    % CALCULATING PARAMETERS OF INTEREST: VELOCITY, Q(T), AND SYS/DIAS VOLUMES
    if (rxi ~= tlx) && (ryi ~= tly)&& (i+tlx ~= brx) && (j+tly ~= bry)    % to exclude top left pixel from calculation
        for i = 1:NoPh
            roi_vel_avg(i) = sum(mask.*vel_unwrapped(:,:,i),'all')/...
                nnz(mask.*vel_unwrapped(:,:,i));           % avg vel [cm/s]
            roi_vel_max(i) = max((mask.*vel_unwrapped(:,:,i)),[],'all'); % max vel [cm/s]
            roi_vel_min(i) = min((mask.*vel_unwrapped(:,:,i)),[],'all'); % max vel [cm/s]
        end
            mean_vel = mean(roi_vel_avg);       % temporal avg of velocity [cm/s]
            roi_disp_min = 10000*min(cumtrapz(del_t, (roi_corrected_vel_avg)));
            roi_disp_max = 10000*max(cumtrapz(del_t, (roi_corrected_vel_avg)));
        

        %Area = round(nnz(mask)*psz(1)*psz(2)/100,2);   % Area [cm^2]
        Area = 1.;
        Q = Area*(roi_vel_avg - mean_vel);             % adjusted Q(t) [ml/s]
        Q_avg = round(mean(Q),3);                      % Q(t)_avg [ml/s]

        [~,index_dia] = find(Q >0);
        t_diastolic = del_t*nnz(index_dia);
        vol_dia = round(t_diastolic*mean(Q(index_dia)),3);      % approximate volume only

        [~,index_sys] = find(Q <0);
        t_systolic = del_t*nnz(index_sys);
        vol_sys = round(t_systolic*mean(Q(index_sys)),3);       % approximate volume only
    end
    % PLOTTING THE SIGNAL OF THE SELECTED PIXEL VS TIME
    fig2 = figure(2);
    fig2.Units = "normalized";
    fig2.Position = [fig1x+width+padx fig1y width height];
        subplot(2,1,1)
    plot(t_phases, squeeze(vel_unwrapped(ryi,rxi,:)), 'bo-')
    hold on
    plot(t_phases, squeeze(vel_unwrapped(ryi,rxi,:))-mean_vel_all(ryi,rxi), 'g')
    plot(t_phases, squeeze(corrected_vel(ryi,rxi,:)), 'r')


    grid on;
    ylabel('Velocity [cm/s]','FontWeight','bold');
    xlabel('t [phases]','FontWeight','bold');
    title(['Phase = ',num2str(ik),', r = ' num2str(ryi), ' c = ', num2str(rxi)])
    hold off
        subplot(2,1,2)
     plot(t_phases,10000*cumtrapz(del_t, squeeze(corrected_vel(ryi,rxi,:))), '.-', 'Color', Blue, 'MarkerSize',10)
    hold on
    yline(0, 'k--')
    hold off
    xlabel('t [phases]','FontWeight','bold');
    ylabel('Displacement(t) [micron]','FontWeight','bold');

    % PLOTTING VELOCITY AND Q(T), AND UPDATING THE PLOT FOR EACH ITERATION
    fig3 = figure(3);
    fig3.Units = "normalized";
    fig3.Position = [fig1x+2*width+2*padx fig1y width height];
    subplot(2,1,1)
    % plot(t_phases,roi_vel_avg, '.-', 'Color', Blue, 'MarkerSize',10)
    % hold on
    plot(t_phases,roi_corrected_vel_avg, 'g', 'MarkerSize',10)
    yline(0, 'k--')
ax = gca;
ax.FontSize = 16; 
    hold off
    xlabel('t [phases]','FontWeight','bold');
    ylabel('V_{avg}(t) [cm/s]','FontWeight','bold');


    subplot(2,1,2)
    plot(t_phases,Q, '.-', 'Color', Blue, 'MarkerSize',10)
    plot(t_phases,10000*cumtrapz(del_t, roi_corrected_vel_avg), '.-', 'Color', Blue, 'MarkerSize',10)
    ax = gca;
ax.FontSize = 16; 
ax.FontWeight="bold";
    hold on
    yline(0, 'k--')
    hold off
        xlabel('t [phases]','FontWeight','bold');
        ylabel('Displacement(t) [micron]','FontWeight','bold');
        text(0.05,0.2, ['max-disp = ',num2str(roi_disp_max),' micron, min-disp = ',num2str(roi_disp_min)],...
        'FontWeight','bold', 'Units','normalized')
str = {['Max-disp = ' num2str(roi_disp_max)], ['Min-disp = ' num2str(roi_disp_min)], ['Displacement= ' num2str((roi_disp_max)-(roi_disp_min))]};
text([2 10 20], [20 20 20], str);


    % PHASE UNWRAP PROCEDURE
    fig1 = figure(1);
    imagesc(s(tly:bry, tlx:brx, :))

    % right click the pixel to initiate phase unwrapp...
    if strcmp(gcf().SelectionType, 'alt') == 1
        disp(count)
        mask(ryi,rxi) = 1;
        mask_unwrapped(ryi,rxi) = 1;
        disp('Phase Angle (- to subtract the VENC, and + to add)');
        prompt = [{'Enter Phase Start'},{'Enter Phase End'}, {'Enter + or -'}];
        answer = inputdlg(prompt, 'Phase Unwrapp Info', [1 50]);

        if ~isempty(answer)
            j = str2double(answer{1});
            J = str2double(answer{2});

            if strcmp(answer{3}, '+') == 1
                while j <= J
                    vel_unwrapped(ryi,rxi,j) = vel(ryi,rxi,j) + 2*VENC;
                    s(ryi,rxi,1) = (vel_unwrapped(ryi,rxi,14)+VENC)/(2*VENC);
                    s(ryi,rxi,2) = (vel_unwrapped(ryi,rxi,14)+VENC)/(2*VENC);
                    s(ryi,rxi,3) = (vel_unwrapped(ryi,rxi,14)+VENC)/(2*VENC);

                    disp(j);
                    j = j + 1;
                end


                if max(vel_unwrapped(ryi,rxi,:)) > VENC
                    ylim([-VENC max(vel_unwrapped(ryi,rxi,:))]);
                elseif min(vel_unwrapped(ryi,rxi,:)) < -VENC
                    ylim([min(vel_unwrapped(ryi,rxi,:)) VENC]);
                else
                    ylim([-VENC +VENC]);
                end
                grid on;
                ylabel('Velocity [cm/s]','FontWeight','bold');
                xlabel('t [phases]','FontWeight','bold');
                title(['Phase = ',num2str(ik),', r = ' num2str(ryi), ' c = ', num2str(rxi)])

            else
                while j <= J
                    vel_unwrapped(ryi,rxi,j) = vel(ryi,rxi,j) - 2*VENC;
                    s(ryi,rxi,1) = (vel_unwrapped(ryi,rxi,14)+VENC)/(2*VENC);
                    s(ryi,rxi,2) = (vel_unwrapped(ryi,rxi,14)+VENC)/(2*VENC);
                    s(ryi,rxi,3) = (vel_unwrapped(ryi,rxi,14)+VENC)/(2*VENC);

                    disp(j);
                    j = j + 1;
                end


                if max(vel_unwrapped(ryi,rxi,:)) > VENC
                    ylim([-VENC max(vel_unwrapped(ryi,rxi,:))]);
                elseif min(vel_unwrapped(ryi,rxi,:)) < -VENC
                    ylim([min(vel_unwrapped(ryi,rxi,:)) VENC]);
                else
                    ylim([-VENC +VENC]);
                end
                grid on;
                ylabel('Velocity [cm/s]','FontWeight','bold');
                xlabel('t [phases]','FontWeight','bold');
                title(['Phase = ',num2str(ik),', r = ' num2str(ryi), ' c = ', num2str(rxi)])
            end
            disp('Phase Unwrapped Done!...');
            fig1 = figure(1);
            s(ryi,rxi,:) = Yellow;             % (fixed wrapped pixel appears in yellow)
            imagesc(s(tly:bry, tlx:brx, :))

            fig4 = figure(4);

        else
            disp('Phase Unwrapped Cancelled!...');
            fig1 = figure(1);
            s(ryi,rxi,:) = Red;                % (ignored wrapped pixel appears in red)
            imagesc(s(tly:bry, tlx:brx, :))

            fig4 = figure(4);
            title(['Phase = ', num2str(ik)]);
        end
    end




    % TERMINATION PROCEDURE
    % to terminate click the pixel in the top left corner
    if (rxi == tlx) && (ryi == tly)
        mask(ryi,rxi) = 0;
        disp(count)
        count = countT + 1;
        disp('Terminated!...');
    end
    %to re-do segmentation
    if (rxi == brx) && (ryi == bry)
        count = countT + 1;

       disp('Cleared! Re-do the segmentation')
       drawPolygon=true;
    end


end
end

figure(7);   
point = drawpoint;                    % (pixel selected)
Point = point.Position;
% in drawpoint it return the x and y cordinate but in an image x and y
% are reversed
yi = round(Point(2));
xi = round(Point(1));
while (xi ~= 1) && (yi ~= 1)
    fig2 = figure(2);
    fig2.Units = "normalized";
    fig2.Position = [fig1x+width+padx fig1y width height];
        subplot(2,1,1)
    plot(t_phases, squeeze(vel_unwrapped(yi,xi,:)), 'bo-')
    hold on
    plot(t_phases, squeeze(vel_unwrapped(yi,xi,:))-mean_vel_all(yi,xi), 'g')
    hold on
    plot(t_phases, squeeze(corrected_vel(yi,xi,:)), 'r')


    grid on;
    ylabel('Velocity [cm/s]','FontWeight','bold');
    xlabel('t [phases]','FontWeight','bold');
    title(['Phase = ',num2str(ik),', r = ' num2str(yi), ' c = ', num2str(xi)])
    hold off
        subplot(2,1,2)
     plot(t_phases,10000*cumtrapz(del_t, squeeze(corrected_vel(yi,xi,:))), '.-', 'Color', Blue, 'MarkerSize',10)
    hold on
    yline(0, 'k--')
    hold off
    xlabel('t [phases]','FontWeight','bold');
    ylabel('Displacement(t) [micron]','FontWeight','bold');
    figure(7);   
    point = drawpoint;                    % (pixel selected)
    Point = point.Position;
    % in drawpoint it return the x and y cordinate but in an image x and y
    % are reversed
    yi = round(Point(2));
    xi = round(Point(1));

    end    
    