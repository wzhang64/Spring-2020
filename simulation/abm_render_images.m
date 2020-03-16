% Author: Dhananjay Bhaskar
% Last Modified: July 21, 2018

posfile = 'Pos_200000.dat';
velfile = 'Velocity_200000.dat';
nbdfile = 'Neighbors_200000.dat';
edgexfile = 'EdgeX_200000.dat';
edgeyfile = 'EdgeY_200000.dat';

toggle_cell_cycle = "off";

% Cells with 4 or more neighbours belong to a cluster
contact_inhibition_threshold = 4;

addpath(strcat(pwd, filesep, 'export_fig'))

boxsize = 10;

adhesion_values = linspace(0.05, 0.25, 11);
polarity_values = linspace(0.005, 0.025, 11);

adh_idx = 2;
pol_idx = 9;

folder_name = strcat(num2str(adh_idx), '_adh_', num2str(pol_idx), '_pol');

if ~exist(strcat(folder_name, filesep, posfile), 'file')
    exit();
end    

disp(folder_name)

pos = csvread(strcat(folder_name, filesep, posfile));
vel = csvread(strcat(folder_name, filesep, velfile));
nbd = csvread(strcat(folder_name, filesep, nbdfile));

edges_exist = 1;
try
    edgex = csvread(strcat(folder_name, filesep, edgexfile));
    edgey = csvread(strcat(folder_name, filesep, edgeyfile));
catch ME
    edges_exist = 0;
end

numcells = size(pos, 2);
avg_vel = 0.0;
avg_num_neighbors = 0.0;

if size(vel, 2) ~= numcells
    disp('Error: Size mismatch!');
    exit();
end

figure
hold on;

% Plot edges
if edges_exist == 1
    for k = 1 : size(edgex, 1)
        plot([edgex(k,1) edgex(k,2)], [edgey(k,1) edgey(k,2)], 'linewidth', 1.5, 'Color', [0 0 0] + 0.5);
    end
end

% Plot cells
for j = 1 : numcells

    num_nbd_cells = nbd(j);

    if num_nbd_cells == 0
        plot(pos(j), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', [1 .6 .6]);
    else
        if toggle_cell_cycle == "off"
            plot(pos(j), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', [.6 .6 1]);
        else
            if num_nbd_cells < contact_inhibition_threshold
                plot(pos(j), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', [.6 .6 1]);
            else
                plot(pos(j), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', [.1 .6 .1], 'MarkerFaceColor', [.4 .8 .4]);
            end
        end
    end

    avg_vel = avg_vel + sqrt(real(vel(j))^2 + imag(vel(j))^2);
    avg_num_neighbors = avg_num_neighbors + num_nbd_cells;

end

avg_vel = avg_vel/numcells;
avg_num_neighbors = avg_num_neighbors/numcells;

%rectangle('Position', [-1*boxsize -1*boxsize 2*boxsize 2*boxsize], 'LineWidth', 3)
hold off;

axis([-1*boxsize boxsize -1*boxsize boxsize]);
set(gca, 'xtick', [])
set(gca, 'ytick', [])
ax = gca;
ax.Visible = 'off';

cond_string = strcat(num2str(adh_idx), " adhesion,", {' '}, num2str(pol_idx), " polarity,", {' '});
stat_string = strcat("avg. speed:", {' '}, string(avg_vel), " avg. # neighbors:", {' '}, string(avg_num_neighbors));

outfilename = strcat(num2str(adh_idx), "_adhesion_", num2str(pol_idx), "_polarity.eps");

title_string = strcat(cond_string, stat_string);
title(title_string);

tightfig;

export_fig outfilename -eps -transparent