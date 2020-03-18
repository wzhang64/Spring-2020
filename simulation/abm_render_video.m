% Author: Dhananjay Bhaskar
% Last Modified: May 23, 2018

lastposfile = 'Pos_200000.dat';
lastvelfile = 'Velocity_200000.dat';

maxiter = 200000;

adh_conditions = ["9"];
pol_conditions = ["3"];

toggle_cell_cycle = "on";
contact_inhibition_threshold = 4;

boxsize = 10;

for adh_idx = 1 : length(adh_conditions)

    for pol_idx = 1 : length(pol_conditions)
        
        adh_cond = char(adh_conditions(adh_idx));
        pol_cond = char(pol_conditions(pol_idx));
        
        folder_name = strcat(adh_cond, '_adh_', pol_cond, '_pol');
        
        if ~exist(strcat(folder_name, filesep, lastposfile), 'file')
            continue;
        end    
        
        disp(folder_name)
        
        vfilename = strcat(adh_cond, '_adhesion_', pol_cond, '_polarity.avi');
        vidWriter = VideoWriter(vfilename);
        open(vidWriter);
        
        iter = 0;
        
        while iter <= maxiter
            
            clf;
            
            posfilename = strcat('Pos_', num2str(iter, '%06.f'), '.dat');
            velfilename = strcat('Velocity_', num2str(iter, '%06.f'), '.dat');
            nbdfilename = strcat('Neighbors_', num2str(iter, '%06.f'), '.dat'); 
            edgexfilename = strcat('EdgeX_', num2str(iter, '%06.f'), '.dat'); 
            edgeyfilename = strcat('EdgeY_', num2str(iter, '%06.f'), '.dat'); 
            
            disp(posfilename);
            
            pos = csvread(strcat(folder_name, filesep, posfilename));
            vel = csvread(strcat(folder_name, filesep, velfilename));
            nbd = csvread(strcat(folder_name, filesep, nbdfilename));
            
            edges_exist = 1;
            try
                edgex = csvread(strcat(folder_name, filesep, edgexfilename));
                edgey = csvread(strcat(folder_name, filesep, edgeyfilename));
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
            
            % Plot edges
            if edges_exist == 1
                for k = 1 : size(edgex, 1)

                    plot([edgex(k,1) edgex(k,2)], [edgey(k,1) edgey(k,2)], 'linewidth', 1, 'Color', [0 0 0] + 0.7);
                    hold on

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
                
                hold on
    
            end

            avg_vel = avg_vel/numcells;
            avg_num_neighbors = avg_num_neighbors/numcells;

            rectangle('Position', [-1*boxsize -1*boxsize 2*boxsize 2*boxsize], 'LineWidth', 2.5)
            hold off
            
            cond_string = strcat(string(adh_cond), " density,", {' '}, string(pol_cond), " polarity,", {' '});
            stat_string = strcat("avg. speed:", {' '}, string(avg_vel), " avg. # neighbors:", {' '}, string(avg_num_neighbors));
            
            title_string = strcat(cond_string, stat_string);
            %title(title_string);
            
            set(gca, 'xtick', []);
            set(gca, 'ytick', []);
            set(gca, 'XTickLabel', []);
            set(gca, 'YTickLabel', []);
            
            frame = getframe(gcf);
            writeVideo(vidWriter, frame);
            
            iter = iter + 100;
            
        end
        
        close(vidWriter);
        
        close all;
        
    end
    
end