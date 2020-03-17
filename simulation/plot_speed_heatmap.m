% Author: Dhananjay Bhaskar
% Last Modified: Jun 20, 2018

max_velocity = -Inf;
min_velocity = Inf;

speed_mat = zeros(11, 11);

adhesion_values = linspace(0.05, 0.25, 11);
polarity_values = linspace(0.005, 0.025, 11);

xlabelvec = {};
ylabelvec = {};

for i = 1 : length(adhesion_values)
    xlabelvec{end+1} = num2str(adhesion_values(i));
end

for i = 1 : length(polarity_values)
    ylabelvec{end+1} = num2str(polarity_values(i));
end

addpath(strcat(pwd, filesep, 'cbrewer'))
CT = cbrewer('seq', 'Reds', 121, 'PCHIP');

for adh_idx = 1 : length(adhesion_values)

    for pol_idx = 1 : length(polarity_values)
        
        folder_name = strcat(num2str(adh_idx), '_adh_', num2str(pol_idx), '_pol');
        
        velfile = 'Velocity_200000.dat';
        
        if ~exist(strcat(folder_name, filesep, velfile), 'file')
            continue;
        end    
        
        disp(folder_name)
        
        vel = csvread(strcat(folder_name, filesep, velfile));
        numcells = size(vel, 2);
            
        frame_avg_speed = 0;
        speeds = [];

        for j = 1 : numcells

            frame_avg_speed = frame_avg_speed + sqrt(real(vel(j))^2 + imag(vel(j))^2);
            speeds = [speeds sqrt(real(vel(j))^2 + imag(vel(j))^2)];
                
        end
            
        frame_avg_speed = frame_avg_speed/numcells;
        average_speed = mean(speeds);
        
        disp(strcat('Adhesion:', num2str(adh_idx), ' Polarity:', num2str(pol_idx), ' Speed:', num2str(average_speed)));
            
        if frame_avg_speed < min_velocity
            min_velocity = frame_avg_speed;
        end

        if frame_avg_speed > max_velocity
            max_velocity = frame_avg_speed;
        end
        
        speed_mat(adh_idx, pol_idx) = average_speed;
        
    end
    
end

disp(min_velocity)
disp(max_velocity)

h = heatmap(flipud(speed_mat'), 'colormap', CT, 'GridVisible', 'off', 'CellLabelColor', 'none');

h.XLabel = 'Adhesion force';
h.YLabel = 'Random cell polarity force';
h.YDisplayLabels = flip(ylabelvec);
h.XDisplayLabels = xlabelvec;

set(h, 'FontSize', 11, 'FontName', 'Myriad Pro Regular')

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 6.25 6.25]);

print('speed_heatmap.eps', '-depsc')