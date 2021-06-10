%% Specify the colormap

% Option 1: define your own colours
cmap = [1 0 0;      % red
        1 0.647 0;  % orange
        0 0 1;      % blue
        0 1 1;];    % cyan
    
% Option 2: use brewermap
% https://au.mathworks.com/matlabcentral/answers/432543-use-of-colours-with-brewermap

%cmap = brewermap(4, 'blues'); % generates 4 shades of blue


%% Plot the lines
set(0, 'DefaultAxesColorOrder', cmap(:,:)) 

figure; hold on;
plot(1:5, 6:10, 'Linewidth',3);
plot(1:5, 3:7, 'Linewidth',3);
plot(1:5, 2:6, 'Linewidth',3);
plot(1:5, 1:5, 'Linewidth',3);
hold off;


%% Another example
% https://au.mathworks.com/matlabcentral/answers/272874-how-to-display-points-from-very-light-red-to-dark-red

figure; scatter(1:100,1:100,50,1:100,'filled','markeredgecolor','k')
colormap(brewermap(2,'reds'))


%% Refs about colours:
https://au.mathworks.com/help/matlab/ref/linespec.html
https://au.mathworks.com/help/matlab/ref/colorspec.html
https://au.mathworks.com/matlabcentral/answers/151011-how-to-plot-a-line-of-a-certian-color
