function [newColors]=colors(varargin)
%   Script displays 31 non-standard colors,
%   the 8 standard Matlab colors,
%   and their respective RGB values.
%   The colors are displayed on white and black backgrounds to
%   show how well they stand out from the background.

%   R. Martinez
%   10.24.03

if isempty(varargin)
    iDisplay=0;
else
    iDisplay=varargin{1};
end
newColors(1:31,1) = {'Aquamarine'; ...
        'Coral'; ...
        'Dark Green'; ...
        'Dark Khaki'; ...
        'Dark Olive Green'; ...
        'Dark Salmon'; ...
        'Dark Sea Green'; ...
        'Deep Pink'; ...
        'Deep Sky Blue'; ...
        'Forest Green'; ...
        'Ghost White'; ...
        'Gold'; ...
        'Khaki'; ...
        'Lavender Blush'; ...
        'Lawn Green'; ...
        'Lemon Chiffon'; ...
        'Maroon'; ...
        'Medium Aquamarine'; ...
        'Midnight Blue'; ...
        'Mint Cream'; ...
        'Navy Blue'; ...
        'Old Lace'; ...
        'Plum'; ...
        'Saddle Brown'; ...
        'Salmon'; ...
        'Sandy Brown'; ...
        'Tan'; ...
        'Tomato'; ...
        'Turquoise'; ...
        'Violet'; ...
        'Wheat'};

newColors(1:31,2) = { [0.26 0.72 0.73]; ...
        [0.97 0.40 0.25]; ...
        [0.15 0.25 0.09]; ...
        [0.72 0.68 0.35]; ...
        [0.29 0.25 0.09]; ...
        [0.88 0.55 0.42]; ...
        [0.55 0.70 0.51]; ...
        [0.96 0.16 0.53]; ...
        [0.23 0.73 1.00]; ...
        [0.31 0.57 0.35]; ...
        [0.97 0.97 1.00]; ...
        [0.83 0.63 0.09]; ...
        [0.68 0.66 0.43]; ...
        [0.99 0.93 0.96]; ...
        [0.53 0.97 0.09]; ...
        [1.00 0.97 0.78]; ...
        [0.51 0.02 0.25]; ...
        [0.20 0.53 0.51]; ...    
        [0.08 0.11 0.33]; ...
        [0.96 1.00 0.98]; ...
        [0.08 0.02 0.40]; ...
        [0.99 0.95 0.89]; ...
        [0.73 0.23 0.56]; ...
        [0.49 0.19 0.09]; ...
        [0.88 0.55 0.42]; ...
        [0.93 0.60 0.30]; ...
        [0.85 0.69 0.47]; ...
        [0.97 0.33 0.19]; ...
        [0.26 0.78 0.86]; ...
        [0.55 0.22 0.79]; ...
        [0.95 0.85 0.66]};

standardColors(1:8,1) = {'yellow'; ...
        'magenta'; ...
        'cyan'; ...
        'red'; ...
        'green'; ...
        'blue'; ...
        'white'; ...
        'black'};

standardColors(1:8,2) = {[1 1 0]; ...
        [1 0 1]; ...
        [0 1 1]; ...
        [1 0 0]; ...
        [0 1 0]; ...
        [0 0 1]; ...
        [1 1 1]; ...
        [0 0 0]};

if iDisplay==1
    %   First, the new colors ...
    
    figure
    plot(0,0)
    axis([0 40 0 65])
    set(gca, 'Color', 'white')
    hold on
    %for i = 1:31
    
    %text(2,64-i*2,newColors{i,1}, 'Color', newColors{i,2}, 'FontWeight', 'bold')
    %text(10,64-i*2,num2str(newColors{i,2}), 'Color', newColors{i,2}, 'FontWeight', 'bold')
    %end
    
    %   Then the standard colors.
    
    for i = 1:8
        text(25,i*5,standardColors{i,1}, 'Color', standardColors{i,2}, 'FontWeight', 'bold')
        text(30,i*5,num2str(standardColors{i,2}), 'Color', standardColors{i,2}, 'FontWeight', 'bold')
    end
    
    hold off
    
    figure
    plot(0,0)
    axis([0 40 0 65])
    set(gca, 'Color', 'black')
    hold on
    for i = 1:31
        text(2,64-i*2,newColors{i,1}, 'Color' ,newColors{i,2}, 'FontWeight', 'bold')
        text(10,64-i*2,num2str(newColors{i,2}), 'Color', newColors{i,2}, 'FontWeight', 'bold')
    end
    
    %   Then the standard colors.
    
    for i = 1:8
        text(25,i*5,standardColors{i,1}, 'Color',standardColors{i,2}, 'FontWeight', 'bold')
        text(30,i*5,num2str(standardColors{i,2}), 'Color', standardColors{i,2}, 'FontWeight', 'bold')
    end
    
    hold off
end
return

