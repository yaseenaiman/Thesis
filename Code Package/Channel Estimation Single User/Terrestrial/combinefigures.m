clear;
close all;
fig1 = openfig('C:\Users\Yaseen_Developer\Desktop\1bit.fig', 'reuse');
fig2 = openfig('C:\Users\Yaseen_Developer\Desktop\2bit.fig', 'reuse');
fig3 = openfig('C:\Users\Yaseen_Developer\Desktop\3bits.fig', 'reuse');
fig4 = openfig('C:\Users\Yaseen_Developer\Desktop\4bit.fig', 'reuse');


% Extract the axes from each figure
ax1 = findobj(fig1, 'type', 'axes');
ax2 = findobj(fig2, 'type', 'axes');
ax3 = findobj(fig3, 'type', 'axes');
ax4 = findobj(fig4, 'type', 'axes');

% Extract the line objects from each axes
line1 = findobj(ax1, 'type', 'line');
line2 = findobj(ax2, 'type', 'line');
line3 = findobj(ax3, 'type', 'line');
line4 = findobj(ax4, 'type', 'line');


% Extract data from each line object
x1 = get(line1, 'XData');
y1 = get(line1, 'YData');

x2 = get(line2, 'XData');
y2 = get(line2, 'YData');

x3 = get(line3, 'XData');
y3 = get(line3, 'YData');


x4 = get(line4, 'XData');
y4 = get(line4, 'YData');

% Create a new figure
combinedFig = figure;

% Plot the data into the new figure
hold on; % Hold on to add multiple plots to the same axes
plot(x1{3}, y1{3},'-+c','LineWidth',1.5);
plot(x1{2}, y1{2},'-^g','LineWidth',1.5);
plot(x1{1}, y1{1},'->','Color',[0.5 0 0.8],'LineWidth',1.5);

plot(x2{3}, y2{3},'-rs','LineWidth',1.5);
plot(x2{2}, y2{2},'-b*','LineWidth',1.5);
plot(x2{1}, y2{1},'-yo','LineWidth',1.5);

plot(x3{3}, y3{3},':p','LineWidth',1.5);
plot(x3{2}, y3{2},'--k^','LineWidth',1.5);
plot(x3{1}, y3{1},'-md','LineWidth',1.5);

plot(x4{3}, y4{3},'--bo','LineWidth',1.5);
plot(x4{2}, y4{2},'-.rs','LineWidth',1.5);
plot(x4{1}, y4{1},'-*','Color',[0 0 0.5],'LineWidth',1.5);


set(gcf,'color','w');
set(gca, 'FontName', 'Times New Roman');  % Specify the font name (e.g., Arial)
set(gca, 'FontSize', 12);       % Specify the font size (e.g., 12)
    
legend('LS, 1-bit ADC','OMP, 1-bit ADC','NOMP, 1-bit ADC', 'LS, 2-bit ADC','OMP, 2-bit ADC','NOMP, 2-bit ADC','LS, 3-bit ADC','OMP, 3-bit ADC','NOMP, 3-bit ADC','LS, 4-bit ADC','OMP, 4-bit ADC','NOMP, 4-bit ADC')
xlabel('SNR [dB]','FontSize',12);
ylabel('MMSE [dB]','FontSize',12);
grid;