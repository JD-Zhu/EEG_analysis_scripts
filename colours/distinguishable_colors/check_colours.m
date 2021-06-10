% call distinguishable_colors() to generate N colours,
% and plot them to see what they look like

N=9; 
colors = distinguishable_colors(N);

X = linspace(0,pi*3,1000); 
Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N); 
C = colors;

axes('NextPlot','replacechildren', 'ColorOrder',C); 
plot(X,Y,'linewidth',5) 
ylim([-1.1 1.1]);
