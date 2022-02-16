%Make wavenumber range figure
%Parameters
Y1 =0; Y2 =0; %location of point source
Y1move = 0; Y2move = .2; %location to move point source
X1 = 0; X2 = .43;

trunc = 20;
theta = linspace(0.5,20,100);
thickLines

figure(1); clf
hold on
plot(theta,zeros(size(theta)), 'LineWidth', 5, 'Color', [0.4660 0.6740 0.1880])
plot(zeros(size(theta)), theta, 'LineWidth', 5, 'Color', 'r')
plot(theta*1/sqrt(3), theta*sqrt(2)/sqrt(3), 'LineWidth', 5, 'Color', 'b')
plot(theta*99/100, -theta*sqrt(199)/100, 'LineWidth', 5, 'Color', 'm')
xlim([0,20])
ylim([-5,20])
xlabel('Re')
ylabel('Im')
axis square
hold off