%Make truncation figure
%Parameters
Y1 =0; Y2 =0; %location of point source
Y1move = 0; Y2move = .2; %location to move point source
X1 = 0; X2 = .43;

trunc = 20;
theta = linspace(0.5,20,100);

bounds = @(k) bounds_freq(k, Y1, Y2, Y1move, Y2move, X1, X2, trunc);

[rmerr, rmb, rderr, rdb]= bounds(theta);
[ierr, imb, iderr, idb] = bounds(1i*theta);
[cerr, cmb, cderr, cdb]= bounds((1+sqrt(2)*1i)*1/sqrt(3)*theta);
[nerr, nmb, nderr, ndb]= bounds((99/100-1i*sqrt(199)/100)*theta);


thickLines;
figure(1); clf;
semilogy(nan, nan,'Color', 'k')
hold on;
plot(nan, nan, 'Color', 'k', 'LineStyle','--')
plot(theta, rmerr,'Color', [0.4660 0.6740 0.1880])
plot(theta, rmb, 'Color', [0.4660 0.6740 0.1880], 'LineStyle','--')
plot(theta, ierr,'Color', 'r')
plot(theta, imb, 'Color', 'r', 'LineStyle','--')
plot(theta, cerr,'Color', 'b')
plot(theta, cmb, 'Color', 'b', 'LineStyle','--')
plot(theta, nerr,'Color', 'm')
plot(theta, nmb, 'Color', 'm', 'LineStyle','--')
%legend('Real k error', 'Real k bound', 'Imaginary k error', 'Imaginary k bound', 'Complex k error', 'Complex k bound')
xlabel('\theta')
xlim([min(theta),max(theta)])
ylim([min(ierr),max(rmb)+.0005])
legend('errors', 'bounds');
hold off;

figure(2); clf;
semilogy(nan, nan,'Color', 'k')
hold on;
plot(nan, nan, 'Color', 'k', 'LineStyle','--')
plot(theta, rderr,'Color', [0.4660 0.6740 0.1880])
plot(theta, rdb, 'Color', [0.4660 0.6740 0.1880], 'LineStyle','--')
plot(theta, iderr,'Color', 'r')
plot(theta, idb, 'Color', 'r', 'LineStyle','--')
plot(theta, cderr,'Color', 'b')
plot(theta, cdb, 'Color', 'b', 'LineStyle','--')
plot(theta, nderr,'Color', 'm')
plot(theta, ndb, 'Color', 'm', 'LineStyle','--')
legend('errors', 'bounds');
%legend('Real k error', 'Real k bound', 'Imaginary k error', 'Imaginary k bound', 'Complex k error', 'Complex k bound')
xlabel('\theta')
xlim([min(theta),max(theta)])
hold off;
