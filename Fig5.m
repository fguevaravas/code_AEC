%Make space figure
%Parameters
Y1 =0; Y2 =0; %location of point source
Y1move = 0; Y2move = .2; %location to move point source
X1 = 0;
X2s = linspace(0.4001,.6,3000);

trunc = 20;
bounds = @(k) bounds_space(k, Y1, Y2, Y1move, Y2move, X1, X2s, trunc);

[rmerr, rmb, rderr, rdb,dists]= bounds(1);
[ierr, imb, iderr, idb] = bounds(1i);
[cerr, cmb, cderr, cdb]= bounds((1+sqrt(2)*1i)*1/sqrt(3));
[nerr, nmb, nderr, ndb]= bounds(99/100-sqrt(199)/100*1i);

thickLines;
figure(1); clf;
semilogy(nan, nan,'Color', 'k')
hold on;
plot(nan, nan, 'Color', 'k', 'LineStyle','--')
plot(dists, rmerr,'Color', [0.4660 0.6740 0.1880])
plot(dists, rmb, 'Color', [0.4660 0.6740 0.1880], 'LineStyle','--')
plot(dists, ierr,'Color', 'r')
plot(dists, imb, 'Color', 'r', 'LineStyle','--')
plot(dists, cerr,'Color', 'b')
plot(dists, cmb, 'Color', 'b', 'LineStyle','--')
plot(dists, nerr,'Color', 'm')
plot(dists, nmb, 'Color', 'm', 'LineStyle','--')
legend('errors','bounds','Location','northwest')
xlabel('|y-x_j|/|x-x_j|')
hold off;

figure(2); clf;
semilogy(nan, nan,'Color', 'k')
hold on;
plot(nan, nan, 'Color', 'k', 'LineStyle','--')
plot(dists, rderr,'Color', [0.4660 0.6740 0.1880])
plot(dists, rdb, 'Color', [0.4660 0.6740 0.1880], 'LineStyle','--')
plot(dists, iderr,'Color', 'r')
plot(dists, idb, 'Color', 'r', 'LineStyle','--')
plot(dists, cderr,'Color', 'b')
plot(dists, cdb, 'Color', 'b', 'LineStyle','--')
plot(dists, nderr,'Color', 'm')
plot(dists, ndb, 'Color', 'm', 'LineStyle','--')
legend('errors','bounds','Location','northwest')
xlabel('|y-x_j|/|x-x_j|')
hold off;