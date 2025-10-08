%prva naloga

%parametri za elipso in elipsoid 
a = 2; %x
b = 1; %y
c = 3; %z

prava_ploscina = pi * a * b;
pravi_volumen = (4/3)*pi*a*b*c;


%parametri za monte carlo
N_vrednosti = round(logspace(2, 6, 20));
napaka_ploscina = zeros(size(N_vrednosti));
napaka_volumen = zeros(size(N_vrednosti));
cas_ploscina = zeros(size(N_vrednosti));
cas_volumen = zeros(size(N_vrednosti));

%integral za ploscino
for i = 1:length(N_vrednosti)
    N = N_vrednosti(i);
    x = -a + 2*a*rand(N,1);
    y = -b + 2*b*rand(N,1);
    tic;
    znotraj = (x.^2 / a^2 + y.^2 / b^2) <= 1;
    priblizek_ploscina = 4*a*b * sum(znotraj) / N;
    napaka_ploscina(i) = abs(priblizek_ploscina-prava_ploscina);
    cas_ploscina(i) = toc;
end

% integral za volumen
for i = 1:length(N_vrednosti)
    N = N_vrednosti(i);
    x = -a + 2*a*rand(N, 1);
    y = -b + 2*b*rand(N, 1);
    z = -c + 2*c*rand(N, 1);
    tic;
    znotraj = (x.^2 / a^2 + y.^2 / b^2 + z.^2 / c^2) <= 1;
    priblizek_volumen = 8*a*b*c * sum(znotraj) / N;
    napaka_volumen(i) = abs(priblizek_volumen - pravi_volumen);
    cas_volumen(i) = toc;
end

%enakomerna mreza
n_vrednosti = [10, 20, 50, 100, 200, 500, 1000];
n_mreza = n_vrednosti.^2;
napaka_ploscina_mreza = zeros(size(n_vrednosti));
cas_mreza = zeros(size(n_vrednosti));

for i = 1:length(n_vrednosti)
    n = n_vrednosti(i);
    x = linspace(-a, a, n);
    y = linspace(-b, b, n);
    [X, Y] = meshgrid(x, y);
    tic;
    znotraj = (X.^2 / a^2 + Y.^2 / b^2) <= 1;
    priblizek_ploscine_mreza = 4*a*b * sum(znotraj(:)) / n^2;
    napaka_ploscina_mreza(i) = abs(priblizek_ploscine_mreza - prava_ploscina);
    cas_mreza(i) = toc;
end

% MC for grid-comparable N values
napaka_mrezaN = zeros(size(n_mreza));
cas_mrezaN = zeros(size(n_mreza));
for i = 1:length(n_mreza)
    N = n_mreza(i);
    x = -a + 2*a*rand(N, 1);
    y = -b + 2*b*rand(N, 1);
    tic;
    znotraj = (x.^2 / a^2 + y.^2 / b^2) <= 1;
    priblizek_ploscinaN = 4*a*b * sum(znotraj) / N;
    napaka_mrezaN(i) = abs(priblizek_ploscinaN - prava_ploscina);
    cas_mrezaN(i) = toc;
end

figure;
loglog(N_vrednosti, napaka_ploscina, 'b-o', 'DisplayName', 'Monte Carlo ploščina');
hold on;
loglog(N_vrednosti, napaka_volumen,'r-o', 'DisplayName', 'Monte Carlo volumen');
loglog(n_mreza, napaka_ploscina_mreza, 'g-s', 'DisplayName', 'Enakomerna mreža');
loglog(n_mreza, napaka_mrezaN, 'm--d', 'DisplayName', 'Monte Carlo mreža');
xlabel('Število točk');
ylabel('Absolutna napaka');
legend('Location', 'northeast');
%title('Napaka in število točk');
grid on;

figure;
loglog(n_mreza, cas_mreza, 'g-s', 'DisplayName', 'Mreža');
hold on;
loglog(n_mreza, cas_mrezaN, 'm--d', 'DisplayName', 'MC Mreža');
xlabel('Število točk');
ylabel('Čas');
legend('Location', 'northeast');
grid on;