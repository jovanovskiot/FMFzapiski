%parametri
intervali = {[-1,1],[-2,2],[-3,3]};
mu_vrednosti = [0, 1, 2, 5, 7]; %vzorcenje
N_vrednosti = round(logspace(2, 7, 10)); %1e4 do 1e6
n_trapez_vrednosti = round(logspace(2, 7, 10)); %1e2 do 1e4

prave_vrednosti = zeros(size(intervali));
napake_trapez = cell(size(intervali));
casi_trapez = cell(size(intervali));
napake_mc = cell(size(intervali));
casi_mc = cell(size(intervali));
napake_vzor = cell(size(intervali));
casi_vzor = cell(size(intervali));

normpdf = @(x) (1/sqrt(2*pi)) * exp(-x.^2/2);
normcdf = @(x) 0.5*(1 + erf(x/sqrt(2)));

koncne_vrednosti_trapez = zeros(size(intervali));
koncne_vrednosti_mc = zeros(size(intervali));
koncne_vrednosti_vzor = cell(size(intervali));

%%
for k = 1:length(intervali)
    a = intervali{k}(1);
    b = intervali{k}(2);
    prave_vrednosti(k) = normcdf(b)-normcdf(a);

    %Trapezna formula
    napake_trapez{k} = zeros(size(n_trapez_vrednosti));
    casi_trapez{k} = zeros(size(n_trapez_vrednosti));
    f = @(x) normpdf(x);
    for i = 1:length(n_trapez_vrednosti)
        n = n_trapez_vrednosti(i);
        x = linspace(a, b, n);
        y = f(x);
        h = (b-a)/(n-1);
        tic;
        I = h/2 * (y(1) + 2*sum(y(2:end-1)) + y(end));
        casi_trapez{k}(i) = toc;
        napake_trapez{k}(i) = abs(I - prave_vrednosti(k));
        if i == length(n_trapez_vrednosti)
            koncne_vrednosti_trapez(k) = I;
        end
    end

    %Monte carlo
    napake_mc{k} = zeros(size(N_vrednosti));
    casi_mc{k} = zeros(size(N_vrednosti));
    for i=1:length(N_vrednosti)
        N = N_vrednosti(i);
        tic;
        x_enakomerno = a + (b-a)*rand(N, 1);
        I = (b-a) * mean(normpdf(x_enakomerno));
        casi_mc{k}(i) = toc;
        napake_mc{k}(i) = abs(I - prave_vrednosti(k));
        if i == length(N_vrednosti)
            koncne_vrednosti_mc(k) = I;
        end
    end
    
    %prednostno vzorcenje
    napake_vzor{k} = zeros(length(mu_vrednosti), length(N_vrednosti));
    casi_vzor{k} = zeros(length(mu_vrednosti), length(N_vrednosti));
    koncne_vrednosti_vzor{k} = zeros(length(mu_vrednosti));
    for m = 1:length(mu_vrednosti)
        mu = mu_vrednosti(m);
        for i = 1:length(N_vrednosti)
            N = N_vrednosti(i);
            tic;
            x_vzor = mu + randn(N, 1);
            v_intervalu = (x_vzor >= a) & (x_vzor <= b);
            f = normpdf(x_vzor);
            g = (1/sqrt(2*pi)) * exp(-(x_vzor - mu).^2 / 2);
            f_g = f ./ g;
            I = mean(v_intervalu .* f_g);
            casi_vzor{k}(m, i) = toc;
            napake_vzor{k}(m, i) = abs(I - prave_vrednosti(k));
            if i == length(N_vrednosti)
                koncne_vrednosti_vzor{k}(m) = I;
            end
        end
    end
end

%%
% Plotting results for one interval (e.g., [-1, 1])
k = 2; % Change index for other intervals
a = intervali{k}(1);
b = intervali{k}(2);
% Error vs N for Trapezoidal and MC
figure;
loglog(n_trapez_vrednosti, napake_trapez{k}, 'b-o', 'DisplayName', 'Trapezna');
hold on;
loglog(N_vrednosti, napake_mc{k}, 'r-s', 'DisplayName', 'Monte Carlo');
xlabel('Število točk');
ylabel('Absolutna napaka');
legend;
%title(['Error Comparison: [', num2str(a), ', ', num2str(b), ']']);
grid on;
% Time vs N for Trapezoidal and MC
figure;
loglog(n_trapez_vrednosti, casi_trapez{k}, 'b-o', 'DisplayName', 'Trapezna');
hold on;
loglog(N_vrednosti, casi_mc{k}, 'r-s', 'DisplayName', 'Monte Carlo');
xlabel('Število točk');
ylabel('Čas (s)');
legend;
%title(['Time Comparison: [', num2str(a), ', ', num2str(b), ']']);
grid on;
% Error vs mu for Importance Sampling (fixed N)
figure;
N_idx = 3; % Example: N = 1e5
semilogy(mu_vrednosti, napake_vzor{k}(:, N_idx), 'b-o');
xlabel('\mu za prednostno vzorčenje');
ylabel('Absolutna napaka');
%title(['Error vs \mu (N=1e5): [', num2str(a), ', ', num2str(b), ']']);
grid on;

%%
% Displaying the final integral values
disp('Končne vrednosti integralov:');
for k = 1:length(intervali)
    fprintf('Interval: [%.1f, %.1f]\n', intervali{k}(1), intervali{k}(2));
    fprintf('  Prava vrednost: %.8f\n', prave_vrednosti(k));
    fprintf('  Trapezna metoda (N = %d): %.8f\n', n_trapez_vrednosti(end), koncne_vrednosti_trapez(k));
    fprintf('  Monte Carlo (N = %d): %.8f\n', N_vrednosti(end), koncne_vrednosti_mc(k));
    fprintf('  Prednostno vzorčenje (N = %d):\n', N_vrednosti(end));
    for m = 1:length(mu_vrednosti)
        fprintf('    mu = %.1f: %.8f\n', mu_vrednosti(m), koncne_vrednosti_vzor{k}(m));
    end
end