%integracija
format long
% funkcija 1: 1 - x^2 od -1 do 1
f1 = @(x) 1 - x.^2;
a1 = -1; b1 = 1;
exact1 = 4/3;

% trapezna
n1 = 100; h1 = (b1 - a1)/n1;
x1 = a1:h1:b1; y1 = f1(x1);
trap1 = (h1/2)*(y1(1) + 2*sum(y1(2:end-1)) + y1(end));

% simpson
n1_simp = n1 + mod(n1,2); h1_simp = (b1 - a1)/n1_simp;
x1_simp = a1:h1_simp:b1; y1_simp = f1(x1_simp);
simp1 = (h1_simp/3)*(y1_simp(1) + 4*sum(y1_simp(2:2:end-1)) + 2*sum(y1_simp(3:2:end-2)) + y1_simp(end));

% MATLAB integral
matlab1 = integral(f1, a1, b1);

% funkcija 2: x^3 e^{-x} od 0 do Inf
f2 = @(x) x.^3 .* exp(-x);
a2 = 0; b2 = Inf; T = 50; exact2 = 6;

% trapezna
n2 = 1000; h2 = T/n2; x2 = 0:h2:T; y2 = f2(x2);
trap2 = (h2/2)*(y2(1) + 2*sum(y2(2:end-1)) + y2(end));

% Simpson
n2_simp = n2 + mod(n2,2); h2_simp = T/n2_simp; x2_simp = 0:h2_simp:T;
y2_simp = f2(x2_simp);
simp2 = (h2_simp/3)*(y2_simp(1) + 4*sum(y2_simp(2:2:end-1)) + 2*sum(y2_simp(3:2:end-2)) + y2_simp(end));

% MATLAB integral
matlab2 = integral(f2, a2, b2);

% rezultati
disp('1 - x^2 from -1 to 1:');
disp(['Tocno: ', num2str(exact1), ', Trapezna: ', num2str(trap1), ...
    ', Simpson: ', num2str(simp1), ', MATLAB: ', num2str(matlab1)]);
disp(['Trapezoidal: ', num2str(trap1), ', Napaka: ', num2str(abs(trap1 - exact1))]);
disp(['Simpson: ', num2str(simp1), ', Napaka: ', num2str(abs(simp1 - exact1))]);
disp(['MATLAB: ', num2str(matlab1), ', Napaka: ', num2str(abs(matlab1 - exact1))]);

disp('x^3 e^{-x} from 0 to Inf:');
disp(['Tocno: ', num2str(exact2), ', Trapezna: ', num2str(trap2), ...
    ', Napaka: ', num2str(abs(trap2 - exact2))]);
disp(['Simpson: ', num2str(simp2), ', Napaka: ', num2str(abs(simp2 - exact2))]);
disp(['MATLAB: ', num2str(matlab2), ', Napaka: ', num2str(abs(matlab2 - exact2))]);