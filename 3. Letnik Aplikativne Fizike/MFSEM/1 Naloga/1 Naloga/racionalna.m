#racionalna
function erf_rac = racionalna(z_vrednosti)
  format long;


  erf_rac = [];
  p = 0.3275911;
  c = 1.421413741;
  a = 0.254829592;
  d = -1.453152027;
  b = -0.284496736;
  e = 1.61405429;

  for z = z_vrednosti
    t = 1/(1+p*z);
    erf_rac_iter = 1 - (a*t + b*t^2 + c*t^3 + d*t^4 + e*t^5)*exp(-z^2);
    erf_rac = [erf_rac, erf_rac_iter];
  endfor

