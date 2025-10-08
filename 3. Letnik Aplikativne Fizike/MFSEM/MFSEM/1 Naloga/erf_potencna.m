function val = erf_potencna(z)

vsota = 0;
n=0;

while true
  clen = (-1)^n*(z^(2*n+1))/(factorial(n)*(2*n+1));
  vsota = vsota + clen;
  n = n+1;
  if abs(clen) <= eps
    break
  endif
endwhile
val = 2/sqrt(pi)*vsota;

