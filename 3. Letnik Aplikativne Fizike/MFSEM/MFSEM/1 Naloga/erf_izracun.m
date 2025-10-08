function val = erf_izracun(z)

if z < 0
  val = -erf_izracun(-z)
endif

if z>0 && z<4
  val = erf_potencna(z)
endif

if z>=4
  val = 1-erfc_asymptotic(z)
endif

