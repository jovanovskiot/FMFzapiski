#naloga 1
format long;
# potencna vrsta
# uporabna za majhne z

#z_vrednosti = linspace(0, 3, 100); #za zacetek
z_vrednosti = [3:0.5:8]
erf_values = [];
erf_vgrajeni = erf(z_vrednosti);

#zanka
for z = z_vrednosti
  erf_approx = erf_potencna(z);
  erf_values = [erf_values, erf_approx];
endfor

# asimptotska vrsta
erf_values_as = []
erfc_values = []
for z = z_vrednosti
  erfc_values = [erfc_values, erfc_asymptotic(z)];
endfor

erf_values_as = 1 - erfc_values;

#racionalna
erf_rac = racionalna(z_vrednosti);

# primerjava relativnih razlik
comp_pot = log10(abs(erf_values-erf_vgrajeni)./erf_vgrajeni);
comp_as = log10(abs(erf_values_as-erf_vgrajeni)./erf_vgrajeni);
comp_rac = log10(abs(erf_rac-erf_vgrajeni)./erf_vgrajeni);;

#konvergenca asimptotske vrste
plot(z_vrednosti, comp_pot, 'r-')#, ... # Potenčna vrsta - rdeča polna
#     z_vrednosti, erf_values_as, 'g--', ...      # Asimptotska vrsta - zelena prekinjena
#     z_vrednosti, erf_rac, 'm-.', ... # Racionalna aproksimacija - magenta črtkano-pikčasta
#     z_vrednosti, erf_vgrajeni, 'b-');               # Vgrajena erf - modra polna

     # Legenda, oznake osi, naslov, mreža (enako kot prej)
#legend('Relativna razlika med potenčno vrsto in vgrajeno erf funkcijo(potenčna vrsta)');
xlabel('z');
ylabel('Logaritem relativne razilka');
title({'Logaritem relativne razlike med potenčno vrsto', 'in vgrajeno erf funkcijo'});
grid on;
