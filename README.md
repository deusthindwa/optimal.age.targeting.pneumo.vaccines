# Optimal age targeting for pneumococcal vaccination in older adults

Estimating the impact of introducing PCV13/15/20 or PPV23 in the over-55 populations of England, Brazil, Malawi and South Africa.

Code by Deus Thindwa and Sam Clifford under the supervision of Stefan Flasche and in collaboration with John Ojal and Peter McIntyre.

The data stored in `/data` are loaded and analysed by the code in `/script`, with resulting figures and tables created in `/output`.

Briefly, the code fits an exponential increase in incidence of _Streptococcus pneumoniae_ by age for each data set and simulates predicted incidences. These are then combined with age-fractionated population data and estimates of vaccine efficacy by time since vaccination (with optional age-dependence) to obtain estimates of the number of averted cases given the vaccine programme (both per vaccinee within the population).
