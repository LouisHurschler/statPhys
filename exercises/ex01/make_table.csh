##
#


set numatoms = 500
set area = 2500
set temperature = 300

set rgas = 8.314E-3

echo | awk '{ printf "%-5s %-8s %-8s -> %-12s : %-12s : %-8s\n", "N","A","T","P","R","PA/(NRT)"}'
echo | awk '{ printf "%-5s %-8s %-8s -> %-12s : %-12s : %-8s\n", "[-]","[nm^2]","[K]","[kJ/(K mol)]","[kJ/(K mol)]","[-]"}'

echo

foreach tmp_numatoms ( 300 500 700 )
  g++ -o gas -DNUMATOMS=$tmp_numatoms -DAREA=$area -DTEMPERATURE=$temperature gas_mac.cc
  set pressure = `gas | tail -1 | sed 's/^.*ave //g' | awk '{print $1}'`
  echo $tmp_numatoms $area $temperature $pressure $rgas |\
       awk '{ ratio = $2*$4/($1*$5*$3); printf "%-5d %-8.1f %-8.1f -> %-12.6f : %-12.6f : %-8.6f\n", $1,$2,$3,$4,$5,ratio}'
end

echo

foreach tmp_area ( 2000 2500 3000 )
  g++ -o gas -DNUMATOMS=$numatoms -DAREA=$tmp_area -DTEMPERATURE=$temperature gas_mac.cc
  set pressure = `gas | tail -1 | sed 's/^.*ave //g' | awk '{print $1}'`
  echo $numatoms $tmp_area $temperature $pressure $rgas |\
       awk '{ ratio = $2*$4/($1*$5*$3); printf "%-5d %-8.1f %-8.1f -> %-12.6f : %-12.6f : %-8.6f\n", $1,$2,$3,$4,$5,ratio}'
end

echo

foreach tmp_temperature ( 200 300 400 )
  g++ -o gas -DNUMATOMS=$numatoms -DAREA=$area -DTEMPERATURE=$tmp_temperature gas_mac.cc
  set pressure = `gas | tail -1 | sed 's/^.*ave //g' | awk '{print $1}'`
  echo $numatoms $area $tmp_temperature $pressure $rgas |\
       awk '{ ratio = $2*$4/($1*$5*$3); printf "%-5d %-8.1f %-8.1f -> %-12.6f : %-12.6f : %-8.6f\n", $1,$2,$3,$4,$5,ratio}'
end










