# this file documents name and runtime parameter changes from flash2 to flash3
#  for [/source/Simulation/SimulationMain/unitTest/Burn]

# use like this
# prompt> sed -f f2tof3.sed file_old.F90 > file_new.F90

s/t_min/tempMin/g
s/t_max/tempMax/g
s/rho_min/rhoMin/g
s/rho_max/rhoMax/g
s/compa/compA/g
s/compb/compB/g

s/shock_burning/useShockBurn/
s/ode_steper/odeStepper/
s/use_table/useBurnTable/

s/bn_enuc_factor/bn_enucDtFactor/
s/enuc_factor/enucDtFactor/

s/bn_conserved_var/bn_conservedVar/
s/conserved_var/conservedVar/
s/bn_tnucmin/bn_nuclearTempMin/
s/bn_tnucmax/bn_nuclearTempMax/
s/bn_dnucmin/bn_nuclearDensMin/
s/bn_dnucmax/bn_nuclearDensMax/
s/bn_ni56max/bn_nuclearNI56Max/
s/tnucmin/nuclearTempMin/
s/tnucmax/nuclearTempMax/
s/dnucmin/nuclearDensMin/
s/dnucmax/nuclearDensMax/
s/ni56max/nuclearNI56Max/

s/iburn/useBurn/

