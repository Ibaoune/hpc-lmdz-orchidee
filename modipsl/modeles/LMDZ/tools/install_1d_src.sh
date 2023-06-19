# This is a script in Bash.

# This script installs the source files of LMDZ 1D. Run this script
# from "libf".

set -xe

# cp --no-clobber phylmd/* phy1d:
cd phylmd
set +x
for my_file in *
  do
  if [[ ! -a ../phy1d/$my_file ]]
      then
#     cp $my_file ../phy1d
      cp -R $my_file ../phy1d 
  fi
done
set -x

cd ../phy1d

 ln -s ../dyn3d/mod_const_mpi.F90 .
 ln -sf ../dyn3d_common/infotrac.F90 .
 ln -s ../dyn3d/conf_gcm.F .
 ln -s ../dyn3d_common/control_mod.F90 .
 ln -sf ../dyn3d_common/q_sat.F .
 ln -s ../dyn3d_common/disvert.F90 .
