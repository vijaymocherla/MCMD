
bin='../bin/'

cd src/

gfortran -O3 init.f90 initialize.f90 io_module.f90 -o $bin'init.x'

gfortran -O3 mc_nvt_lj.f90 montecarlo.f90 io_module.f90 -o $bin'mc_nvt_lj.x'

gfortran -O3 md_nve_lj.f90 moldyn.f90 io_module.f90 -o $bin'md_nve_lj.x'

gfortran -O3 radialdist.f90 io_module.f90 -o $bin'radialdist.x'

gfortran -O3 diffusion.f90 io_module.f90 -o $bin'diff.x'

gfortran -O3 corr.f90 io_module.f90 -o $bin'corr.x'
