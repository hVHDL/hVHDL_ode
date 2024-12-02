# hVHDL_ode
numerical integrators and associated packages for dynamic simulation of control systems

Currently implemented are simplest basic integrators Runge-Kutta 1, 2 and 4th order and Adams-Moulton 2nd or 4th order

developed using open source NVC simulator

a oneliner to get the vhdl_ls.toml for syntax checking, running testbenches and plotting the resulting waveforms

> python vunit_run_ode.py -p 32 --export-json compiles.json ; python /from_vunit_export.py compiles.json ; python vunit_run_ode.py -p 32 ; python python/test_plot.py lcr_simulation_rk4_tb.dat

the "from_vunit_export.py" is taken from vhdl_ls repository example/

note, currently time is not properly set as part of the small steps of the rk integrators but it will be added
