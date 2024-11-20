import cftime
import os
import fabmos.transport.tmm
import fabmos

tm_config_dir = "."  # directory with a TM configuration from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
calendar = "360_day"  # any valid calendar recognized by cftime, see https://cfconventions.org/cf-conventions/cf-conventions.html#calendar

script_dir = os.path.dirname(__file__)

fabm_yaml = os.path.join(
    script_dir, "fabm.yaml"
)

domain = fabmos.transport.tmm.create_domain(tm_config_dir)

sim = fabmos.transport.tmm.Simulator(
    domain, calendar=calendar, fabm_config=fabm_yaml
)

sim.fabm.get_dependency("mole_fraction_of_carbon_dioxide_in_air").set(390.0)
sim.fabm.get_dependency("surface_air_pressure").set(101325.0)
#sim.fabm.get_dependency("absorption_of_silt").set(0.02)
#sim.fabm.get_dependency("bottom_stress").set(0.0)

out = sim.output_manager.add_netcdf_file(
    "output.nc", interval=1, interval_units=fabmos.TimeUnit.DAYS, save_initial=True
)
out.request(*sim.fabm.state_variables, time_average=True)

start = cftime.datetime(2000, 1, 1, calendar=calendar)
stop = cftime.datetime(2001, 1, 1, calendar=calendar)
sim.start(start, timestep=3600, transport_timestep=12 * 3600)
while sim.time < stop:
    sim.advance()
sim.finish()
