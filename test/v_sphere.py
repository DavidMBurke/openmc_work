
import openmc
import openmc.deplete
from pathlib import Path
import math
import openmc_depletion_plotter


openmc.config['chain_file'] = './chain_file.xml'

v = openmc.Material()
v.add_element('V', 1, percent_type='ao')
v.set_density('g/cm3', 6.1)

sphere_thickness = 1 #cm
inner_rad = 114 #cm
outer_rad = inner_rad + sphere_thickness
sphere_vol = (4/3) * math.pi * (outer_rad**3 - inner_rad ** 3)
v.volume = sphere_vol
v.depletable = True
materials = openmc.Materials([v])
materials.export_to_xml()


# GEOMETRY

inner_sphere = openmc.Sphere(r=inner_rad)
outer_sphere = openmc.Sphere(r=outer_rad, boundary_type='reflective')

sphere_region = +inner_sphere & -outer_sphere
sphere_cell = openmc.Cell(region=sphere_region)
sphere_cell.fill = v
void_region = -inner_sphere
void_cell = openmc.Cell(region=void_region)
void_cell.fill = None
geometry = openmc.Geometry([sphere_cell, void_cell])

# SOURCE
source = openmc.IndependentSource()
source.space = openmc.stats.Point((0,0,0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14.1e6], [1])
source.particles = 'neutron'

# SETTINGS

settings = openmc.Settings(
    batches = 2,
    inactive = 0,
    particles = 10000,
    source = source,
    run_mode = 'fixed source'
)

model = openmc.model.Model(geometry, materials, settings)

timesteps_and_source_rates = [
    (365*24*60*60, 1e20),
    (365*24*60*60, 1e20),
    (365*24*60*60, 1e20),
    (365*24*60*60, 0),
    (365*24*60*60, 0),
    (365*24*60*60, 0)
]

timesteps = [item[0] for item in timesteps_and_source_rates]
source_rates = [item[1] for item in timesteps_and_source_rates]


model.deplete(
    timesteps,
    source_rates=source_rates,
    method = "predictor",
    operator_kwargs={
        "normalization_mode": "source-rate",
        "chain_file": openmc.config['chain_file'],
        "reduce_chain_level": 5,
        "reduce_chain": True
    }
)
results = openmc.deplete.ResultsList.from_hdf5("depletion_results.h5")
results.plot_activity_vs_time()
