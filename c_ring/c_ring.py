!rm tallies.xml tallies.out materials.xml geometry.xml settings.xml

import openmc
import openmc.deplete
from pathlib import Path
import math
import openmc_depletion_plotter
import re
import dot_out_generator
from dot_out_generator import convert_time_units, get_short_lived_limits, get_long_lived_limit, create_full_run_tallies, get_single_depletion_tallies, get_waste_class, make_flux_file

v = openmc.Material()
v.add_element('C', 1, percent_type='ao')
v.set_density('g/cm3', 2.26)

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


timesteps_and_source_rates = [
    (365*24*60*60, 1e20)
]
for i in range(12):
    timesteps_and_source_rates.append((30*24*60*60, 0))
for i in range(9):
    timesteps_and_source_rates.append((365*24*60*60, 0))
for i in range(9):
    timesteps_and_source_rates.append((10*365*24*60*60, 0))
for i in range(9):
    timesteps_and_source_rates.append((100*365*24*60*60, 0))
    
    
timesteps = [item[0] for item in timesteps_and_source_rates]
source_rates = [item[1] for item in timesteps_and_source_rates]

# VIZ
color_assignment = {sphere_cell: 'blue', void_cell: 'red'}

plot = geometry.plot(basis='xz', color_by='cell', colors=color_assignment)
plot.figure.savefig('xz-cell.png')

plot = geometry.plot(basis='xy', color_by='cell',  colors=color_assignment)
plot.figure.savefig('xy-cell.png')

plot = geometry.plot(basis='yz', color_by='cell',  colors=color_assignment)
plot.figure.savefig('yz-cell.png')

#Tallies

chain = openmc.deplete.Chain.from_xml(openmc.config['chain_file'])
nuclides = [nuclide.name for nuclide in chain.nuclides]

tallies = create_full_run_tallies(sphere_cell)

#Deplete

model = openmc.model.Model(geometry, materials, settings, tallies)

operator = openmc.deplete.Operator(
    model,
    normalization_mode='source-rate',
    chain_file=openmc.config['chain_file'],
    reduce_chain_level=5,
    reduce_chain=True
)

integrator = openmc.deplete.PredictorIntegrator(
    operator,
    timesteps,
    source_rates=source_rates
)

integrator.integrate()

#model.deplete(
#    timesteps,
#    source_rates=source_rates,
#    method = "predictor",
#    operator_kwargs={
#        "normalization_mode": "source-rate",
#        "chain_file": openmc.config['chain_file'],
#        "reduce_chain_level": 5,
#        "reduce_chain": True
#    }
#)
    
        
        
results = openmc.deplete.Results("depletion_results.h5")
prev_time = 0
total_mass = v.get_mass()
with open('v.out', 'w') as file:
    for i, time in enumerate(results.get_times(time_units = "h")):
        
        # Tallies for each depletion
        statepoint_file = f"openmc_simulation_n{i}.h5"
        damage_energy, heat_energy, dose_rate = get_single_depletion_tallies(statepoint_file)
            
        # Composition for each depletion
        d_time = time - prev_time
        prev_time = time
        converted_d_time, converted_d_time_unit = convert_time_units(d_time)
        converted_time, converted_time_unit = convert_time_units(time)
        print(f"\nDepletion step {i}. Time interval: {converted_d_time} {converted_d_time_unit}. Total time elapsed: {converted_time} {converted_time_unit}.", file = file)
        
        print("Nuclide   Atoms         Mass (grams)     Activity (Bq)  g-Energy (eV)  half life", file = file)
        long_lived_concentration = 0;
        short_lived_concentration_A = 0;
        short_lived_concentration_B = 0;
        short_lived_concentration_C = 0;
        for nuc in nuclides:
            try:
                 # Get the number of atoms for the nuclide at this depletion step
                _, num_atoms = results.get_atoms(v, nuc)
                
                if num_atoms[i] > 0:
                    mass_grams = (openmc.data.atomic_mass(nuc) * num_atoms[i]) / 6.02214076e23
                    activity = openmc.data.decay_constant(nuc) * num_atoms[i]
                    if openmc.data.decay_photon_energy(nuc) is not None:
                        g_energy_eV = openmc.data.decay_photon_energy(nuc).integral() * num_atoms[i]
                    else:
                        g_energy_eV = 0
                    h_life = openmc.data.half_life(nuc)
                    h_life_text = "Stable"
                    if h_life != None:
                        h_life_num, h_life_unit = convert_time_units(h_life)
                        h_life_text = str(f"{h_life_num:3g}{h_life_unit}")
                    space1 = " " * (10-len(nuc))
                    space2 = " " * (14-len(str(f"{num_atoms[i]:5g}")))
                    space3 = " " * (17-len(str(f"{mass_grams:4g}")))
                    space4 = " " * (15-len(str(f"{activity:3g}")))
                    space5 = " " * (15-len(str(f"{g_energy_eV:3g}")))
                    print(f"{nuc}{space1}{num_atoms[i]:5g}{space2}{mass_grams:4g}{space3}{activity:3g}{space4}{g_energy_eV:3g}{space5}{h_life_text}", file = file)
                    activity_bq_per_cm_3 = activity / v.volume
                    activity_cu_per_m_3 = (activity_bq_per_cm_3 / 3.7e10) * 1e6
                    long_lived_limit = get_long_lived_limit(nuc)
                    sll_A, sll_B, sll_C = get_short_lived_limits(nuc)
                    if sll_A > 0:
                        contribution = activity_cu_per_m_3 / sll_A
                        if contribution > 0:
                            short_lived_concentration_A += contribution
                    if sll_B > 0:
                        contribution = activity_cu_per_m_3 / sll_B
                        if contribution > 0:
                            short_lived_concentration_B += contribution
                    if sll_C > 0:
                        contribution = activity_cu_per_m_3 / sll_C
                        if contribution > 0:
                            short_lived_concentration_C += contribution
                    if long_lived_limit > 0:
                        long_lived_contribution = activity_cu_per_m_3 / long_lived_limit
                        if long_lived_contribution > 0:
                            long_lived_concentration += long_lived_contribution
            except Exception as e:
                pass
        waste_class = get_waste_class(short_lived_concentration_A, short_lived_concentration_B, short_lived_concentration_C, long_lived_concentration)
        print(f"Nuclear waste classification: {waste_class}\n", file = file)

            
make_flux_file("statepoint.2.h5", "Graphite Sphere Shell Model")