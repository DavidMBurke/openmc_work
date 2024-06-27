!rm tallies.xml tallies.out materials.xml geometry.xml settings.xml


import openmc
import openmc.deplete
from pathlib import Path
import math
import openmc_depletion_plotter
import re

#Vanadium test material
v = openmc.Material()
v.add_element('V', 1, percent_type='ao')
v.set_density('g/cm3', 6.1)

ring_thickness = 1 #cm
inner_rad = 114 #cm
outer_rad = inner_rad + ring_thickness
#Volume calc assuming ring as tall as thick
ring_vol = ring_thickness * math.pi * (outer_rad**2 - inner_rad ** 2)
v.volume = ring_vol
v.depletable = True
materials = openmc.Materials([v])
materials.export_to_xml()

# GEOMETRY
outer_sphere = openmc.Sphere(r=outer_rad, boundary_type='vacuum')
inner_cyl = openmc.ZCylinder(r=inner_rad)
outer_cyl = openmc.ZCylinder(r=outer_rad, boundary_type='reflective')
bottom_plane = openmc.ZPlane(z0=0)
top_plane = openmc.ZPlane(z0=1)

ring_region = +inner_cyl & -outer_cyl & +bottom_plane & -top_plane
ring_cell = openmc.Cell(region=ring_region)
ring_cell.fill = v
void_region = -outer_sphere & ~ring_region
void_cell = openmc.Cell(region=void_region)
void_cell.fill = None
geometry = openmc.Geometry([ring_cell, void_cell])


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
color_assignment = {ring_cell: 'blue', void_cell: 'red'}

plot = geometry.plot(basis='xz', color_by='cell', colors=color_assignment)
plot.figure.savefig('xz-cell.png')

plot = geometry.plot(basis='xy', color_by='cell',  colors=color_assignment)
plot.figure.savefig('xy-cell.png')

plot = geometry.plot(basis='yz', color_by='cell',  colors=color_assignment)
plot.figure.savefig('yz-cell.png')

#Tallies

chain = openmc.deplete.Chain.from_xml(openmc.config['chain_file'])
nuclides = [nuclide.name for nuclide in chain.nuclides]


energy_filter = openmc.EnergyFilter.from_group_structure('CCFE-709')
ccfe_tally = openmc.Tally(name="ccfe_tally")
ccfe_tally.filters.append(energy_filter)
ccfe_tally.scores.append('flux')
tallies = openmc.Tallies([ccfe_tally])

dose_tally = openmc.Tally(name="dose_rate")
dose_tally.filters = [openmc.CellFilter(ring_cell)]
dose_tally.scores = ['flux']
dose_tally.nuclides = ['total']
tallies.append(dose_tally)

heat_tally = openmc.Tally(name='heat')
heat_tally.filters = [openmc.CellFilter(ring_cell)]
heat_tally.scores = ['heating']
tallies.append(heat_tally)

damage_tally = openmc.Tally(name='damage_energy')
damage_tally.filters = [openmc.CellFilter(ring_cell)]
damage_tally.scores = ['damage-energy']
tallies.append(damage_tally)

tallies.export_to_xml()

#Deplete

model = openmc.model.Model(geometry, materials, settings, tallies)

# Deplete

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
    
def convert_time_units(time):
    if time > 24 * 365:
        return time / (24 * 365), 'y'
    elif time > 24:
        return time / 24, 'd'
    else:
        return time, 'h'

def get_short_lived_limits(nuc):
    isotope = nuc;
    half_life = openmc.data.half_life(nuc)
    if half_life is not None and half_life < (5 * 365 * 24 * 60 * 60):
        isotope = "sub_5_year" #nuclide w/ less than 5 year half life
    A, B, C = 0, 0, 0
    match isotope:
        case "H3":
            A = 40
        case "Co60":
            A = 700
        case "Ni63": #assuming in activated metal
            A = 35
        case "Sr90":
            A = .04
        case "Cs137":
            A = 1
        case "sub_5_year":
            A = 700
        case _:
            pass
    match isotope:
        case "Ni63": #assuming in activated metal
            B = 700
        case "Sr90":
            B = 150
        case "Cs137":
            B = 44
        case _:
            pass
    match isotope:
        case "Ni63": #assuming in activated metal
            C = 7000
        case "Sr90":
            C = 7000
        case "Cs137":
            C = 4600
        case _:
            pass
    return A, B, C
    
def get_long_lived_limit(nuc):
    match nuc:
        case 'C14':
            return 80 #assumes C14 in activated metal
        case 'Ni59':
            return 220
        case 'Nb94':
            return .2
        case 'Tc99':
            return 3
        case 'I129':
            return .08
        case 'Pu241':
            return 3500
        case 'Cm242':
            return 20000
        #alpha-emitting transuranic elements with half life > 5 year
        case 'Np237' | 'Pu238' | 'Pu239' | 'Pu240' | 'Pu242' | 'Pu244' | 'Am241' | 'Am243' | 'Cm243' | 'Cm244' | 'Cm245' | 'Cm246' | 'Cm247' | 'Cm248' | 'Cm250' | 'Bk247' | 'Cf249' | 'Cf250' | 'Cf251':
            return 100
        case _:
            return 0
        
results = openmc.deplete.Results("depletion_results.h5")
prev_time = 0
total_mass = v.get_mass()
with open('v.out', 'w') as file:
    for i, time in enumerate(results.get_times(time_units = "h")):
        
        # Tallies for each depletion
        statepoint_file = f"openmc_simulation_n{i}.h5"
        with openmc.StatePoint(statepoint_file) as sp:
            
            damage_energy_tally = sp.get_tally(name='damage_energy')
            damage_energy_value = 0
            if damage_energy_tally.num_realizations > 0:
                damage_energy_value = damage_energy_tally.get_values(scores=['damage-energy']).item()
            print("damage_energy: ", damage_energy_value)
            
            heat_energy_tally = sp.get_tally(name='heat')
            heat_energy = 0
            if heat_energy_tally.num_realizations > 0:
                heat_energy = heat_energy_tally.get_values(scores=['heating']).item()
            print("heat_energy: ", heat_energy)
            
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
                
                if nuc == "H3":
                    print("# atoms", num_atoms[i])
                if num_atoms[i] > 0:
                    mass_grams = v.get_mass(nuc)
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
                    space3 = " " * (17-len(str(f"{mass_grams:3g}")))
                    space4 = " " * (15-len(str(f"{activity:3g}")))
                    space5 = " " * (15-len(str(f"{g_energy_eV:3g}")))
                    print(f"{nuc}{space1}{num_atoms[i]:5g}{space2}{mass_grams:3g}{space3}{activity:3g}{space4}{g_energy_eV:3g}{space5}{h_life_text}", file = file)
                    activity_bq_per_cm_3 = activity / v.volume
                    activity_cu_per_m_3 = (activity_bq_per_cm_3 / 3.7e10) * 1e6
                    if nuc == "H3":
                        print(f"activity per m^3: {activity_cu_per_m_3}")
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
        waste_class = "A"
        print(f"Short lived concentration: {short_lived_concentration_A} by column 1")
        if short_lived_concentration_A > 1: 
            waste_class = "B"
            print(f"Short lived concentration: {short_lived_concentration_B} by column 2")
        if short_lived_concentration_B > 1:
            waste_class = "C"
            print(f"Short lived concentration: {short_lived_concentration_C} by column 3")
        if short_lived_concentration_C > 1:
            waste_class = "Not generally acceptable for near-surface disposal."
        if long_lived_concentration > .1 and (waste_class == "A" or waste_class == "B"):
            waste_class = "C"
        if long_lived_concentration > 1:
            waste_class = "Not generally acceptable for near-surface disposal."
        print(f"Long lived concentration: {long_lived_concentration}")
        print(f"Nuclear waste classification: {waste_class}\n")

            
with open("fluxes", 'w') as file:
    sp = openmc.StatePoint('statepoint.2.h5')
    flux_values = sp.get_tally(name='ccfe_tally').get_values(scores=['flux'])
    sp.close()
    for flux in flux_values:
        print(f"{flux[0].item():.18e}", file = file)