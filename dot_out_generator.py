import openmc

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

def create_full_run_tallies(cell):
    
    energy_filter = openmc.EnergyFilter.from_group_structure('CCFE-709')
    ccfe_tally = openmc.Tally(name="ccfe_tally")
    ccfe_tally.filters.append(energy_filter)
    ccfe_tally.scores.append('flux')
    tallies = openmc.Tallies([ccfe_tally])

    dose_tally = openmc.Tally(name="dose_rate")
    dose_tally.filters = [openmc.CellFilter(cell)]
    dose_tally.scores = ['flux']
    dose_tally.nuclides = ['total']
    tallies.append(dose_tally)

    heat_tally = openmc.Tally(name='heat')
    heat_tally.filters = [openmc.CellFilter(cell)]
    heat_tally.scores = ['heating']
    tallies.append(heat_tally)

    damage_tally = openmc.Tally(name='damage_energy')
    damage_tally.filters = [openmc.CellFilter(cell)]
    damage_tally.scores = ['damage-energy']
    tallies.append(damage_tally)

    tallies.export_to_xml()
    return tallies

def get_single_depletion_tallies(statepoint_file):
    damage_energy = 0
    heat_energy = 0
    dose_rate = 0
    with openmc.StatePoint(statepoint_file) as sp:

        damage_energy_tally = sp.get_tally(name='damage_energy')
        if damage_energy_tally.num_realizations > 0:
            damage_energy = damage_energy_tally.get_values(scores=['damage-energy']).item()
        print("damage energy: ", damage_energy)

        heat_energy_tally = sp.get_tally(name='heat')
        if heat_energy_tally.num_realizations > 0:
            heat_energy = heat_energy_tally.get_values(scores=['heating']).item()
        print("heat energy: ", heat_energy)

        dose_tally = sp.get_tally(name='dose_rate')
        if dose_tally.num_realizations > 0:
            dose_rate = dose_tally.get_values(scores=['flux']).item()
        print("dose rate: ", dose_rate)
        
    return damage_energy, heat_energy, dose_rate


def get_waste_class(slc_a, slc_b, slc_c, llc): #short and long lived concentrations
    waste_class = "A"
    print(f"Short lived concentration: {slc_a} by column 1")
    if slc_a > 1: 
        waste_class = "B"
        print(f"Short lived concentration: {slc_b} by column 2")
    if slc_b > 1:
        waste_class = "C"
        print(f"Short lived concentration: {slc_c} by column 3")
    if slc_c > 1:
        waste_class = "Not generally acceptable for near-surface disposal."
    if llc > .1 and (waste_class == "A" or waste_class == "B"):
        waste_class = "C"
    if llc > 1:
        waste_class = "Not generally acceptable for near-surface disposal."
    print(f"Long lived concentration: {llc}")
    
def make_flux_file(statepoint_file, name):
    with open("fluxes", 'w') as file:
        sp = openmc.StatePoint(statepoint_file)
        flux_values = sp.get_tally(name='ccfe_tally').get_values(scores=['flux'])
        sp.close()
        for flux in reversed(flux_values):
            print(f"{flux[0].item():.18e}", file = file)
        print(name, file = file)