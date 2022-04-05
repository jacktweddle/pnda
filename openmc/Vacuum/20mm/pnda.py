import openmc
import numpy as np
import matplotlib.pyplot as plt

# Material creation

poly = openmc.Material(name='polyethylene')
poly.set_density('g/cm3', 1.07)
poly.add_element('C', 85.6, percent_type='wo')
poly.add_element('H', 14.4, percent_type='wo')
poly.add_s_alpha_beta('c_H_in_CH2')

air = openmc.Material(name='air')
air.set_density('g/cm3', 0.001225)
air.add_element('N', 0.784431)
air.add_element('O', 0.210748)
air.add_element('Ar', 0.0046)

steel = openmc.Material(name='Steel')
steel.set_density('g/cm3', 8.0)
steel.add_nuclide('C12', 0.02, percent_type='wo')
steel.add_nuclide('Si28', 0.2297, percent_type='wo')
steel.add_nuclide('Si29', 0.0121, percent_type='wo')
steel.add_nuclide('Si30', 0.0082, percent_type='wo')
steel.add_nuclide('P31', 0.0115, percent_type='wo')
steel.add_nuclide('S32', 0.0071, percent_type='wo')
steel.add_nuclide('S33', 0.0001, percent_type='wo')
steel.add_nuclide('S34', 0.0003, percent_type='wo')
steel.add_nuclide('Cr50', 0.3965, percent_type='wo')
steel.add_nuclide('Cr52', 57.9515, percent_type='wo')
steel.add_nuclide('Cr53', 0.919, percent_type='wo')
steel.add_nuclide('Cr54', 0.2331, percent_type='wo')
steel.add_nuclide('Mn55', 0.5, percent_type='wo')
steel.add_nuclide('Fe54', 1.9809, percent_type='wo')
steel.add_nuclide('Fe56', 32.245, percent_type='wo')
steel.add_nuclide('Fe57', 0.758, percent_type='wo')
steel.add_nuclide('Fe58', 0.1026, percent_type='wo')
steel.add_nuclide('Ni58', 3.1079, percent_type='wo')
steel.add_nuclide('Ni60', 1.2384, percent_type='wo')
steel.add_nuclide('Ni61', 0.0547, percent_type='wo')
steel.add_nuclide('Ni62', 0.1773, percent_type='wo')
steel.add_nuclide('Ni64', 0.0466, percent_type='wo')

helium3 = openmc.Material(name='Helium3')
helium3.set_density('g/cm3', 0.000165)
helium3.add_nuclide('He3', 100, percent_type='wo')

mats = openmc.Materials([poly, air, steel, helium3])
mats.cross_sections = '/Users/jacktweddle/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents – Jack’s MacBook Pro/STFC/Birmingham Project/OpenMC/endfb80_hdf5.tmp/cross_sections.xml'
mats.export_to_xml()

# Geometry definition

# Neutron detector

inner_casing_1 = openmc.XCylinder(y0=0, z0=0, r=1.549)
outer_casing_1 = openmc.XCylinder(y0=0, z0=0, r=1.60)
tube_front_1 = openmc.XPlane(15)
tube_back_1 = openmc.XPlane(30.39)
det_front_1 = openmc.XPlane(16.27)
det_back_1 = openmc.XPlane(18.66)

steel_casing_region_1 = +inner_casing_1 & -outer_casing_1 & +tube_front_1 & -tube_back_1
det_region_1 = -inner_casing_1 & +det_front_1 & -det_back_1
tube_region_1_1 = -inner_casing_1 & +tube_front_1 & -det_front_1
tube_region_2_1 = -inner_casing_1 & +det_back_1 & -tube_back_1

steel_casing_cell_1 = openmc.Cell(region=steel_casing_region_1)
steel_casing_cell_1.fill = steel
det_cell_1 = openmc.Cell(region=det_region_1)
det_cell_1.fill = helium3
tube_cell_1_1 = openmc.Cell(region=tube_region_1_1)
tube_cell_1_1.fill = air
tube_cell_2_1 = openmc.Cell(region=tube_region_2_1)
tube_cell_2_1.fill = air

# Moderator and bounding sphere

mod_cyl = openmc.ZCylinder(0, 0, r=5)
mod_bottom = openmc.ZPlane(-1)
mod_top = openmc.ZPlane(1)
mod_region = -mod_cyl & +mod_bottom & -mod_top
mod_cell = openmc.Cell(region=mod_region)
mod_cell.fill = poly

bounding_sphere = openmc.Sphere(r=100, boundary_type='vacuum')
air_region = -bounding_sphere & ~mod_region & ~det_region_1 & ~steel_casing_region_1 & ~tube_region_1_1 & ~tube_region_2_1
air_cell = openmc.Cell(region=air_region)
# air_cell.fill = air

universe = openmc.Universe(cells=[mod_cell, det_cell_1, steel_casing_cell_1, tube_cell_1_1, tube_cell_2_1, air_cell])
geom = openmc.Geometry(universe)
geom.export_to_xml()


colours = {mod_cell: 'red', steel_casing_cell_1: 'blue', tube_cell_1_1: 'brown', tube_cell_2_1: 'purple', det_cell_1: 'green', air_cell: 'grey'}
plt.show(universe.plot(width=(50, 50), basis='xy', colors=colours))
plt.show(universe.plot(width=(50, 50), basis='yz', colors=colours))
plt.show(universe.plot(width=(50, 50), basis='xz', colors=colours))


# Simulation settings

sett = openmc.Settings()
sett.batches = 10
sett.inactive = 0
sett.particles = 1000000
sett.run_mode = 'fixed source'
sett.output = {'tallies': True}

# Source settings

source = openmc.Source()
source.space = openmc.stats.Point((-20, 0, 0))
#source.angle = openmc.stats.Isotropic()
source.angle = openmc.stats.Monodirectional(reference_uvw=[1.0, 0.0, 0.0])
source.energy = openmc.stats.Discrete(14.1E6, 1)
source.time = openmc.stats.Uniform(0, 1E-3)  # 1 millisecond pulse
source.particle = 'neutron'
sett.source = source
sett.export_to_xml()

# Detection filters

time_bins = np.linspace(0, 2E-3, 10000)
energy_bins = np.linspace(0, 14.1E6, 1000)

particle_filter = openmc.ParticleFilter(['neutron'])
time_filter = openmc.TimeFilter(time_bins)
energy_filter = openmc.EnergyFilter(energy_bins)
surface_filter_1 = openmc.SurfaceFilter(det_front_1)

mesh = openmc.RegularMesh()
mesh.dimension = [100, 100]
mesh.lower_left = [-25, -25]
mesh.upper_right = [25, 25]
mesh_filter = openmc.MeshFilter(mesh)

# Tally settings

tallies = openmc.Tallies()

tally1 = openmc.Tally(name='detector1')
tally1.filters = [particle_filter, time_filter, surface_filter_1]
tally1.scores = ['current']
tallies.append(tally1)

tally2 = openmc.Tally(name='mesh_flux')
tally2.filters = [particle_filter, mesh_filter]
tally2.scores = ['flux']
tallies.append(tally2)

tallies.export_to_xml()

openmc.run()
