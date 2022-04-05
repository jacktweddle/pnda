import openmc
import numpy as np
import matplotlib.pyplot as plt

"""
PNDA simulation replicating that produced by LLNL in MCNP as a means of benchmarking
OpenMC and assisting in the design and analysis of the PNDA experiment to be conducted in the
NILE bunker
"""

# Material creation

hdpe = openmc.Material(name='HDPE')
hdpe.set_density('g/cm3', 1.07)
hdpe.add_element('H', 11.7, percent_type='wo')
hdpe.add_element('O', 22.4, percent_type='wo')
hdpe.add_element('C', 61.73, percent_type='wo')
hdpe.add_nuclide('B10', 0.98, percent_type='wo')
hdpe.add_nuclide('B11', 4.02, percent_type='wo')
hdpe.add_s_alpha_beta('c_H_in_CH2')

aluminium = openmc.Material(name='Aluminium')
aluminium.set_density('g/cm3', 2.7)
aluminium.add_element('Al', 97.9, percent_type='wo')
aluminium.add_element('Si', 0.6, percent_type='wo')
aluminium.add_element('Mg', 1, percent_type='wo')
aluminium.add_element('Cr', 0.2, percent_type='wo')
aluminium.add_element('Cu', 0.28, percent_type='wo')

cadmium = openmc.Material(name='Cadmium')
cadmium.set_density('g/cm3', 4.047)
cadmium.add_element('Cd', 100, percent_type='wo')

helium3 = openmc.Material(name='Helium3')
helium3.set_density('g/cm3', 0.000165)
helium3.add_nuclide('He3', 100, percent_type='wo')

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

concrete = openmc.Material(name='Concrete')
concrete.set_density('g/cm3', 2.3)
concrete.add_nuclide('H1', 0.767, percent_type='wo')
concrete.add_nuclide('H2', 0.000115, percent_type='wo')
concrete.add_nuclide('Li6', 0.0000887, percent_type='wo')
concrete.add_nuclide('Li7', 0.00109, percent_type='wo')
concrete.add_nuclide('B10', 0.000339, percent_type='wo')
concrete.add_nuclide('B11', 0.000136, percent_type='wo')
concrete.add_nuclide('C12', 0.000697, percent_type='wo')
concrete.add_nuclide('O16', 46.6, percent_type='wo')
concrete.add_nuclide('Na23', 0.501, percent_type='wo')
concrete.add_element('Mg', 0.872, percent_type='wo')
concrete.add_nuclide('Al27', 4.18, percent_type='wo')
concrete.add_nuclide('Si28', 22.7, percent_type='wo')
concrete.add_nuclide('Si29', 1.15, percent_type='wo')
concrete.add_nuclide('Si30', 0.763, percent_type='wo')
concrete.add_nuclide('P31', 0.0513, percent_type='wo')
concrete.add_element('S', 0.152, percent_type='wo')
concrete.add_element('K', 1.6, percent_type='wo')
concrete.add_element('Ca', 14.2, percent_type='wo')
concrete.add_nuclide('Ti46', 0.0161, percent_type='wo')
concrete.add_nuclide('Ti47', 0.0147, percent_type='wo')
concrete.add_nuclide('Ti48', 0.148, percent_type='wo')
concrete.add_nuclide('Ti49', 0.011, percent_type='wo')
concrete.add_nuclide('Ti50', 0.0108, percent_type='wo')
concrete.add_element('V', 0.0108, percent_type='wo')
concrete.add_nuclide('Cr50', 0.00106, percent_type='wo')
concrete.add_nuclide('Cr52', 0.0205, percent_type='wo')
concrete.add_nuclide('Cr53', 0.00233, percent_type='wo')
concrete.add_nuclide('Cr54', 0.00058, percent_type='wo')
concrete.add_nuclide('Mn55', 0.0325, percent_type='wo')
concrete.add_nuclide('Fe54', 0.349, percent_type='wo')
concrete.add_nuclide('Fe56', 5.47, percent_type='wo')
concrete.add_nuclide('Fe57', 0.126, percent_type='wo')
concrete.add_nuclide('Fe58', 0.0167, percent_type='wo')
concrete.add_nuclide('Ni58', 0.0142, percent_type='wo')
concrete.add_nuclide('Ni60', 0.000543, percent_type='wo')
concrete.add_nuclide('Ni61', 0.000235, percent_type='wo')
concrete.add_nuclide('Ni62', 0.000747, percent_type='wo')
concrete.add_nuclide('Ni64', 0.000189, percent_type='wo')
concrete.add_nuclide('Cu63', 0.0115, percent_type='wo')
concrete.add_nuclide('Cu65', 0.00513, percent_type='wo')
concrete.add_nuclide('Zn65', 0.00221, percent_type='wo')
concrete.add_element('Zr', 0.00392, percent_type='wo')
concrete.add_nuclide('Mo92', 0.00168, percent_type='wo')
concrete.add_nuclide('Mo94', 0.00105, percent_type='wo')
concrete.add_nuclide('Mo95', 0.0018, percent_type='wo')
concrete.add_nuclide('Mo96', 0.00189, percent_type='wo')
concrete.add_nuclide('Mo97', 0.00108, percent_type='wo')
concrete.add_nuclide('Mo98', 0.00273, percent_type='wo')
concrete.add_nuclide('Mo100', 0.00109, percent_type='wo')
concrete.add_nuclide('Cd106', 0.0000162, percent_type='wo')
concrete.add_nuclide('Cd108', 0.0000115, percent_type='wo')
concrete.add_nuclide('Cd110', 0.000162, percent_type='wo')
concrete.add_nuclide('Cd111', 0.000166, percent_type='wo')
concrete.add_nuclide('Cd112', 0.000312, percent_type='wo')
concrete.add_nuclide('Cd113', 0.000158, percent_type='wo')
concrete.add_nuclide('Cd114', 0.000372, percent_type='wo')
concrete.add_nuclide('Cd116', 0.00000696, percent_type='wo')
concrete.add_nuclide('Ba138', 0.0665, percent_type='wo')

water = openmc.Material(name='Water')
water.set_density('g/cm3', 0.997)
water.add_element('H', 11.2098, percent_type='wo')
water.add_element('O', 88.7902, percent_type='wo')
water.add_s_alpha_beta("c_H_in_H2O")

air = openmc.Material(name='Air')
air.set_density('g/cm3', 0.001225)
air.add_element('N', 0.784431)
air.add_element('O', 0.210748)
air.add_element('Ar', 0.0046)

mats = openmc.Materials([hdpe, aluminium, cadmium, helium3, steel, water, air, concrete])
mats.cross_sections = '/Users/jacktweddle/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents – Jack’s MacBook Pro/STFC/Birmingham Project/OpenMC/endfb80_hdf5.tmp/cross_sections.xml'
mats.export_to_xml()

# Geometry definition

# Neutron detector

inner_casing = openmc.XCylinder(y0=-21.607, z0=0, r=1.549)
outer_casing = openmc.XCylinder(y0=-21.607, z0=0, r=1.60)
tube_front = openmc.XPlane(14)
tube_back = openmc.XPlane(29.39)
det_front = openmc.XPlane(15.27)
det_back = openmc.XPlane(17.66)

steel_casing_region = +inner_casing & -outer_casing & +tube_front & -tube_back
det_region = -inner_casing & +det_front & -det_back
tube_region_1 = -inner_casing & +tube_front & -det_front
tube_region_2 = -inner_casing & +det_back & -tube_back

steel_casing_cell = openmc.Cell(region=steel_casing_region)
steel_casing_cell.fill = steel
det_cell = openmc.Cell(region=det_region)
det_cell.fill = helium3
tube_cell_1 = openmc.Cell(region=tube_region_1)
tube_cell_1.fill = air
tube_cell_2 = openmc.Cell(region=tube_region_2)
tube_cell_2.fill = air

# Neutron source pipe

pipe = openmc.YCylinder(0, 0, r=5.08, boundary_type='vacuum')
pipe_bottom = openmc.YPlane(-40, boundary_type='vacuum')
pipe_top = openmc.YPlane(-30.353)
pipe_region = -pipe & -pipe_top & +pipe_bottom
pipe_cell = openmc.Cell(region=pipe_region)
pipe_cell.fill = air

# Moderator cylinder

R = 2.578 # Console 9

outer_mod_cyl = openmc.YCylinder(0, 0, r=R+0.1)
inner_mod_cyl = openmc.YCylinder(0, 0, r=R)
mod_cyl_outer_bottom = pipe_top
mod_cyl_inner_bottom = openmc.YPlane(-30.253)
mod_cyl_outer_top = openmc.YPlane(-12.353)
# mod_cyl_inner_top = openmc.YPlane(-12.453)

mod_region = -inner_mod_cyl & -mod_cyl_outer_top & +mod_cyl_inner_bottom
mod_cyl_region = -outer_mod_cyl & -mod_cyl_outer_top & +mod_cyl_outer_bottom & ~mod_region

mod_cell = openmc.Cell(region=mod_region)
mod_cell.fill = water
mod_cyl_cell = openmc.Cell(region=mod_cyl_region)
mod_cyl_cell.fill = aluminium

# Inner region (air)

air_left = openmc.XPlane(-30.353)
air_right = openmc.XPlane(30.353)
air_bottom = pipe_top
air_top = openmc.YPlane(30.353)
air_back = openmc.ZPlane(-30.353)
air_front = openmc.ZPlane(30.353)
air_region = +air_left & -air_right & +air_bottom & -air_top & +air_back & -air_front & ~mod_region & ~mod_cyl_region & ~steel_casing_region & ~det_region & ~tube_region_1 & ~tube_region_2
air_cell = openmc.Cell(region=air_region)
air_cell.fill = air

# Cadmium lining

cadmium_left = openmc.XPlane(-30.607)
cadmium_right = openmc.XPlane(30.607)
cadmium_bottom = openmc.YPlane(-30.607)
cadmium_top = openmc.YPlane(30.607)
cadmium_back = openmc.ZPlane(-30.607)
cadmium_front = openmc.ZPlane(30.607)
cadmium_region = +cadmium_left & -cadmium_right & +cadmium_bottom & -cadmium_top & +cadmium_back & -cadmium_front & ~pipe_region & ~mod_region & ~mod_cyl_region & ~air_region & ~steel_casing_region & ~det_region & ~tube_region_1 & ~tube_region_2
cadmium_cell = openmc.Cell(region=cadmium_region)
cadmium_cell.fill = cadmium

# Poly box

poly_left = openmc.XPlane(-33.02, boundary_type='vacuum')
poly_right = openmc.XPlane(33.02, boundary_type='vacuum')
poly_bottom = openmc.YPlane(-33.02, boundary_type='vacuum')
poly_top = openmc.YPlane(33.02, boundary_type='vacuum')
poly_back = openmc.ZPlane(-33.02, boundary_type='vacuum')
poly_front = openmc.ZPlane(33.02, boundary_type='vacuum')
poly_region = +poly_left & -poly_right & +poly_bottom & -poly_top & +poly_back & -poly_front & ~pipe_region & ~mod_cyl_region & ~mod_region & ~air_region & ~cadmium_region & ~steel_casing_region & ~det_region & ~tube_region_1 & ~tube_region_2
poly_cell = openmc.Cell(region=poly_region)
poly_cell.fill = hdpe

# Concrete bounding sphere

sphere_inner = openmc.Sphere(0, 0, 0, r=548.64)
sphere_outer = openmc.Sphere(0, 0, 0, r=609.6, boundary_type='vacuum')
concrete_region = +sphere_inner & -sphere_outer
concrete_cell = openmc.Cell(region=concrete_region)
concrete_cell.fill = concrete

universe = openmc.Universe(cells=[concrete_cell, pipe_cell, mod_cell, mod_cyl_cell, steel_casing_cell, det_cell, tube_cell_1, tube_cell_2, air_cell, cadmium_cell, poly_cell])
geom = openmc.Geometry(universe)
geom.export_to_xml()

colours = {concrete_cell: 'grey', pipe_cell: 'red', mod_cell: 'green', mod_cyl_cell: 'blue', steel_casing_cell: 'grey', det_cell: 'yellow', air_cell: 'white', cadmium_cell: 'purple', poly_cell: 'black'}
plt.show(universe.plot(width=(100, 100), basis='xy', colors=colours))
plt.show(universe.plot(width=(100, 100), basis='xz', colors=colours))

# Simulation settings

sett = openmc.Settings()
sett.batches = 20
sett.inactive = 0
sett.particles = 100000000
sett.run_mode = 'fixed source'
sett.output = {'tallies': True}

# Source settings

source = openmc.Source()
source.space = openmc.stats.Point((0, -35.894, 0))
'''
source.angle = openmc.stats.Isotropic()
#source.angle = openmc.stats.Monodirectional(reference_uvw=[0.0, 1.0, 0.0])
'''
aperture = 30.0
mu = openmc.stats.Uniform(np.cos(aperture/2), 1.0)
phi = openmc.stats.Uniform(0.0, 2*np.pi)
source.angle = openmc.stats.PolarAzimuthal(mu, phi, reference_uvw=(0., 1., 0.))
source.energy = openmc.stats.Discrete(2.5E6, 1)
source.time = openmc.stats.Uniform(0, 0.66E-3)  # 0.66 millisecond pulse
source.particle = 'neutron'
sett.source = source
sett.export_to_xml()

# Detection filters

time_bins = np.linspace(0, 3E-3, 10000)

particle_filter = openmc.ParticleFilter(['neutron'])
time_filter = openmc.TimeFilter(time_bins)
surface_filter = openmc.SurfaceFilter(det_front)

# Tally settings

tallies = openmc.Tallies()

tally = openmc.Tally(name='detector')
tally.filters = [particle_filter, time_filter, surface_filter]
tally.scores = ['current']
tallies.append(tally)
tallies.export_to_xml()

openmc.run()
