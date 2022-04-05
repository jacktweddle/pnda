import openmc
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np

#####################
##Material creation##
#####################

conc_mat = openmc.Material(name='concrete')
conc_mat.set_density('g/cm3', 2.3)
conc_mat.add_element('H',  0.168038, percent_type='ao')
conc_mat.add_element('O',  0.532000, percent_type='ao')
conc_mat.add_element('Na', 0.021365, percent_type='ao')
conc_mat.add_element('Al', 0.021343, percent_type='ao')
conc_mat.add_element('Si', 0.203231, percent_type='ao')
conc_mat.add_element('Ca', 0.018595, percent_type='ao')
conc_mat.add_element('Fe', 0.004246, percent_type='ao')
conc_mat.add_s_alpha_beta('c_H_in_H2O')

poly = openmc.Material(name='polyethylene')
poly.set_density('g/cm3', 0.96)
poly.add_element('C', 85.6, percent_type='wo')
poly.add_element('H', 14.4, percent_type='wo')
poly.add_s_alpha_beta('c_H_in_CH2')

hdpe = openmc.Material(name='HDPE')
hdpe.set_density('g/cm3', 1.07)
hdpe.add_element('C', 85.6, percent_type='wo')
hdpe.add_element('H', 14.4, percent_type='wo')
hdpe.add_s_alpha_beta('c_H_in_CH2')

air = openmc.Material(name='air')
air.set_density('g/cm3', 0.001225)
air.add_element('N', 0.784431)
air.add_element('O', 0.210748)
air.add_element('Ar', 0.0046)

steel = openmc.Material(name='steel')
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

mirrobor = openmc.Material(name='mirrobor')
mirrobor.set_density('g/cm3', 1.36)
mirrobor.add_nuclide('B10', 12.1, percent_type='wo')
mirrobor.add_nuclide('B11', 48.5, percent_type='wo')
mirrobor.add_nuclide('C12', 30.8, percent_type='wo')
mirrobor.add_nuclide('O16', 5, percent_type='wo')
mirrobor.add_nuclide('H1', 1, percent_type='wo')

helium3 = openmc.Material(name='He3')
helium3.set_density('g/cm3', 0.000165)
helium3.add_nuclide('He3', 100, percent_type='wo')

cadmium = openmc.Material(name='Cadmium')
cadmium.set_density('g/cm3', 4.047)
cadmium.add_element('Cd', 100, percent_type='wo')

mats = openmc.Materials([conc_mat, poly, hdpe, air, steel, mirrobor, helium3, cadmium])
mats.cross_sections = '/Users/jacktweddle/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents – Jack’s MacBook Pro/STFC/Birmingham Project/OpenMC/endfb80_hdf5.tmp/cross_sections.xml'
mats.export_to_xml()

############
##Geometry##
############

# Unshielded bunker, floor + outer air

floor_sur = openmc.ZPlane(0)
roof_sur = openmc.ZPlane(300)
roof_top_sur = openmc.ZPlane(400)
floor_bottom = openmc.ZPlane(-200, boundary_type='vacuum')
wall_left = openmc.XPlane(-180)
wall_right = openmc.XPlane(180)
wall_left_out = openmc.XPlane(-280)
wall_right_out = openmc.XPlane(280)
wall_down = openmc.YPlane(-455)
wall_up = openmc.YPlane(95)
wall_down_out = openmc.YPlane(-555)
wall_up_out = openmc.YPlane(555)
sur_door_cor_1 = openmc.YPlane(455)
sur_door_cor_2 = openmc.YPlane(325)
sur_cor_1 = openmc.XPlane(50)

outside_edge_left = openmc.XPlane(-1000, boundary_type='vacuum')
outside_edge_right = openmc.XPlane(1000, boundary_type='vacuum')
outside_edge_down = openmc.YPlane(-1000, boundary_type='vacuum')
outside_edge_up = openmc.YPlane(1000, boundary_type='vacuum')
outside_edge_top = openmc.XPlane(-1000, boundary_type='vacuum')

# Door

door_steel_right = openmc.XPlane(-273)
door_poly_left = openmc.XPlane(-278.5)
door_poly_right = openmc.XPlane(-274.5)

# Back shielding

back_shield_back = openmc.YPlane(-605)

# Roof shielding

roof_shield_back = openmc.YPlane(-355)
roof_shield_front = openmc.YPlane(145)
roof_shield_left = openmc.XPlane(-255)
roof_shield_right = openmc.XPlane(255)
roof_shield_top = openmc.ZPlane(450)

# Beamstop

beamstop_y_right = openmc.XPlane(-130)
beamstop_x_right = openmc.XPlane(50)
beamstop_x_back = openmc.YPlane(-305)
beamstop_front = openmc.YPlane(-275)
beamstop_top = openmc.ZPlane(180)

# Concrete shielding block

concrete_block_front = openmc.YPlane(-85)
concrete_block_right = openmc.XPlane(-10)
concrete_block_top = openmc.ZPlane(100)

# Poly blocks + thermal irradiation position

th_ir_left = openmc.XPlane(-160)
th_ir_right = openmc.XPlane(-135)
th_ir_front = openmc.YPlane(50)
th_ir_back = openmc.YPlane(75)

poly_l_right = openmc.XPlane(-102.5)
poly_l_top = openmc.ZPlane(180)
poly_r_left = openmc.XPlane(-59.5)
poly_r_front = openmc.YPlane(-30)

# Mirrobor

mb_l_right = openmc.XPlane(-102)
mb_r_left = openmc.XPlane(-60)
mb_back_front = openmc.YPlane(94.5)

# Trench

trench_front = openmc.YPlane(-185)
trench_bottom = openmc.ZPlane(-80)
trench_shield_back = openmc.YPlane(-100)
trench_shield_l_l_left = openmc.XPlane(-340)
trench_shield_l_l_right = openmc.XPlane(-310)
trench_shield_l_r_left = openmc.XPlane(-150)
trench_shield_l_r_right = openmc.XPlane(-120)
trench_shield_r_l_left = openmc.XPlane(120)
trench_shield_r_l_right = openmc.XPlane(150)
trench_shield_r_r_left = openmc.XPlane(310)
trench_shield_r_r_right = openmc.XPlane(340)

# Moderator cylinder

mod_cyl = openmc.ZCylinder(x0=-81, y0=-57.5, r=5)
mod_bottom = openmc.ZPlane(119)
mod_top = openmc.ZPlane(121)

# Shielding box - 0.254cm cadmium, 2.413cm poly

box_left_in = openmc.XPlane(-89)
box_right_in = openmc.XPlane(-47)
box_front_in = openmc.YPlane(-65.5)
box_back_in = openmc.YPlane(-49.5)
box_top_in = openmc.ZPlane(128)
box_bottom_in = openmc.ZPlane(112)

box_left_middle = openmc.XPlane(-89.254)
box_right_middle = openmc.XPlane(-46.746)
box_front_middle = openmc.YPlane(-65.754)
box_back_middle = openmc.YPlane(-49.246)
box_top_middle = openmc.ZPlane(128.254)
box_bottom_middle = openmc.ZPlane(111.746)

box_left_out = openmc.XPlane(-91.667)
box_right_out = openmc.XPlane(-44.333)
box_front_out = openmc.YPlane(-68.167)
box_back_out = openmc.YPlane(-46.833)
box_top_out = openmc.ZPlane(130.667)
box_bottom_out = openmc.ZPlane(109.333)

# Neutron detector

inner_casing = openmc.XCylinder(y0=-57.5, z0=120, r=1.549)
outer_casing = openmc.XCylinder(y0=-57.5, z0=120, r=1.60)
# Tube length = 15.39 cm, detector region length = 2.39 cm
tube_front = openmc.XPlane(-68.977)
tube_back = openmc.XPlane(-53.587)
det_front = openmc.XPlane(-67.707)
det_back = openmc.XPlane(-65.317)

# Regions

# Corridor

door_cor_r = -sur_door_cor_1 & +sur_door_cor_2 & +floor_sur & -roof_sur & -wall_right & +wall_left_out
corr_r = +wall_up & -sur_door_cor_2 & +floor_sur & -roof_sur & -wall_right & +sur_cor_1

# Door

door_poly_region = +door_poly_left & -door_poly_right & +sur_door_cor_2 & -sur_door_cor_1 & +floor_sur & -roof_sur
door_steel_region = +wall_left_out & -door_steel_right & +sur_door_cor_2 & -sur_door_cor_1 & +floor_sur & -roof_sur & ~door_poly_region

# Back shielding

back_shield_region = +back_shield_back & -wall_down_out & + wall_left & -wall_right & +floor_sur & -roof_sur

# Roof shielding

roof_shield_region = +roof_top_sur & -roof_shield_top & +roof_shield_left & -roof_shield_right & +roof_shield_back & -roof_shield_front

# Beamstop

beamstop_x = +beamstop_y_right & -beamstop_x_right & +beamstop_x_back & -beamstop_front & +floor_sur & -beamstop_top
beamstop_y = +wall_left & -beamstop_y_right & +wall_down & -beamstop_front & +floor_sur & -beamstop_top
beamstop_region = beamstop_x | beamstop_y

# Concrete base block

concrete_block_region = +concrete_block_front & -wall_up & +wall_left & -concrete_block_right & +floor_sur & -concrete_block_top

# Poly blocks + thermal irradiation position

th_ir_region = +th_ir_left & -th_ir_right & +th_ir_front & -th_ir_back & +concrete_block_top & -poly_l_top
poly_l_region = +wall_left & -poly_l_right & +concrete_block_front & -wall_up & +concrete_block_top & -poly_l_top & ~th_ir_region
poly_r_region = +poly_r_left & -concrete_block_right & +poly_r_front & -wall_up & +concrete_block_top & -poly_l_top

# Mirrobor

mb_l_region = +poly_l_right & -mb_l_right & +concrete_block_front & -wall_up & +concrete_block_top & -poly_l_top
mb_r_region = +mb_r_left & -poly_r_left & +poly_r_front & -wall_up & +concrete_block_top & -poly_l_top
mb_back_region = +mb_l_right & -mb_r_left & +mb_back_front & -wall_up & +concrete_block_top & -poly_l_top

mirrobor_region = mb_l_region | mb_r_region | mb_back_region

# Trench

trench_shield_l_l_region = +trench_shield_l_l_left & -trench_shield_l_l_right & +trench_bottom & -floor_sur & +trench_front & -trench_shield_back
trench_shield_l_r_region = +trench_shield_l_r_left & -trench_shield_l_r_right & +trench_bottom & -floor_sur & +trench_front & -trench_shield_back
trench_shield_r_l_region = +trench_shield_r_l_left & -trench_shield_r_l_right & +trench_bottom & -floor_sur & +trench_front & -trench_shield_back
trench_shield_r_r_region = +trench_shield_r_r_left & -trench_shield_r_r_right & +trench_bottom & -floor_sur & +trench_front & -trench_shield_back
trench_shield_left_wall_region = +wall_left_out & -wall_left & +trench_bottom & -floor_sur & +trench_front & -trench_shield_back
trench_shield_right_wall_region = +wall_right & -wall_right_out & +trench_bottom & -floor_sur & +trench_front & -trench_shield_back

trench_region = +trench_front & -concrete_block_front & +outside_edge_left & -outside_edge_right & +trench_bottom & -floor_sur & ~trench_shield_l_l_region & ~trench_shield_l_r_region & ~trench_shield_r_l_region & ~trench_shield_r_r_region

# Moderator cylinder

mod_region = -mod_cyl & +mod_bottom & -mod_top

# Neutron detector

steel_casing_region = +inner_casing & -outer_casing & +tube_front & -tube_back
det_region = -inner_casing & +det_front & -det_back
tube_region_1 = -inner_casing & +tube_front & -det_front
tube_region_2 = -inner_casing & +det_back & -tube_back

# Mod cyl box

mod_cyl_box_air_region = +box_left_in & -box_right_in & +box_front_in & -box_back_in & +box_bottom_in & -box_top_in & ~mod_region & ~steel_casing_region & ~det_region & ~tube_region_1 & ~tube_region_2
mod_cyl_box_cd_region = +box_left_middle & -box_right_middle & +box_front_middle & -box_back_middle & +box_bottom_middle & -box_top_middle & ~mod_cyl_box_air_region & ~mod_region & ~steel_casing_region & ~det_region & ~tube_region_1 & ~tube_region_2
mod_cyl_box_poly_region = +box_left_out & -box_right_out & +box_front_out & -box_back_out & +box_bottom_out & -box_top_out & ~mod_cyl_box_air_region & ~mod_cyl_box_cd_region & ~mod_region & ~steel_casing_region & ~det_region & ~tube_region_1 & ~tube_region_2

# Outside regions

outside_left_r = +outside_edge_left & -wall_left_out & +floor_sur & -roof_top_sur & -outside_edge_up & +outside_edge_down
outside_right_r = -outside_edge_right & +wall_right_out & +floor_sur & -roof_top_sur & -outside_edge_up & +outside_edge_down
above_r = +outside_edge_left & -outside_edge_right & +outside_edge_top & +roof_top_sur & -outside_edge_up & +outside_edge_down & ~roof_shield_region
outside_up_r = -outside_edge_right & +outside_edge_left & +floor_sur & -roof_top_sur & -outside_edge_up & +wall_up_out
outside_down_r = -outside_edge_right & +outside_edge_left & +floor_sur & -roof_top_sur & -wall_down_out & +outside_edge_down & ~back_shield_region
outside_r = outside_left_r | outside_right_r | above_r | outside_up_r | outside_down_r

# Walls + floor

bunker_inner_region = +floor_sur & -roof_sur & + wall_left & - wall_right & +wall_down & -wall_up & ~beamstop_region & ~concrete_block_region & ~th_ir_region & ~poly_l_region & ~poly_r_region & ~mirrobor_region & ~mod_region & ~steel_casing_region & ~det_region & ~tube_region_1 & ~tube_region_2 & ~mod_cyl_box_air_region & ~mod_cyl_box_cd_region & ~mod_cyl_box_poly_region
bunker_wall_region = +floor_sur & -roof_top_sur & + wall_left_out & - wall_right_out & +wall_down_out & -wall_up_out
bunker_wall_r = bunker_wall_region & ~bunker_inner_region & ~door_cor_r & ~corr_r & ~beamstop_region & ~concrete_block_region & ~th_ir_region & ~poly_l_region & ~poly_r_region & ~mirrobor_region & ~mod_region & ~steel_casing_region & ~det_region & ~tube_region_1 & ~tube_region_2 & ~mod_cyl_box_air_region & ~mod_cyl_box_cd_region & ~mod_cyl_box_poly_region
floor_r = +outside_edge_left & -outside_edge_right & -floor_sur & +floor_bottom & -outside_edge_up & +outside_edge_down & ~trench_region & ~trench_shield_l_l_region & ~trench_shield_l_r_region & ~trench_shield_r_l_region & ~trench_shield_r_r_region & ~trench_shield_left_wall_region & ~trench_shield_right_wall_region

# Cells
bunker_cell = openmc.Cell(region=bunker_inner_region)
bunker_wall = openmc.Cell(region=bunker_wall_r)
bunker_cell.fill = air
bunker_wall.fill = conc_mat

floor = openmc.Cell(region=floor_r)
floor.fill = conc_mat

corridor_1 = openmc.Cell(region=door_cor_r)
corridor_2 = openmc.Cell(region=corr_r)
corridor_1.fill = air
corridor_2.fill = air

outside = openmc.Cell(region=outside_r)
outside.fill = air

door_poly = openmc.Cell(region=door_poly_region)
door_poly.fill = poly

door_steel = openmc.Cell(region=door_steel_region)
door_steel.fill = steel

back_shield = openmc.Cell(region=back_shield_region)
back_shield.fill = conc_mat

roof_shield = openmc.Cell(region=roof_shield_region)
roof_shield.fill = conc_mat

beamstop = openmc.Cell(region=beamstop_region)
beamstop.fill = conc_mat

concrete_block = openmc.Cell(region=concrete_block_region)
concrete_block.fill = conc_mat

thermal_irradiation_position = openmc.Cell(region=th_ir_region)
thermal_irradiation_position.fill = air

poly_l = openmc.Cell(region=poly_l_region)
poly_l.fill = poly

poly_r = openmc.Cell(region=poly_r_region)
poly_r.fill = poly

mirrobor_cell = openmc.Cell(region=mirrobor_region)
mirrobor_cell.fill = mirrobor

trench_shield_l_l = openmc.Cell(region=trench_shield_l_l_region)
trench_shield_l_r = openmc.Cell(region=trench_shield_l_r_region)
trench_shield_r_l = openmc.Cell(region=trench_shield_r_l_region)
trench_shield_r_r = openmc.Cell(region=trench_shield_r_r_region)
trench_shield_left_wall = openmc.Cell(region=trench_shield_left_wall_region)
trench_shield_right_wall = openmc.Cell(region=trench_shield_right_wall_region)
trench = openmc.Cell(region=trench_region)

trench_shield_l_l.fill = conc_mat
trench_shield_l_r.fill = conc_mat
trench_shield_r_l.fill = conc_mat
trench_shield_r_r.fill = conc_mat
trench_shield_left_wall.fill = conc_mat
trench_shield_right_wall.fill = conc_mat
trench.fill = air

mod_cylinder = openmc.Cell(region=mod_region)
mod_cylinder.fill = hdpe

mod_cyl_box_air = openmc.Cell(region=mod_cyl_box_air_region)
mod_cyl_box_air.fill = air

mod_cyl_box_cd = openmc.Cell(region=mod_cyl_box_cd_region)
mod_cyl_box_cd.fill = cadmium

mod_cyl_box_poly = openmc.Cell(region=mod_cyl_box_poly_region)
mod_cyl_box_poly.fill = poly

steel_casing = openmc.Cell(region=steel_casing_region)
steel_casing.fill = steel
detector = openmc.Cell(region=det_region)
detector.fill = helium3
tube_1 = openmc.Cell(region=tube_region_1)
tube_1.fill = air
tube_2 = openmc.Cell(region=tube_region_2)
tube_2.fill = air

# makes a universe to contain all the cells
universe = openmc.Universe(cells=[bunker_cell, bunker_wall, floor, corridor_1, corridor_2, outside, door_poly, door_steel,
                                  back_shield, roof_shield, beamstop, concrete_block, thermal_irradiation_position,
                                  poly_l, poly_r, mirrobor_cell, trench_shield_l_l, trench_shield_l_r, trench_shield_r_l,
                                  trench_shield_r_r, trench_shield_left_wall, trench_shield_right_wall,
                                  trench, mod_cylinder, steel_casing, detector, tube_1, tube_2,
                                  mod_cyl_box_air, mod_cyl_box_cd, mod_cyl_box_poly])
geom = openmc.Geometry(universe)
geom.export_to_xml()

# shows the plots,
# note the new color scheme is based on materials not cells

'''
color_assignment = {bunker_wall: 'blue', bunker_cell: 'red', floor: 'grey', corridor_1: 'red',
                    corridor_2: 'red', outside: 'green', door_poly: 'black', door_steel: 'white', back_shield: 'purple', roof_shield: 'brown',
                    beamstop: 'pink', concrete_block: 'black', thermal_irradiation_position: 'black', poly_l: 'white', poly_r: 'white', mirrobor_cell: 'black',
                    trench_shield_l_l: 'black', trench_shield_l_r: 'black', trench_shield_r_l: 'black', trench_shield_r_r: 'black',
                    trench_shield_left_wall: 'black', trench_shield_left_wall: 'black', trench: 'white', mod_cylinder: 'purple',
                    steel_casing: 'grey', detector: 'green', tube_1: 'blue', tube_2: 'blue'}
# note the additional argument color_by, normally this defaults to 'cell'
plt.show(universe.plot(width=(2200, 2200), basis='xz', colors=color_assignment))
plt.show(universe.plot(width=(2200, 2200), basis='xy', colors=color_assignment))
plt.show(universe.plot(width=(2200, 2200), basis='yz', colors=color_assignment))



# Can also plot by colouring by material
color_assignment = {conc_mat: 'blue'}

plt.show(universe.plot(width=(1200, 1200), basis='xz', color_by='material', colors=color_assignment))
# plt.show(universe.plot(width=(1200, 1200), basis='xy', color_by='material', colors=color_assignment))
# plt.show(universe.plot(width=(1200, 1200), basis='yz', color_by='material', colors=color_assignment))
'''

#######################
##Simulation settings##
#######################

sett = openmc.Settings()
sett.batches = 10
sett.inactive = 0
sett.particles = 10000
sett.run_mode = 'fixed source'
sett.output = {'tallies': True}

source = openmc.Source()
source.space = openmc.stats.Point((-81, 2, 120))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete(14.1E6, 1)
source.time = openmc.stats.Uniform(0, 1E-3)  # 1 millisecond pulse
source.particle = 'neutron'
sett.source = source
sett.export_to_xml()

'''
# Create mesh which will be used for tally
mesh = openmc.RegularMesh()
mesh_height = 100   # number of cells in the X and Z dimensions
mesh_width = mesh_height
mesh_length = mesh_height
mesh.dimension = [mesh_width, mesh_length, 1]  # only 1 cell in the Z dimension
mesh.lower_left = [-500, -500, 110]   # physical limits (corners) of the mesh
mesh.upper_right = [1000, 1000, 130]
'''

###########
##Tallies##
###########

# Detection filters

time_bins = np.linspace(0, 2E-3, 10000)
energy_bins = np.linspace(0, 14.1E6, 1000)

particle_filter = openmc.ParticleFilter(['neutron'])
time_filter = openmc.TimeFilter(time_bins)
energy_filter = openmc.EnergyFilter(energy_bins)
surface_filter_1 = openmc.SurfaceFilter(det_front)

mesh = openmc.RegularMesh()
mesh.dimension = [250, 250, 1]
mesh.lower_left = [-180, -455, 120]
mesh.upper_right = [180, 455, 120]
mesh_filter = openmc.MeshFilter(mesh)

# Tally settings

tallies = openmc.Tallies()

tally1 = openmc.Tally(name='detector')
tally1.filters = [particle_filter, time_filter, surface_filter_1]
tally1.scores = ['current']
tallies.append(tally1)

tally2 = openmc.Tally(name='mesh_flux')
tally2.filters = [particle_filter, mesh_filter]
tally2.scores = ['flux']
tallies.append(tally2)

tallies.export_to_xml()

openmc.run()

'''
tallies = openmc.Tallies()
# Create mesh filter for tally
mesh_filter = openmc.MeshFilter(mesh)
mesh_tally = openmc.Tally(name='tallies_on_mesh')
mesh_tally.filters = [mesh_filter]
mesh_tally.scores = ['flux']
tallies.append(mesh_tally)
'''
