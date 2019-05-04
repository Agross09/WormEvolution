#!/usr/bin/env python3

# Copyright 2018 Alexis Pietak and Joel Grodstein
# See "LICENSE" for further details.

#    import pdb; pdb.set_trace()

import numpy as np
import sim as sim
import eplot as eplt
import edebug as edb
import re

#####################################################
# For an overview of the network data structures, see sim.py

# This is the function we're using for lab #1
def setup_lab1 (p):
    p.adaptive_timestep = False # Slow but accurate mode

    num_cells=4     # 4 cells

    # Gap junctions connect cells (we'll learn more soon).
    # No gap junctions means that the 4 cells are independent.
    n_GJs = 0

    # We've declared how many cells and GJs the simulation has. Create them
    # all, with default values.
    sim.init_big_arrays (num_cells, n_GJs, p)

    # To make it easier to set per-ion values below.
    Na = sim.ion_i['Na']; K = sim.ion_i['K']; Cl = sim.ion_i['Cl']

    # We want this lab to use a little different default values than the usual
    # ones. Specifically, a different diffusion constant for the K ion channels.
    # Note that this is the baseline value; simulation #3 will double it.
    sim.Dm_array[K,:] = 10.0e-18

    # Simulation #2 asks you to change the initial cell-internal concentrations.
    # Here's where you can do that.
    # The defaults (set in sim.init_big_arrays) are [Na+]=10, [K+]=125,
    # [Cl-]=55, [other-]=80 (units are all moles/m^3). Note that they add up to
    # charge neutral. Make sure to keep the cell interior pretty close to charge
    # neutral, or else the simulation will explode!
    # For example, if we want to double [Na+], it would go from 10 to 20. We
    # will also increase [Cl] by 10 moles/m^3 to preserve charge neutrality.
    # You should do similar things for the other cells, as per the instructions.
    #sim.cc_cells[Na,1] += 10   # Cell #1: [Na] is 20 rather than 10
    #sim.cc_cells[Cl,1] += 10   # [Cl] changes to keep charge neutrality.

    # Simulation #3 asks you to alter the diffusion constants D_Na, D_K and
    # D_Cl. By default, sim.init_big_arrays() sets all three to 1e-18, and then
    # we overrode D_K above to 10e-18. Here's where we override them for
    # simulation #3.
    # sim.Dm_array[Na,1] = 2.0e-18   # Cell #1 doubles D_Na

def setup_lab2 (p):
    p.adaptive_timestep = False # Slow but accurate mode
    p.sim_dump_interval = 1000  # Not much debug printout during the sim

    num_cells=4
    n_GJs = 0
    sim.init_big_arrays (num_cells, n_GJs, p)

    # By default, we have D_Na = D_K = 1e-18
    Na = sim.ion_i['Na']; K = sim.ion_i['K']
    sim.Dm_array[Na,0]= 20.0e-18    # +44 mV (reversal potential of Na+)
    sim.Dm_array[Na,1]=  5.0e-18    # +13mV
    sim.Dm_array[K,2] = 10.0e-18    # -57mV
    sim.Dm_array[K,3] = 160.0e-18   # -82mV (reversal potential of K+)

# This function shows that the QSS Vmem is insensitive to initial Vmem.
# 4 cells. Each one has Dm_Na set so that QSS Vmem will be -57mV.
# Tweak initial [Na] in cells 1-3 so that each starts out 5mV higher than
# the previous one. Nonetheless, they all converge to the same 30mV very
# quickly.
def setup_lab2b (p):
    p.adaptive_timestep = False # Slow but accurate mode
    # ... fill in the rest

# This is a *very* simple NN. There is only one layer, and no activation
# function. It merely validates that we can in fact compute a weighted sum.
# Cells 0, 1 and 2 are the input layer, each with a different Vmem.
# Cells 3, 4 and 5 are the output layer: each one weights a different
# combination of inputs. The weights are the GJ permeabilities; the inputs are
# the input-layer Vmem. It all works in QSS.
# def setup_lab_QSS_weighted_sum(p):
#     p.adaptive_timestep = False # Slow but accurate mode
#     num_cells=6
#     n_GJs = 7
#     sim.init_big_arrays (num_cells, n_GJs, p)
#     eplt.set_network_shape ([2, 3]) # For future pretty plotting.
#     GJ_scale = 1        # experiment with making all GJs bigger/smaller

#     # Cells 0, 1 and 2 are the input layer, each with a different Vmem.
#     Na = sim.ion_i['Na']; K=sim.ion_i['K']
#     sim.Dm_array[Na, 0]= 20e-18         # Vm =  44mV
#     sim.Dm_array[[Na,K], 1]= 20e-18         # Vm =  -1mV
#     sim.Dm_array[K,  2]= 20e-18         # Vm = -66mV

#     # Cell 3 = mix (0, 2) using GJ #0, #1.
#     sim.gj_connects[0] = (0, 3, 1*GJ_scale)
#     sim.gj_connects[1] = (2, 3, 1*GJ_scale)
#     sim.Dm_array[:,3]= 0    # No ion channels of its own, so it senses well

#     # Cell 4 = mix (0, 2) using GJ #2, #3.
#     sim.gj_connects[2] = (0, 4, ?*GJ_scale)
#     sim.gj_connects[3] = (2, 4, ?*GJ_scale)
#     sim.Dm_array[:,4]= 0    # No ion channels of its own, so it senses well

#     # Cell 5 = mix (0, 1, 2) using GJ #4, #5 and #6
#     sim.gj_connects[4] = (0, 5, ?*GJ_scale)
#     sim.gj_connects[5] = (1, 5, ?*GJ_scale)
#     sim.gj_connects[6] = (2, 5, ?*GJ_scale)
#     sim.Dm_array[:,5]= 0    # No ion channels of its own, so it senses well

# Finally, the worm. A few things to note:
# - It uses Alexis' parameters and chemicals.
# - We can switch it between having ion channels in the head/tail only, or in
#   the entire body.
def setup_lab_worm(params, worm):
    # Parameters:
    # N and kM are for the Hill model, where [M] controls the K channels.
    # 'scale' controls how strong the gap junctions are
    kM=worm.Km; N=worm.N; scale=worm.GJ_scale      # Lab sim #1
   # kM = 1; N = 2; scale =1  # Lab sim #2
    #kM=.2; N=5; scale=.1        # Lab sim #3

    eplt.set_network_shape([10, 1])

    num_cells= worm.num_cells
    params.sim_dump_interval = 500   # Not much debug printout during the sim

    n_GJs = num_cells-1

    sim.init_big_arrays (num_cells, n_GJs, params, ['M'])
    Na=sim.ion_i['Na']; K=sim.ion_i['K']; Cl=sim.ion_i['Cl']
    P=sim.ion_i['P']; M=sim.ion_i['M']

    # Alexis uses feedback on every cell; usually I only do it at the head/tail.
    fb_cells = [0,num_cells-1]

    # The ligand-gated K channel has diff = 1.7 e-17 m2/s when fully on.
    sim.Dm_array[K, fb_cells] = 1.7e-17
    sim.Dm_array[:,1:-1]= 0 # No ion channels in the body.

    # Ligand-gated K channels at the head and tail.
    for i in fb_cells:
        sim.ion_magic[K,i] = sim.magic_Hill_inv (M, N, kM, i)

    # Straight-line GJ connections along the worm. All are simple & ungated.
    sim.gj_connects['from']  = range(n_GJs)
    sim.gj_connects['to']    = sim.gj_connects['from'] + 1
    sim.gj_connects['scale'] = scale

    # Physical parameters of the 'mystery' ion.
    sim.z_array[M] = -1     # fixed ion valence
    sim.GJ_diffusion[M] = 1e-18 # diffusion rate through GJs

    # Initial concentrations: external {Na=145 mM, K=5mM, P=10mM}
    sim.cc_env[Na] = 145; sim.cc_env[K]=5; sim.cc_env[P]=10
    sim.cc_env[Cl] = 140    # To reach charge neutrality.

    # Initial concentrations: internal {Na=12, K=139, P=135}
    sim.cc_cells[Na] = 12; sim.cc_cells[K]=139; sim.cc_cells[P]=135
    sim.cc_cells[Cl] = 15   # To reach charge neutrality (off by 1 for M).

    # Ligand creation/decay/gating: creation at .1μM/s, and decaying
    # at (.1/s)*(its concentration)). This yields concentration=1uM.
    # However, we'll spread it from head to tail to seed the gradient.
    spread = 1
    spread = np.linspace (spread,-spread,num_cells)
    sim.cc_cells[M, :] = 1 + spread
    sim.cc_cells[Na,:] += spread

    # Ion-channel diffusion coef: Na=K=M = 1e-18 m2/s, and P=0. Already done.

    # Change the gap-junction length from 100nm to 15nm.
    params.gj_len = 15e-9

    # GJ scaled diff: Na=1.33e-17 m2/s, K=1.96e-17, M=1e-14, P=5e-17
    sim.GJ_diffusion[Na]=1.33e-17; sim.GJ_diffusion[K]=1.96e-17
    sim.GJ_diffusion[M] =1e-14;    sim.GJ_diffusion[P]=0

    #Na/K-ATPase pump with a maximum rate of 1.0x10-7 mol/(m^2 s)

    print ('Head D_K={}, tail D_K={}, scale={}'.format(sim.Dm_array[K,0],
                                                      sim.Dm_array[K,-1],scale))
    print ('Head magic={}'.format (sim.ion_magic[K,0]))
    print ('Tail magic={}'.format (sim.ion_magic[K,-1]))

def setup_and_sim(worm, time_to_run):
    import sys

    # If no command-line arguments are given, prompt for them.
    # if (len (sys.argv) <= 1):
    #     args = 'main.py ' + input('Arguments: ')
    #     sys.argv = re.split (' ', args)

    # regress_mode = (len(sys.argv) == 4) and (sys.argv[3]=="--regress")
    # if (regress_mode):
    #     sys.argv = sys.argv[0:3]

    # if (len(sys.argv) != 3):
    #    raise SyntaxError('Usage: python3 main.py test-name-to-run sim-end-time')

   # end_time = float (sys.argv[2])
    end_time = time_to_run

    # Run whatever test got chosen on the command line.
    GP = sim.Params()

    #I changed this to call setup_lab_worm every time and pass it a worm
    print(str(worm.Km))
    setup_lab_worm(GP, worm)


    # if (regress_mode):  # Works even if setup_...() overrules these params!
    #     GP.adaptive_timestep = False    # Force regression
    #     GP.sim_dump_interval = 1
    #     GP.sim_long_dump_interval = 10

    # Initialize Vmem -- typically 0, or close to that.
    Vm = sim.compute_Vm (sim.cc_cells, GP)
    assert (Vm<.5).all()    # I.e.,roughly charge neutral
    np.set_printoptions (formatter={'float': '{: 6.3g}'.format}, linewidth=90)
    print ('Initial Vm   ={}mV\n'.format(1000*Vm))

    print ('Starting main simulation loop')
    t_shots, cc_shots = sim.sim (end_time, GP)
    print ('Simulation is finished.')

    # Now, the simulation is over. Do any end-of-sim analysis.
    #edb.analyze_equiv_network (GP)

    # If your network contains cells that are interconnected by GJs, then
    # pretty_plot() can make a nice picture of the interconnections and each
    # cell's Vm.
    #Vmem post simulation
    Vm = sim.compute_Vm (sim.cc_cells, GP)
    eplt.pretty_plot (Vm*1e3)

    # We often want a printed dump of the final simulation results.
    np.set_printoptions (formatter={'float': '{:.6g}'.format}, linewidth=90)
    edb.dump (end_time, sim.cc_cells, edb.Units.mV_per_s, True)
    #edb.dump (end_time, sim.cc_cells, edb.Units.mol_per_m2s, True)

    # We often want graphs of various quantities vs. time.
    # If so, then comment out the quit() here, and then modify the code below
    # as desired.
    # quit()
    Na = sim.ion_i['Na']; K = sim.ion_i['K']; Cl=sim.ion_i['Cl']
    # P=sim.ion_i['P']; M=sim.ion_i['M']
    Vm_shots = [sim.compute_Vm (c,GP)*1000 for c in cc_shots]
    n_cells = sim.cc_cells.shape[1]
    eplt.plot_Vmem_graph(t_shots,Vm_shots, np.arange(n_cells),'Vmem(mV)')
    eplt.plot_Vmem_graph(t_shots,[s[Na] for s in cc_shots],np.arange(n_cells),'[Na] (mol/m3')
    eplt.plot_Vmem_graph(t_shots,[s[K]  for s in cc_shots],np.arange(n_cells),'[K] (mol/m3')
    eplt.plot_Vmem_graph(t_shots,[s[Cl]  for s in cc_shots],np.arange(n_cells),'[Cl] (mol/m3')
    #eplt.plot_Vmem_graph(t_shots,[s[M] for s in cc_shots],np.arange(n_cells),'[M] (mol/m3')