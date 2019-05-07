# Copyright 2018 Alexis Pietak and Joel Grodstein
# See "LICENSE" for further details.

'''
Big picture: what does the network look like?
    The network is just a bunch of cells and GJs. However, there is no "class
    cell" or "class GJ".
    Instead, our data structure is most just a bunch of arrays.
    Some (like the membrane voltage, Vm) are 1D arrays[N_CELLS], with 1 element
    per cell. Others (like the ion concentration, cc_cells) are 2D
    arrays[N_IONS,N_CELLS], and can thus track (e.g.,) the concentration of
    every ion in every cell independently. By keeping the information in arrays,
    we get to use canned matrix-math routines in numpy for efficient simulation.

Which arrays track the cells?
    The first set of arrays are all [N_IONS,N_CELLS]; they have one row per ion,
    and one column per cell.
	cc_cells[i,c] = current concentration of ion #i in cell #c. It is given
	  in moles/m3, which is equivalent to mmoles/liter.
	Dm_array[i,c] = membrane ion-channel diffusion constant of cell #c for
	  ion #i (m2/s). One might think that Dm would just be the diffusion
	  constant of an ion channel itself, and would remain constant and be
	  mostly identical across all cells. Instead, Dm merges this with the
	  fraction of the cell's surface area that is covered by channels. So
	  if 1/3 of a cell membrane were covered by ion channels, then we would
	  cut that cell's Dm down by 3x. Also note that most ion channels use
          facilitated diffusion, and so these diffusion constants may be
          unrelated to an ion's "normal" diffusion constant in water.
    The next set are all arrays[N_CELLS]:
	z_array[i] = valence of ion #i. This array remains constant.
	cc_env [i] = concentration of ion #i in the ECF (moles/m3). In fact, it
	  remains constant also, since we assume the ECF is too big to be
	  affected by anything.
	Vm[i] = current Vmem of cell #i in volts. It is a derived array -- we
	  recompute it frequently from cc_cells and z_array.

    ion_i is a dict that maps an ion name into its row index in the 2D
	arrays above; e.g., ion_i['K'] => 1.

What about gap junctions (GJs)?
    To get a network more complex than just isolated individual cells, we need
    GJs. These, too, are just stored in a big array:
	gj_connects is a structured array [N_GJs]. Each element is one GJ.
          ['from']  and ['to'] specify the indices of its input and output cells
          (i.e., the column in the big arrays above).
    In addition to the connectivity of the GJ, we must also set how quickly
    the various ions travel through them. This happens in two parts.
    - First, GJ_diffusion[n_ions] gives the basic diffusion constant for each
      ion type through GJs. As with ion channels, the fraction of the cell
      membrane occupied by GJs is already accounted for by scaling down this
      number.
    - gj_connects['scale'] lets you, for any individual GJ, scale all of the
      ions' diffusion constants identically. By default, gj_connects['scale']
      is set to 1.

The 'magic' system for gating
    The 'magic' system is how we implement gated channels and gated GJs. We can
    gate based on either ligand concentration or Vmem. Why is this 'magic'?
    Because, for instructional purposes, we allow gatings that are completely
    non-physical. E.g., the concentration of a ligand in one cell can affect
    an ion-channel diffusion constant in a completely different cell!

    At the top level, we have two big arrays. Both get initialized to a default
    in init_big_arrays(), and then you can modify them if needed in setup_xxx().
	ion_magic[N_IONS, N_CELLS] describes the gating of all ion channels
	GJ_magic[N_GJS] describes the gating of all GJs.
    These arrays are both structured arrays of type magic_dtype, which has
    fields (in order):
      type: either 0 for ungated, 1 for Vmem gating, 2 for a Hill-model buffer,
	    or 3 for a Hill-model inverter.
      kM, N: the usual Hill model parameters
      cell, ion, cell2
   If 'type' is 0, the respective ion channel or GJ is ungated.
   If 'type'=1, we use voltage gating. First, get v1 and v2 from the Vmem of
     cells #cell and cell #cell2 respectively. We then compute
	scale = 1 / 1+exp(N*(v1-v2-kM))
     If either cell or cell2 are the special index -1, then that Vmem is assumed
     to be 0 (e.g., to just look at Vmem from one cell rather than a difference
     between two cells).
   If 'type' is 2 or 3, we use ligand-based gating: take ligand 'ion' from
     cell #'cell', and then use the usual Hill inverter model:
	scale = 1 / (1 + ([ion]/kM)^N).
     Types 2 and 3 are both a ligand-based Hill-inverter model; type 2 follows
     the Hill inverter model with "scale=1-scale" to make a buffer.
     'Cell2' is unused for both type=2 and type=3.

Generation and decay
    Bitsey implements a poor-man's GRN capability. For now, we merely implement
    generation and decay; between these and the magic system, we can implement
    Hill buffers and inverters. Furthermore, wire-OR then gives us multi-input
    gates.

    We implement this with several arrays:
    gen_cells[N_IONS,N_CELLS] (units of moles/(m3*s)) and gen_magic[same size]
	give generation parameters. For any ion i in cell c, the generation rate
	(moles/(m3*s)) is d[ion]/dt=gen_cells[i,c] * eval_magic(gen_magic[i,c])
    decay_rates[N_IONS]	gives the decay rate for each ion in units of 1/s.
	In any cell, we then have d[ion]/dt = [ion] * decay_rate, with resulting
	units of moles/(m3*s)
'''
import numpy as np
import scipy
import math	# To get pi.
import sim_toolbox as stb
import operator
import edebug as edb

#    import pdb; pdb.set_trace()

# Stuff all of the fundamental constants and commonly used parameters into a
# class. Any instance of this class will thus have all the constants.
class Params(object):
    def __init__(self):
        self.F = 96485  # Faraday constant [C/mol]
        self.R = 8.314  # Gas constant [J/K*mol]
        self.eo = 8.854e-12  # permittivity of free space [F/m]
        self.kb = 1.3806e-23  # Boltzmann constant [m2 kg/ s2 K1]
        self.q = 1.602e-19  # electron charge [C]
        self.tm = 7.5e-9  # thickness of cell membrane [nm]
        self.cm = 0.05  # patch capacitance of membrane [F/m2]
        self.T = 310  # temperature in Kelvin
        self.deltaGATP = -37000 # free energy released in ATP hydrolysis [J/mol]
        self.cATP = 1.5  # ATP concentration (mol/m3)
        self.cADP = 0.15  # ADP concentration (mol/m3)
        self.cPi = 0.15  # Pi concentration (mol/m3)
        self.alpha_NaK =1.0e-7 # max rate constant Na-K ATPase/unit surface area
        self.KmNK_Na = 12.0  # NaKATPase enzyme ext Na half-max sat value
        self.KmNK_K = 0.2  # NaKATPase enzyme ext K half-max sat value
        self.KmNK_ATP = 0.5  # NaKATPase enzyme ATP half-max sat value
        self.cell_r = 5.0e-6  # radius of single cell
        self.gj_len = 100e-9  # distance between two GJ connected cells [m]
        self.cell_sa = (4 * math.pi * self.cell_r ** 2)  # cell surface area
        self.cell_vol = ((4 / 3) * math.pi * self.cell_r ** 3)  # cell volume

        # Simulation control.
        self.sim_dump_interval=10
        self.sim_long_dump_interval=100
        self.no_dumps = False

        # Numerical-integration parameters. These place a limit on how much
        # any cell's Vmem, or any ion concentration, can change in one timestep.
        self.sim_integ_max_delt_Vm = .0001	# Volts/step
        # .001 means that no [ion] can change by more than .1% in one timestep.
        self.sim_integ_max_delt_cc = .001
        self.adaptive_timestep = True	# So that the params above get used.
        self.use_implicit = False	# Use the implicit integrator

        ###JJJ My new neutral ion "A"
        self.concA_fixed_amount = 100	# moles/m3

def init_big_arrays (n_cells, n_GJs, p, extra_ions=[]):
    global cc_cells, cc_env, Dm_array, z_array, ion_i,Vm, \
           GJ_diffusion, ion_magic, GJ_magic, gj_connects, GP, \
           gen_cells, gen_magic, decay_rates
    GP = p

    # ion properties (Name, base membrane diffusion [m2/s], valence
    #	initial concentration inside cell [mol/m3],
    #	fixed concentration outside cell [mol/m3], 
    # These are temporary structures. We use them to provide initial values for
    # the big arrays we are about to build, and to specify the order of which
    # row represents which ion in those arrays.
    Na={'Name':'Na', 'D_mem':1e-18, 'D_GJ':1e-18, 'z':1, 'c_in':10, 'c_out':145}
    K ={'Name':'K',  'D_mem':1e-18, 'D_GJ':1e-18, 'z':1, 'c_in':125,'c_out':5}
    Cl={'Name':'Cl', 'D_mem':1e-18, 'D_GJ':1e-18, 'z':-1,'c_in':55, 'c_out':140}
    P= {'Name':'P',  'D_mem':0,     'D_GJ':1e-18, 'z':-1,'c_in':80, 'c_out':10}

    # stack the above individual dictionaries into a list to make it easier to
    # process them in the loop below.
    ions_vect = [Na, K, Cl, P]

    # Any particular sim may want to declare extra ions.
    for ion in extra_ions:
        ions_vect.append ({'Name':ion, 'D_mem':0.0, 'D_GJ':1e-18,
                           'z':0, 'c_in':0,  'c_out':0})
    n_ions = len(ions_vect)

    cc_cells = np.empty ((n_ions, n_cells))
    Dm_array = np.empty ((n_ions, n_cells))
    z_array  = np.empty ((n_ions))
    cc_env   = np.empty ((n_ions))
    GJ_diffusion = np.empty ((n_ions))

    ion_i    = {}

    # Push the parameters of the above ions into the various arrays.
    for row, ion_obj in enumerate(ions_vect):
        cc_cells[row,:] = ion_obj['c_in']	# initial cell conc
        cc_env  [row]   = ion_obj['c_out']	# fixed environmental conc
        Dm_array [row] = ion_obj['D_mem']	# initial membrane diff coeff
        z_array[row] = ion_obj['z']		# fixed ion valence
        GJ_diffusion[row] = ion_obj['D_GJ']	# diffusion rate through GJs
        ion_i[ion_obj['Name']] = row		# map ion name -> its row

    # Initialize the magic arrays to their default no-magic state.
    magic_dtype=np.dtype ([('type','i4'),('kM','f4'), ('N','f4'), \
                           ('cell','i4'), ('ion','i4'), ('cell2','i4')])
    ion_magic = np.zeros ((n_ions, n_cells), dtype=magic_dtype)
    GJ_magic  = np.zeros ((n_GJs), dtype=magic_dtype)
    gen_magic = np.zeros ((n_ions, n_cells), dtype=magic_dtype)

    # Create default arrays for GJs, and for generation, decay rates.
    gj_connects=np.zeros((n_GJs), dtype=[('from','i4'),('to','i4'),('scale','f4')])
    gen_cells   = np.zeros ((n_ions, n_cells))
    decay_rates = np.zeros ((n_ions))

# Register a function to be called after every simulation step.
post_hook_func=None
def register_post_hook (func):
    global post_hook_func
    post_hook_func = func

# The main "do-it" simulation function.
# Takes the current cc_cells[n_ions,n_cells], does all of the physics work, and
# returns an array of concentration slew rates [n_ions,n_cells]; i.e.,
# moles/m3 per second.
# In normal operation, sim_slopes() is only called from sim.sim() and thus is
# not needed outside of sim.py. However, the solve() functions in main.py call
# sim_slopes() directly -- sim_slopes() == vector of zeroes indicates steady
# state.
def sim_slopes (t, Cc_cells):
    global cc_env, Dm_array, z_array, ion_i, Vm, gj_connects, GP, \
           gen_cells, gen_magic, decay_rates, cc_cells
    cc_cells = Cc_cells
    Vm = compute_Vm (cc_cells, GP)

    # General note: our units of flux are moles/(m2*s). The question: m2 of
    # what area? You might think that for, e.g., ion channels, it should be per
    # m2 of ion-channel area -- but it's not. All fluxes are per m2 of cell-
    # membrane area. Thus, the (e.g.,) diffusion rate through ion channels must
    # be scaled down by the fraction of membrane area occupied by channels.
    # The same goes for ion pumps and GJs.
    slew_cc = np.zeros (cc_cells.shape)	# Per-ion cell fluxes

    # Run the Na/K-ATPase ion pump in each cell.
    # Returns two 1D arrays[N_CELLS] of fluxes; units are moles/(m2*s)
    pump_Na,pump_K,_ = stb.pumpNaKATP(cc_cells[ion_i['Na']],cc_env[ion_i['Na']],
                                      cc_cells[ion_i['K']], cc_env[ion_i['K']],
                                      Vm, GP.T, GP, 1.0)

    # Kill the pumps on worm-interior cells (based on Dm=0 for all ions)
    keep_pumps = np.any (Dm_array>0, 0)	# array[n_cells]
    pump_Na *= keep_pumps
    pump_K  *= keep_pumps

    # Update the cell-interior [Na] and [K] after pumping (assume env is too big
    # to change its concentration).
    slew_cc[ion_i['Na']] = pump_Na
    slew_cc[ion_i['K']]  = pump_K

    # for each ion: (sorted to be in order 0,1,2,... rather than random)
    for ion_name,ion_index in sorted (ion_i.items(),key=operator.itemgetter(1)):
        # GHK flux across membranes into the cell
        # It returns array[N_CELLS] of moles/(m2*s)
        f_ED = GHK (ion_index, Vm)
        f_ED *= eval_magic (ion_magic[ion_index,:])
        slew_cc[ion_index] += f_ED

    # Gap-junction computations.
    deltaV_GJ = (Vm[gj_connects['to']] - Vm[gj_connects['from']]) # [n_GJs]

    # Get the gap-junction Norton-equivalent circuits for all ions at once.
    # Units of Ith are mol/(m2*s); units of Gth are mol/(m2*s) per Volt.
    (GJ_Ith, GJ_Gth) = GJ_norton(GP)	# Both are [n_ions,n_GJs].
    f_gj = GJ_Ith + deltaV_GJ*GJ_Gth	# [n_ions,n_GJs]

    magic = eval_magic (GJ_magic)		# [n_GJs]
    f_gj *= magic				# [n_ions,n_GJs]

    # Update cells with GJ flux:
    # Note that the simple slew_cc[ion_index, gj_connects['to']] += f_gj
    # doesn't actually work in the case of two GJs driving the same 'to'
    for ion_name,ion_index in sorted (ion_i.items(),key=operator.itemgetter(1)):
        # cell. Instead, we use np.add.at().
        np.add.at (slew_cc[ion_index,:], gj_connects['from'], -f_gj[ion_index])
        np.add.at (slew_cc[ion_index,:], gj_connects['to'],    f_gj[ion_index])

    # The current slew_cc units are moles/(m2*s), where the m2 is m2 of
    # cell-membrane area. To convert to moles/s entering the cell, we multiply
    # by the cell's surface area. Then, to convert to moles/m3 per s entering
    # the cell, we divide by the cell volume.
    slew_cc *= (GP.cell_sa / GP.cell_vol)

    # Next, do generation and decay.
    for ion_name,ion_index in sorted (ion_i.items(),key=operator.itemgetter(1)):
        gen = gen_cells[ion_index,:] * eval_magic(gen_magic[ion_index,:])
        decay = cc_cells[ion_index,:] * decay_rates[ion_index]
        slew_cc[ion_index] += gen - decay

    global post_hook_func
    if (post_hook_func != None):
        post_hook_func(t, GP, cc_cells, slew_cc)

    return (slew_cc)	# Moles/m3 per second.

# Given: per-cell, per-ion charges in moles/m3.
# First: sum them per-cell, scaled by valence to get "signed-moles/m3"
# Next: multiply by F to convert moles->coulombs. Multiply by cell volume/
# surface area to get coulombs/m2, and finally divide by Farads/m2.
# The final scaling factor is F * p.cell_vol / (p.cell_sa*p.cm),
# or about 3200 mV per (mol/m3)
def compute_Vm (Cc_cells, p):
    global cc_cells
    cc_cells = Cc_cells

    # Calculate Vmem from scratch via the charge in the cells.
    rho_cells = (cc_cells * z_array[:,np.newaxis]).sum(axis=0) * p.F
    return (rho_cells * p.cell_vol / (p.cell_sa*p.cm))

def GHK (ion_index, Vm):
    global cc_env, Dm_array, z_array, GP, cc_cells
    num_cells = cc_cells.shape[1]
    f_ED = stb.electroflux(cc_env[ion_index] * np.ones(num_cells),
                           cc_cells[ion_index],
                           Dm_array[ion_index],
                           GP.tm * np.ones(num_cells),
                           z_array[ion_index] * np.ones(num_cells),
                           Vm,
                           GP.T,
                           GP,
                           rho=np.ones(num_cells)
                           )
    return (f_ED)

def sim (end_time, p):
    if (p.use_implicit):
        return (sim_implicit (end_time, p))
    global cc_cells, Vm
    # Save snapshots of core variables for plotting.
    t_shots=[]; cc_shots=[]; last_shot=-100;

    # run the simulation loop:
    i=0; t=0
    time_step = .005
    while (t < end_time):
        slew_cc = sim_slopes(t,cc_cells)

        # Compute Vmem slew (in Volts/s). Essentially, it's just slew_Q/C.
        # Slew_cc is slew-flux in moles/m3 per second. We first convert to
        # moles/(m2 of cell-membrane cross-sec area).
        # Then sum (slew-moles * valence) to get a "slew signed moles."
        # Finally, multiply by F to get slew-Coulombs/(m2*s), and divide by
        # cap/m2 to get slew-Vmem/s.
        mult = (p.cell_vol / p.cell_sa) * (p.F/ p.cm)
        slew_Vm = (slew_cc * z_array[:,np.newaxis]).sum(axis=0) * mult

        # Timestep control.
        # max_volts / (volts/sec) => max_time
        max_t_Vm = p.sim_integ_max_delt_Vm / (np.absolute (slew_Vm).max())
        # (moles/m3*sec) / (moles/m3) => fractional_change / sec
        if (p.adaptive_timestep):
            frac_cc = np.absolute(slew_cc)/(cc_cells+.00001)
            max_t_cc = p.sim_integ_max_delt_cc / (frac_cc.max())
            n_steps = max (1, int (min (max_t_Vm, max_t_cc) / time_step))
            #print ('At t={}: max_t_Vm={}, max_t_cc={} => {} steps'.format(t, max_t_Vm, max_t_cc, n_steps))
            #print ('steps_Vm=', (.001/(time_step*np.absolute (slew_Vm))).astype(int))
        else:
            n_steps = 1

        cc_cells +=  slew_cc * n_steps * time_step

        # Calculate Vmem from scratch via the charge in the cells.
        Vm = compute_Vm (cc_cells,p)

        # Dump out status occasionally during the simulation.
        # Note that this may be irregular; numerical integration could, e.g.,
        # repeatedly do i += 7; so if sim_dump_interval=10 we would rarely dump!
        if ((i % p.sim_dump_interval == 0) and not p.no_dumps):
            long = (i % p.sim_long_dump_interval == 0)
            edb.dump (t, cc_cells, edb.Units.mV_per_s, long) # mol_per_m2s
            #edb.analyze_equiv_network (p)
            #edb.dump_magic ()

        i += n_steps
        t = i*time_step

        if (t>9000000):		# A hook to stop & debug during a sim.
            edb.debug_print_GJ (p, cc_cells, 1)
            import pdb; pdb.set_trace()
            print (sim_slopes (2000, cc_cells))

        # Save information for plotting at sample points. Early on (when things
        # are changing quickly) save lots of info. Afterwards, save seldom so
        # as to save memory (say 100 points before & 200 after)
        boundary=min (50,end_time);
        before=boundary/100; after=(end_time-boundary)/200
        interval = (before if t<boundary else after)
        if (t > last_shot+interval):
            t_shots.append(t)
            cc_shots.append(cc_cells.copy())
            last_shot = t

    return (t_shots, cc_shots)

# Replacement for sim(); it uses scipy.integrate.solve_ivp()
def sim_implicit (end_time, p):
    import scipy
    global cc_cells, Vm
    num_ions, num_cells = cc_cells.shape

    def wrap (t, y):
        global cc_cells

        print ('----------------\nt={:.9g}'.format(t))
        slew_cc = sim_slopes (t, y.reshape(num_ions,num_cells)) # moles/(m3*s)
        slew_cc = slew_cc.reshape (num_ions*num_cells)
        np.set_printoptions (formatter={'float':'{:6.2f}'.format},linewidth=120)
        print ('y={}'.format(y))
        np.set_printoptions (formatter={'float':'{:7.2g}'.format},linewidth=120)
        print ('slews={}'.format(slew_cc))
        return (slew_cc)

    # Save information for plotting at sample points. Early on (when things
    # are changing quickly) save lots of info. Afterwards, save seldom so
    # as to save memory. So, 100 points in t=[0,50], then 200 in [50, end_time].
    boundary=min (50,end_time)
    t_eval = np.linspace (0,boundary,50,endpoint=False)
    if (end_time>50):
        t_eval = np.append (t_eval, np.linspace (boundary, end_time, 200))

    # run the simulation loop:
    y0 = cc_cells.reshape (num_ions*num_cells)
    bunch = scipy.integrate.solve_ivp (wrap, (0,end_time), y0, method='BDF', \
                                       t_eval=t_eval)

    print ('{} func evals, status={} ({}), success={}'.format \
             (bunch.nfev, bunch.status, bunch.message, bunch.success))
    t_shots = t_eval.tolist()
    # bunch.y is [n_ions*n_cells, n_timepoints]
    cc_shots = [y.reshape((num_ions,num_cells)) for y in bunch.y.T]
    cc_cells = cc_shots[-1]
    return (t_shots, cc_shots)

# Builds and returns a Norton equivalent model for all GJs.
# Specifically, two arrays GJ_Ith and GJ_Gth of [n_ions,n_GJ].
# Ith[i,g] is the diffusive flux of ion #i in the direction of GJ[g].from->to,
# and has units (mol/m2*s)
# Gth*(Vto-Vfrom) is the drift flux of particles in the from->to direction;
# Gth has units (mol/m2*s) per Volt.
def GJ_norton (p):
    global cc_cells
    n_GJ = gj_connects.size
    n_ions = cc_env.size

    GJ_Ith = np.empty ((n_ions, n_GJ))
    GJ_Gth = np.empty ((n_ions, n_GJ))

    # Compute ion drift and diffusion through GJs. Assume fixed GJ spacing
    # of gj_len between connected cells.
    # First, compute d_conc/dx (assume constant conc in cells, and constant
    # gradients in the GJs).
    GJ_from = gj_connects['from']	# Arrays of [n_GJ]
    GJ_to   = gj_connects['to']
    D_scale = gj_connects['scale']

    for ion_index in range(n_ions):
        deltaC_GJ = (cc_cells[ion_index,GJ_to]-cc_cells[ion_index,GJ_from]) \
                      / p.gj_len

        # Assume that ion concentration for any ion is constant within a cell,
        # and then transitions linearly across a GJ. Then c_ave[g] is the conc
        # of the current ion, in the middle of GJ #g. Why do we care? Because
        # when we compute flux = velocity * concentration, then this is the
        # concentration that we will use (regardless of which direction the
        # drift current is actually flowing).
        c_avg = (cc_cells[ion_index,GJ_to] + cc_cells[ion_index,GJ_from]) / 2

        # Finally, electrodiffusive gj flux:
        # f_gj[i] is flux (moles/(m2*s)), in the direction from GJ input to
        # output. Note that D/kT gives the drift mobility.
        D = GJ_diffusion[ion_index] * D_scale	# array[n_GJ], since D_scale is
        alpha = (c_avg * p.q * z_array[ion_index]) * (D/(p.kb*p.T*p.gj_len))
        GJ_Ith[ion_index,:] = -D*deltaC_GJ
        GJ_Gth[ion_index,:] = -alpha

    return (GJ_Ith, GJ_Gth)

# Takes an array of magic_dtype. It is either [N_CELLS] (for ion-channel
# magic) or [N_GJs] (for GJ magic).
# Returns a same-sized array of scalar scale factors in [0,1].
def eval_magic (magic_arr):
    global Vm, cc_cells

    type = magic_arr['type']	# Magic_arr is a structured array; break it
    kM   = magic_arr['kM']	# into a separate simple array for each field.
    N    = magic_arr['N']
    cell = magic_arr['cell']
    ion  = magic_arr['ion']
    cell2= magic_arr['cell2']

    # The default is type==0, which results in scale=1 (i.e., no scaling)
    use_Vmem = np.flatnonzero (type==1) # indices of cells using Vmem
    buf_ion  = np.flatnonzero (type==2) # Hill buffer with ions
    Hill_ion = np.flatnonzero (type>=2) # Hill inv or buf with ions

    # Some advanced-indexing trickery for the Hill-model items (i.e., the
    # magic_arr[] elements whose input is a concentration and not Vm).
    # Create inps[] such that inps[i] is the ion concentration used by
    # magic_arr[i]'s Hill model.
    # Say that cells #3 and #5 are using ion concentrations. Then, using the
    # ion[] and cell[] arrays just above,
    #	magic_arr[3] takes its input from ion number ion[3] (in cell # cell[3])
    #	magic_arr[5] takes its input from ion number ion[5] (in cell # cell[5])
    # and we must set inps[3] = cc_cells[ion[3],cell[3]]
    #		      inps[5] = cc_cells[ion[5],cell[5]]
    inps  = np.empty (magic_arr.size)	# Temporary array to build inputs
    inps[Hill_ion] = cc_cells[ion[Hill_ion], cell[Hill_ion]]

    # Implement a Hill-inverter function -- for those magic_arr items that are
    # either Hill inverters or buffers! Then fix up the buffers by doing 1-x.
    scale = np.ones  (magic_arr.size)	# Final array to return to the caller
    scale[Hill_ion] = 1 / (1 + ((inps[Hill_ion]/kM[Hill_ion])**N[Hill_ion]))
    scale[buf_ion] = 1 - scale[buf_ion]

    # And similar, but even trickier, for items that use Vmem gating (type==1).
    # Implement scale = 1 / 1+exp(N*(v1-v2-kM))
    # Note: "cell = -1" means the cell gets Vm=0
    V1 = np.empty(magic_arr.size)
    V2 = np.empty(magic_arr.size)
    V1[use_Vmem] = Vm[cell [use_Vmem]]
    V2[use_Vmem] = Vm[cell2[use_Vmem]]
    V1_is0V_idx = np.flatnonzero (cell ==-1) # the i such that cell[i]== -1
    V2_is0V_idx = np.flatnonzero (cell2==-1)
    V1[V1_is0V_idx] = 0
    V2[V2_is0V_idx] = 0
    inps [use_Vmem] = V1[use_Vmem] - V2[use_Vmem] - kM[use_Vmem]
    scale[use_Vmem] = 1 / (1 + np.exp(N[use_Vmem]*inps[use_Vmem]))

    return (scale)

def magic_Hill_buf (input_ion, N, kM, input_cell):
    return ((2, kM, N, input_cell, input_ion, 0))
def magic_Hill_inv (input_ion, N, kM, input_cell):
    return ((3, kM, N, input_cell, input_ion, 0))