Wavefunction optimization of bulk MnSe
=====

.. _MnSe:

Aim
------------

The aim of this tutorial is to provide users with publication quality examples on calculating total energies of bulk materials using a literature derived high-throughput strategy. 


Prerequisites
----------------

This workflow assumes that you have completed the installation of Quantum Espresso (QE), QMCPACK and Nexus. Otherwise, please refer to `qmcpack.org <https://qmcpack.org>`_ for further instructions on the installation of these software. 
A basic knowledge of Nexus suite (e.g. loading modules, basic elements of Nexus workflows) is also required. Please refer to `nexus-workflows.readthedocs.io <https://nexus-workflows.readthedocs.io/>`_.

Calculations that would need to be performed in QE
^^^^^^^^^^^^^^^^

* Kinetic energy cutoff convergence study on the pseudopotentials used in the materials. Currently (2024), we mainly utilize `ccECP` potentials in most QMC calculations, which require very large (400-700 Ry) kinetic energy cutoff for accurate calculations. 

* K-point convergence study in DFT: We assume that the k-point sampling that converges the DFT total energy will be sufficient to provide similar accuracy in QMC calculations. 

Calculations that would need to be performed in QMCPACK
^^^^^^^^^^^^^^^^

* Hybrid representation calculations: Spline representation of the wavefunction can be demanding on the memory which is often more scarce for GPUs. Currently (2025), most efficient QMCPACK application stores the wavefunction memory untis that are most easily accessible to the GPU units. Therefore, hybrid representation could allow separating core and valence regions of a material to achieve a memory reduction typically around 90 \%.

Calculation steps
----------------

* Convergence tests on bulk MnSe shows that a kinetic energy cutoff of 500 Ry and a kpoint grid of (:math:`8\times8\times4`) is sufficient to achieve a resolution of 1 meV/atom on the DFT total energy. 

* DFT+U scans:
    * A pair of SCF and NSCF calculation at every U value scanned. NSCF calculation should correspond to a reasonable sized supercell (16 atoms in this example).
    * First variance and then mixed energy/reweighed variance minimization (0.95/0.05) using 2-body Jastrows and followed by the mixed cost 3-body jastrow optimization
    * Twist-averaged DMC calculations at :math:`2\times2\times2` grid. 

* Using the optimal U value, extrapolate to the infinite limit:
    * A single self-consistent field calculation of the primitive cell of the material calculated using converged DFT parameters.
    * Non self-consistent field calculations of the primitive cell at different general (non-diagonal) :math:`3\times3` supercell (tiling) vectors, but fixed k-point grid of :math:`2\times2\times2`. :math:`2\times2\times2` kgrid ensures that the calculations are performed on real twists (wavefunction has no imaginary component), thus the "real" version of the qmcpack compiler can also be used in this context. This route also reduces the memory cost of calculations by half. 
    * Jastrow optimization for each supercell size (similar to the DFT+U step)
    * Twist-averaged DMC calculations at :math:`2\times2\times2` grid for each supercell. 

.. .. figure:: ../../prep/BN_workflow.png
..    :alt: Bilayer BN workflow
..    :width: 100%
..    :align: center

..    Schematic of DFT-VMC-DMC calculation workflow for the bilayer binding energy of BN


Contents of working directory
----------------
.. code-block:: text
  
  /MnSe_tutorial
  ├── pseudos/ 
  │   ├── Mn.ccECP-soft.upf 
  │   ├── Mn.ccECP-soft.xml 
  │   ├── Se.ccECP.upf 
  │   └── Se.ccECP.xml 
  ├── README 
  ├── opt_u.py
  ├── run_library.py 
  ├── run_dmc.py 
  └── structures/ 
      └── MnSe.poscar

Complete Nexus scripts
----------------

Wavefunction optimization script (run_u.py)
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: python

  #!/usr/bin/env python

  # user library imports
  from run_library import get_dft_settings, get_qmc_settings
  # nexus imports
  from nexus import run_project, read_structure, obj
  from nexus import generate_physical_system
  from nexus import generate_pwscf
  from nexus import generate_pw2qmcpack
  from nexus import generate_qmcpack
  from structure import optimal_tilematrix

  structure = 'structures/MnSe.poscar'
  tiling_volume                = 4
  primitive_mag_moment         = 10
  tiling_vector, tiling_wigner = optimal_tilematrix(read_structure(structure), volfac=tiling_volume)
  u_values                     = [1e-6, 1, 2, 3, 4, 5, 6, 7,8]

  system_shared = obj(
      Mn        = 15,
      Se        = 6,
  )

  qmc_shared = obj(
      pseudos      = 'Mn.ccECP-soft.xml  Se.ccECP.xml'.split()
  )

  system_prim = generate_physical_system(
      structure = read_structure(structure),
      net_spin  = primitive_mag_moment, 
      **system_shared
  )

  system_tiled = generate_physical_system(
      structure = read_structure(structure),
      tiling    = tiling_vector,
      kgrid     = (2,2,2),
      net_spin  = primitive_mag_moment * tiling_volume,
      **system_shared
  )

  for u in u_values:
      dft_shared = obj(
          kgrid     = (8,8,4),
          ecutwfc   = 500,
          pseudos   = 'Mn.ccECP-soft.upf Se.ccECP.upf'.split(),
          start_mag = obj(Mn=1),
          hubbard   = {'U' : {'Mn-3d':u}},
      )

      scf_shared, nscf_shared, conv_shared = get_dft_settings(**dft_shared)

      scf_path = 'scf_u_{}'.format(u)
      scf_run = generate_pwscf(
              system  = system_prim,
              path    = scf_path,
              
              **scf_shared
          )
      nscf_path = 'nscf_u_{}_v_{}'.format(u, tiling_volume)
      nscf_run = generate_pwscf(
          system = system_tiled,
          path = nscf_path,
          dependencies = (scf_run, 'charge_density'),
          **nscf_shared
      )        
      conv_run = generate_pw2qmcpack(
          path         = nscf_path,   
          dependencies = (nscf_run, 'orbitals'),
          **conv_shared
      )        

      dmc_path = 'dmc_u_{}_v_{}'.format(u, tiling_volume)
          
      # Optimize jastrows using the first u value listed in U-values
      if u == u_values[0]:
          j2_path = 'j2_u_{}_v_{}'.format(u, tiling_volume)
          j3_path = 'j3_u_{}_v_{}'.format(u, tiling_volume)
          j2_shared, j3_shared, dmc_shared  = get_qmc_settings(system = system_tiled, **qmc_shared)
          
          j2_run = generate_qmcpack(path = j2_path,
                                    dependencies = (conv_run, 'orbitals'),
                                    **j2_shared)

          j3_run = generate_qmcpack(path = j3_path,
                                    dependencies = (conv_run, 'orbitals'),
                                    **j3_shared)
          
      else:
          _, _, dmc_shared = get_qmc_settings(system = system_tiled, **qmc_shared)

      dmc_run = generate_qmcpack(path = dmc_path,
                                  dependencies = [(j3_run, 'jastrow'),(conv_run, 'orbitals')],
                                  **dmc_shared)
  run_project()

DMC script (run_dmc.py)
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: python

  #!/usr/bin/env python

  # user library imports
  from run_library import get_dft_settings, get_qmc_settings
  # nexus imports
  from nexus import run_project, read_structure, obj
  from nexus import generate_physical_system
  from nexus import generate_pwscf
  from nexus import generate_pw2qmcpack
  from nexus import generate_qmcpack
  from structure import optimal_tilematrix

  structure = 'structures/MnSe.poscar'
  tiling_volumes               = [8, 12, 16]
  primitive_mag_moment         = 10
  u                            = 3
  system_shared = obj(
      Mn        = 15,
      Se        = 6,
  )

  dft_shared = obj(
      kgrid    = (8,8,4),
      ecutwfc  = 500,
      pseudos  = 'Mn.ccECP-soft.upf Se.ccECP.upf'.split(),
      start_mag = obj(Mn=1)
  )

  qmc_shared = obj(
      pseudos      = 'Mn.ccECP-soft.xml  Se.ccECP.xml'.split()
  )

  scf_shared, nscf_shared, conv_shared = get_dft_settings(**dft_shared)

  system_prim = generate_physical_system(
      structure = read_structure(structure),
      net_spin  = primitive_mag_moment, 
      **system_shared
  )

  scf_path = 'scf_u_{}'.format(u)
  scf_run = generate_pwscf(
          system = system_prim,
          path = scf_path,
          **scf_shared
      )

  for v in tiling_volumes:
      tiling_vector, tiling_wigner = optimal_tilematrix(read_structure(structure), volfac=v)

      system_tiled = generate_physical_system(
          structure = read_structure(structure),
          tiling    = tiling_vector,
          kgrid     = (2,2,2),
          net_spin  = primitive_mag_moment * v,
          **system_shared
      )

      nscf_path = 'nscf_u_{}_v_{}'.format(u, v)
      
      nscf_run = generate_pwscf(
          system = system_tiled,
          path = nscf_path,
          dependencies = (scf_run, 'charge_density'),
          **nscf_shared
      )        

      conv_run = generate_pw2qmcpack(
          path         = nscf_path,   
          dependencies = (nscf_run, 'orbitals'),
          **conv_shared
      )        

      dmc_path = 'dmc_u_{}_v_{}'.format(u, v)        
      j2_path = 'j2_u_{}_v_{}'.format(u, v)
      j3_path = 'j3_u_{}_v_{}'.format(u, v)
      j2_shared, j3_shared, dmc_shared  = get_qmc_settings(system = system_tiled, **qmc_shared)
      
      j2_run = generate_qmcpack(path = j2_path,
                                  dependencies = (conv_run, 'orbitals'),
                                  **j2_shared)

      j3_run = generate_qmcpack(path = j3_path,
                                  dependencies = (conv_run, 'orbitals'),
                                  **j3_shared)
      
      dmc_run = generate_qmcpack(path = dmc_path,
                                  dependencies = [(j3_run, 'jastrow'),(conv_run, 'orbitals')],
                                  **dmc_shared)
  run_project()


Workflow library script (run_library.py)
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Here, we reuse the same user library script used in the :ref:`hBN tutorial <hBN_wf_script>`.

.. code-block:: python

  #!/usr/bin/env python
  # nexus imports
  from nexus import Job, obj
  from nexus import settings
  from nexus import linear, loop, vmc, dmc
  from qmcpack_input import spindensity
  # general settings for nexus
  settings(
      pseudo_dir    = './pseudos',
      status_only   = 0,                    # only show status of runs
      generate_only = 0,                    # only make input files
      sleep         = 3.0,                    # check on runs every 3 secondsa
      machine       = 'ws16'                # local machine is 16 core workstation
      )

  def get_dft_settings(kgrid     = None, 
                       ecutwfc   = None,
                       pseudos   = None,
                       start_mag = None, 
                       hubbard   = None,):
      
      if settings.machine == 'ws16':
          dft_job = Job(cores=16, app='/pw.x')
          conv_job = Job(cores=1, app='pw2qmcpack.x')
      else:
          print('Error: Unknown computer for DFT, using {}'.format(settings.machine))
          exit()

      qe_shared = obj(
          job          = dft_job,
          input_type   = 'generic',
          ecutwfc      = ecutwfc,                             # DFT planewave energy cutoff
          input_DFT    = 'PBE',                               # DFT functional
          conv_thr     = 1e-8,                                # SCF convergence threshold
          wf_collect   = True,                                # write orbitals
          pseudos      = pseudos,                             # QE Pseudopotentials
          start_mag    = start_mag                            # Starting magnetization
          hubbard      = hubbard                              # Hubbard-U parameters (QE ver. > 7.2)
      )

      scf_shared = obj(
          nosym        = False,            # use symmetry
          identifier   = 'scf',            # identifier/file prefix
          calculation  = 'scf',            # perform scf calculation
          kgrid        = kgrid,        # Converged DFT k-grid
          occupations  = 'smearing',       # Occupation scheme
          smearing     = 'gauss',          # Smearing type
          degauss      = 0.001,            # Smearing width
          **qe_shared
      )

      nscf_shared = obj(
          nosym        = True,             # don't use symmetry
          identifier   = 'nscf',           
          calculation  = 'nscf',           # perform nscf calculation
          occupations  = 'fixed',          # Use fixed occupations for insulators
          **qe_shared
      )    

      conv_shared = obj(
              identifier   = 'conv',          
              job          = conv_job,
              write_psir   = False,
      )

      return scf_shared, nscf_shared, conv_shared

  def get_qmc_settings(system      = None,
                       hybrid_rcut = None,
                       hybrid_lmax = None, 
                       meshfactor  = None,
                       pseudos     = None):
      
      num_kpt = len(system.structure.kpoints)
      
      if settings.machine == 'ws16':
          num_cores   = 16
          dmc_threads = int(num_cores/num_kpt)
          assert dmc_threads > 0, "Number of processors ({}) should be larger than or equal to the number of kpoints ({})".format(num_cores, num_kpt)
          opt_job = Job(cores = num_cores, threads = num_cores, app='qmcpack_complex')
          dmc_job = Job(cores = num_cores, threads = dmc_threads, app='qmcpack_complex')
      else:
          print('Error: Unknown computer for QMC, using {}'.format(settings.machine))
          exit()

      system.structure.change_units('B')
      rwigner = system.structure.rwigner()

      qmc_settings = obj(
          system          = system, 
          input_type      = 'basic',
          pseudos         = pseudos,
          driver          = 'batched',
          hybrid_rcut     = hybrid_rcut,
          hybrid_lmax     = hybrid_lmax,
          meshfactor      = meshfactor,
          lr_handler      = 'ewald',
          lr_dim_cutoff   = 30,
          spin_polarized  = True,
      )

      opt_parameters = obj(
          num_varmin_j2   = 12,
          num_emin_j2     = 8,
          num_emin_j3     = 6,
          j2_init         = "rpa",
          num_j1_jastrows = 10,
          num_j2_jastrows = 10,
          num_j3_jastrows = 3,
          j3_rcut         = 4.0 if rwigner > 4.0 else rwigner,
          timestep        = 1.0
      )

      opt_settings = obj(
          job             = opt_job,
          twistnum        = 0,
          # warmupsteps     = 200,
          # samples         = 128000,
          # blocks          = 100,
          # steps           = 1,
          # timestep        = 1.0,
          # substeps        = 10,
      )
      opt_settings = opt_settings.set(qmc_settings)


      varmin = linear(
          energy               = 0.0, 
          unreweightedvariance = 1.0,
          reweightedvariance   = 0.0, 
          minwalkers           = 1e-4,
          shift_i              = 0.05,
          shift_s              = 1.0,
          warmupsteps          = 200,
          blocks               = 100,
          steps                = 1,
          timestep             = 1.0,
          minmethod            = "OneShiftOnly",
          substeps             = 10,        
      )    

      emin = varmin.copy()
      emin.minwalkers             = 0.5
      emin.energy                 = 0.95
      emin.unreweightedvariance   = 0.0
      emin.reweightedvariance     = 0.05
      emin.shift_i                = 0.01
      emin.shift_s                = 1.0

      
      j2_settings     = obj(
          calculations = [loop(max=opt_parameters.num_varmin_j2, qmc=varmin), 
                          loop(max=opt_parameters.num_emin_j2,   qmc=emin)],
          jastrows     = [('J1','bspline',opt_parameters.num_j1_jastrows, rwigner),                        # 1 body bspline jastrow
                          ('J2','bspline',opt_parameters.num_j2_jastrows, 'init', opt_parameters.j2_init)], # 2 body bspline jastrow
          **opt_settings
      )

      j3_settings     = obj(
          calculations = [loop(max=opt_parameters.num_emin_j3, qmc=emin)],
          jastrows     = [('J3', 'polynomial', opt_parameters.num_j3_jastrows,3, opt_parameters.j3_rcut)],
          **opt_settings
      )    

      dmc_parameters = obj(
          vmcdt       = 0.3,
          vmcwarmup   = 25,
          vmcblocks   = 100,
          vmcsubsteps = 4,
          dmc_eq_dt   = 0.02,
          dmc_eq_blocks = 100,
          dmcdt       = 0.005,
          dmcblocks   = 500,
          dmcwarmup   = 100,
          dmcsteps    = 10,
          vmc_walkers_per_rank = 240,
          dmc_walkers_per_rank = 240,
          nonlocalmoves = False, 
      )

      vmc_dmc = obj(
          warmupsteps = dmc_parameters.vmcwarmup,
          blocks      = dmc_parameters.vmcblocks,
          steps       = 1,
          timestep    = dmc_parameters.vmcdt,
          substeps    = dmc_parameters.vmcsubsteps,
          walkers_per_rank = dmc_parameters.vmc_walkers_per_rank
      )
      dmc_eq  = obj(
          warmupsteps = dmc_parameters.dmcwarmup,
          blocks      = dmc_parameters.dmc_eq_blocks,
          steps       = dmc_parameters.dmcsteps,
          timestep    = dmc_parameters.dmc_eq_dt,
          walkers_per_rank = dmc_parameters.dmc_walkers_per_rank,
          nonlocalmoves = dmc_parameters.nonlocalmoves, 
      )
      dmc_stat = obj(
          warmupsteps = dmc_parameters.dmcwarmup,
          blocks      = dmc_parameters.dmcblocks,
          steps       = dmc_parameters.dmcsteps,
          timestep    = dmc_parameters.dmcdt,
          walkers_per_rank = dmc_parameters.dmc_walkers_per_rank,
          nonlocalmoves = dmc_parameters.nonlocalmoves, 
      )
      
      dmc_settings = obj(
          job           = dmc_job,
          calculations  = [vmc(**vmc_dmc), dmc(**dmc_eq), dmc(**dmc_stat)],
          estimators    = [spindensity(dr=3*[0.3])],
          **qmc_settings    
      )

      return j2_settings, j3_settings, dmc_settings


Work through of the Nexus scripts
----------------

The workflows in this example are managed by :code:`opt_u.py`, :code:`run_dmc.py`, while DFT and QMC settings are generated using functions imported from :code:`run_library.py`. 
Therefore, all scripts need to be in the same directory to complete the workflow. 
If you plan to use modified versions of the scripts in your own work repeatedly, you can alternatively place :code:`run_library.py` 
in a directory defined under :code:`PYTHONPATH` environment variable to make it accessible to Python interpreter. 

This workflow can provide a good starting point to write a general workflow that can be used to calculate DMC ground state energies of bulk materials. 
The workflow in this example is based on the workflow implemented in :cite:`Saritas2017` which aims to calculate the formation energies of uncorrelated solids using DMC. 
Basic structure of the scripts, especially the division of the workflow and the settings to separate files is very similar to the workflow explained in the :ref:`hBN tutorial <hBN>`, therefore will not be covered here again. 

[CONTINUE]





.. bibliography::


.. To retrieve a list of random ingredients,
.. you can use the ``lumache.get_random_ingredients()`` function:

.. .. autofunction:: lumache.get_random_ingredients

.. The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
.. or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
.. will raise an exception.

.. .. autoexception:: lumache.InvalidKindError

.. For example:

.. >>> import lumache
.. >>> lumache.get_random_ingredients()
.. ['shells', 'gorgonzola', 'parsley']

