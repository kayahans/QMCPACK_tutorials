Binding energy of bilayer hexagonal BN 
=================

.. _hBN:

Aim
------------

The aim of this tutorial is to provide users with publication quality examples on calculating binding energies of layered materials. 

Prerequisites
----------------
This workflow assumes that you have completed the installation of Quantum Espresso (QE), QMCPACK and Nexus. Otherwise, please refer to `qmcpack.org <https://qmcpack.org>`_ for further instructions on the installation of these software. 
A basic knowledge of Nexus suite (e.g. loading modules, basic elements of Nexus workflows) is also required. Please refer to `nexus-workflows.readthedocs.io <https://nexus-workflows.readthedocs.io/>`_.

Calculations that would need to be performed in QE
^^^^^^^^^^^^^^^^

* Kinetic energy cutoff convergence study on the pseudopotentials used in the materials. Currently (2025), we mainly utilize `ccECP` potentials from `pseudopotential-library <https://pseudopotentiallibrary.org>`_. These pseudopotentials typically require very large (400-700 Ry) kinetic energy cutoff for accurate calculations. 

* K-point convergence study in DFT: We assume that the k-point sampling that converges the DFT total energy will be sufficient to provide similar accuracy in QMC calculations. 

Calculations that would need to be performed in QMCPACK
^^^^^^^^^^^^^^^^

* Hybrid representation calculations: Spline representation of the wavefunction can be demanding on the memory which is often more scarce for GPUs. Currently (2024), most efficient QMCPACK application stores the wavefunction memory untis that are most easily accessible to the GPU units. Therefore, hybrid representation could allow separating core and valence regions of a material to achieve a memory reduction typically around 90 \%.

Calculation steps
----------------

* Convergence tests on BN shows that a kinetic energy cutoff of 400 eV and a kpoint grid of (8x8x1) is sufficient to achieve a resolution of 1 meV/atom on the DFT total energy. 

* Self-consistent field calculation of the primitive cell of the material calculated using converged DFT parameters as above for interlayer separations of 2.5, 3.0, 3.25, 3.5, 4.0, 4.5 and 5.0 Angstroms. 

* Non self-consistent field calculations of the primitive cell for each interlayer separation, which then will be unfolded to supercells (V = 4, 9, 16; M = (2x2x1), (3x3x1), (4x4x1)) for finite size extrapolation in QMCPACK. Here, V is the volume factor compared to primitive supercell, while M is the diagonal supercell matrix used to tile the primitive cell.  Different choices can be explored here such that non-diagonal supercell matrices (M) can also be utilized.  However, for 2D materials we expect that the coulomb interaction screened more weakly compared to 3D materials, hence using diagonal supercells could provide better cancellation of errors reducing the finite size effects. 


* Jastrow optimization for each supercell size at a single interlayer separation

    * Variance minimization using 2-body Jastrows

    * Mixed energy/reweighed variance minimization (0.95/0.05) using 2-body Jastrows

    * Mixed energy/reweighed variance minimization (0.95/0.05) using 3-body Jastrows


* Twist-averaged DMC calculations

.. figure:: ../../prep/BN_workflow.png
   :alt: Bilayer BN workflow
   :width: 100%
   :align: center

   Schematic of DFT-VMC-DMC calculation workflow for the bilayer binding energy of BN


Contents of working directory
----------------
.. code-block:: text
  
  /BN_tutorial
  ├── pseudos/ 
  │   ├── B.ccECP.upf 
  │   ├── B.ccECP.xml 
  │   ├── N.ccECP.upf 
  │   └── N.ccECP.xml 
  ├── README 
  ├── run_library.py 
  ├── run.py 
  └── structures/ 
      ├── hBN_d_2500.xsf 
      ├── hBN_d_3000.xsf 
      ├── hBN_d_3250.xsf 
      ├── hBN_d_3500.xsf 
      ├── hBN_d_4000.xsf 
      ├── hBN_d_4500.xsf 
      ├── hBN_d_5000.xsf 
      └── hBN_mono.xsf 

Complete Nexus scripts
----------------

Workflow running script (run.py)
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

  # material specific DFT settings start
  structures = {2.5:  'structures/hBN_d_2500.xsf',
                3.0:  'structures/hBN_d_3000.xsf',
                3.25: 'structures/hBN_d_3250.xsf',
                3.5:  'structures/hBN_d_3500.xsf',
                4.0:  'structures/hBN_d_4000.xsf',
                4.5:  'structures/hBN_d_4500.xsf',
                5.0:  'structures/hBN_d_5000.xsf',
                'mono' : 'structures/hBN_mono.xsf'}
  interlayer_separations  = list(structures.keys())
  tiling_vectors          = [(2,2,1), (3,3,1), (4,4,1)]
  tiling_kgrids           = {(2,2,1):(4,4,1), 
                             (3,3,1):(2,2,1), 
                             (4,4,1):(2,2,1)}

  system_shared = obj(
              B        = 3,
              N        = 5,
              net_spin = 0
  )

  dft_shared = obj(
                  kgrid    = (8,8,1),
                  ecutwfc  = 400,
                  pseudos  = 'B.ccECP.upf N.ccECP.upf'.split()
  )
  qmc_shared = obj(
                  hybrid_rcut  = obj(B=1.1, N=1.1),
                  hybrid_lmax  = obj(B=5, N=5),
                  meshfactor   = 0.5,
                  pseudos      = 'B.ccECP.xml  N.ccECP.xml'.split()
  )
  # material specific DFT settings end 

  scf_shared, nscf_shared, conv_shared = get_dft_settings(**dft_shared) # Load DFT settings 
  
  # Binding energy workflow start 
  for d in interlayer_separations:
      if isinstance(d, (int, float)):
          d_name = int(d*1000)
      else:
          d_name = d

      scf_path = 'scf_{}'.format(d_name)
      prim_system = generate_physical_system(
              structure = read_structure(structures[d]),
              **system_shared
          )
      scf_run = generate_pwscf(
              system = prim_system,
              path = scf_path,
              **scf_shared
          )
      for t in tiling_vectors:
          nscf_path = 'nscf_{}_{}'.format(d_name, t[0])
          tiled_system = generate_physical_system(
              structure = read_structure(structures[d]),
              tiling   = t,
              kgrid    = tiling_kgrids[t],
              **system_shared
          )

          nscf_run = generate_pwscf(
              system = tiled_system,
              path = nscf_path,
              **nscf_shared
          )        

          conv_run = generate_pw2qmcpack(
              path         = nscf_path,   
              dependencies = (nscf_run, 'orbitals'),
              **conv_shared
          )        

          dmc_path = 'dmc_{}_{}'.format(d_name, t[0])
          
          # Optimize jastrows using the first structure listed in 
          if d == interlayer_separations[0]: # In this example, this is d == 2.5 since dictionary keys are always ordered in Python 3.7+
              j2_path = 'j2_{}_{}'.format(d_name, t[0])
              j3_path = 'j3_{}_{}'.format(d_name, t[0])
              j2_shared, j3_shared, dmc_shared  = get_qmc_settings(system = tiled_system, **qmc_shared) # Load QMC settings 
              
              j2_run = generate_qmcpack(path = j2_path,
                                        dependencies = (conv_run, 'orbitals'),
                                        **j2_shared)

              j3_run = generate_qmcpack(path = j3_path,
                                        dependencies = (conv_run, 'orbitals'),
                                        **j3_shared)
          else:
              _, _, dmc_shared = get_qmc_settings(system = tiled_system, **qmc_shared) # Load DMC settings only (no J2/J3)

          dmc_run = generate_qmcpack(path = dmc_path,
                                      dependencies = [(j3_run, 'jastrow'),(conv_run, 'orbitals')],
                                      **dmc_shared)
    # Binding energy workflow end  
    run_project()

Workflow library script (run_library.py)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

  #!/usr/bin/env python

  # nexus imports
  from nexus import Job, obj
  from nexus import settings
  from nexus import linear, loop, vmc, dmc
  from qmcpack_input import spindensity

  structures = {2.5:  'structures/hBN_d_2500.xsf',
                3.0:  'structures/hBN_d_3000.xsf',
                3.25: 'structures/hBN_d_3250.xsf',
                3.5:  'structures/hBN_d_3500.xsf',
                4.0:  'structures/hBN_d_4000.xsf',
                4.5:  'structures/hBN_d_4500.xsf',
                5.0:  'structures/hBN_d_5000.xsf',
                'mono' : 'structures/hBN_mono.xsf'}
  interlayer_separations  = list(structures.keys())
  tiling_vectors          = [(2,2,1), (3,3,1), (4,4,1)]
  tiling_kgrids           = {(2,2,1):(4,4,1), 
                             (3,3,1):(2,2,1), 
                             (4,4,1):(2,2,1)}

  system_shared = obj(
              B        = 3,
              N        = 5,
              net_spin = 0
  )

  dft_shared = obj(
                  kgrid    = (8,8,1),
                  ecutwfc  = 400,
                  pseudos  = 'B.ccECP.upf N.ccECP.upf'.split()
  )
  qmc_shared = obj(
                  hybrid_rcut  = obj(B=1.1, N=1.1),
                  hybrid_lmax  = obj(B=5, N=5),
                  meshfactor   = 0.5,
                  pseudos      = 'B.ccECP.xml  N.ccECP.xml'.split()
  )

  scf_shared, nscf_shared, conv_shared = get_dft_settings(**dft_shared)

  for d in interlayer_separations:
      if isinstance(d, (int, float)):
          d_name = int(d*1000)
      else:
          d_name = d

      scf_path = 'scf_{}'.format(d_name)
      prim_system = generate_physical_system(
              structure = read_structure(structures[d]),
              **system_shared
          )
      scf_run = generate_pwscf(
              system = prim_system,
              path = scf_path,
              **scf_shared
          )
      for t in tiling_vectors:
          nscf_path = 'nscf_{}_{}'.format(d_name, t[0])
          tiled_system = generate_physical_system(
              structure = read_structure(structures[d]),
              tiling   = t,
              kgrid    = tiling_kgrids[t],
              **system_shared
          )

          nscf_run = generate_pwscf(
              system = tiled_system,
              path = nscf_path,
              **nscf_shared
          )        

          conv_run = generate_pw2qmcpack(
              path         = nscf_path,   
              dependencies = (nscf_run, 'orbitals'),
              **conv_shared
          )        

          dmc_path = 'dmc_{}_{}'.format(d_name, t[0])
          
          # Optimize jastrows using the first structure listed in 
          if d == interlayer_separations[0]:
              j2_path = 'j2_{}_{}'.format(d_name, t[0])
              j3_path = 'j3_{}_{}'.format(d_name, t[0])
              j2_shared, j3_shared, dmc_shared  = get_qmc_settings(system = tiled_system, **qmc_shared)
              
              j2_run = generate_qmcpack(path = j2_path,
                                        dependencies = (conv_run, 'orbitals'),
                                        **j2_shared)

              j3_run = generate_qmcpack(path = j3_path,
                                        dependencies = (conv_run, 'orbitals'),
                                        **j3_shared)
              
          else:
              _, _, dmc_shared = get_qmc_settings(system = tiled_system, **qmc_shared)

          dmc_run = generate_qmcpack(path = dmc_path,
                                      dependencies = [(j3_run, 'jastrow'),(conv_run, 'orbitals')],
                                      **dmc_shared)
    run_project()
    
Work through of the Nexus scripts
----------------

The workflow in this example is managed by `run.py`, while DFT and QMC settings are generated using functions imported from `run_library.py`. Therefore, both scripts need to be in the same directory to complete the workflow. 
If you plan to use modified versions of the scripts in your own work repeatedly, you can alternatively place `run_library.py` in a directory defined under :code:`PYTHONPATH` environment variable to make it accessible to Python interpreter. 

We start with `run.py`. This script has 4 main components: 
(1) Importing Nexus and user library modules (from :code:`run_library`), 
(2) materials specific DFT settings and inputs (e.g. structure files, kgrid, kinetic energy cutoffs), 
(3) A :code:`for` loop running over the :code:`interlayer_separations` and 
(4) an inner :code:`for` loop running over the :code:`tiling_vectors`. 

After the Nexus module imports, in the DFT settings section of the script, all the structures to be used are defined in the :code:`structures` dictionary.
Some of the QMC settings that are shared across all calculations are also stored in here, such as hybrid representation parameters and DMC pseudopotential files. 
For more info on the hybrid representation and how to choose hybrid representation parameters, please refer to `QMCPACK manual <https://qmcpack.readthedocs.io/en/develop/intro_wavefunction.html#hybrid-orbital-representation>`_.
In this example, we find that DFT ground state energy of a bilayer system converges at a kgrid of :math:`8\times8\times1` and kinetic energy cutoff of 400 eV.
We use :code:`tiling_vectors` which are simply diagonal supercell matrices of :math:`2\times2\times1`, :math:`3\times3\times1` and :math:`4\times4\times1` that are represented as tuples. 
To ensure convergence in DMC, we use combinations of :code:`tiling_vectors` and :code:`tiling_kgrids`, such that the product of these two matrices exceed or match the converged kgrid in the DFT calculation.
This assumes that the one-body errors in DMC are very similar to DFT, which is a reasonable approximation if the valence band structure is not expected to modified significantly.
The keys of this :code:`structures` dictionary are later used to denote interlayer separations of the bilayer hBN. These xsf files are processed by Nexus in :code:`read_structure(structures[d])` line. 
:code:`generate_pwscf` function creates the :code:`scf_run` object using Quantum Espresso code which is of type :code:`Simulation`. 
Each initialization of a :code:`Simulation` object, such as :code:`scf_run`, enables storing that object in the internal memory of Nexus that is inaccessible to the user, 
therefore these objects do not need to be stored in an array to be passed to the :code:`run_project` function at the end of the script. 
However, storing this object in the local memory as :code:`scf_run`, enables setting up dependencies between these :code:`Simulation` objects, such as :code:`'charge_density'`, :code:`orbitals` and :code:`jastrow`. 
The :code:`tiling_vectors` loop then initiates the NSCF calculations with separate tiling vectors and kgrids, which are later unfolded by QMCPACK and later used as supercell calculations. 
Each tiling of the primitive cell requires defining a separate :code:`PhysicalSystem` object using :code:`generate_physical_system`. 
:code:`PhysicalSystem` object carries the information on both the primitive and the supercell (tiling) of the material studied. 
:code:`nscf_run` and :code:`conv_run` objects define the NSCF calculations and Quantum Espresso to QMCPACK conversion calculations respectively. 
After these steps, in the last if/else block Jastrow optimizations are performed and DMC calculation settings are imported. 
Jastrow optimizations are performed only for a single bilayer separation (:code:`interlayer_separations[0]`, which is 2.5 :math:`\AA`). 
For well optimized wavefunctions (e.g. variance/energy ratios close to or below 0.01 as a rule of thumb), the localization error is minimal, 
therefore Jastrow parameters can be accurately reused across similar geometries (e.g. bilayer separations, equation of states) without modification if the cutoff radii of the Jastrow parameters are smaller than the Wigner-Seitz radii of these geometries. 
We optimize Jastrow parameters from scratch at every supercell size (:code:`tiling_vectors`), because at increased system size, larger Jastrow cutoff radii are allowed, and especially two-body Jastrow parameters with larger cutoff radii could improve the quality of the optimized wavefunction. 
First, two-body jastrow parameters are optimized and then these optimal parameters are used to initialize combined two and three body Jastrow parameters in :code:`j3_run`.
Finally, using converted the Slater wavefunction from :code:`conv_run` and optimized Jastrow parameters from :code:`j3_run`, DMC calculations are executed. 

Next, we work through the `run_library.py` script. This script has 3 main components: 
(1) A call to :code:`settings` function in Nexus,
(2) :code:`get_dft_settings` function,
(3) :code:`get_qmc_settings` function.

The call to :code:`settings` function initializes the workspace for Nexus and controls the workflow execution. :code:`status_only` and :code:`generate_only` can be set to either 0 or 1 to check the workflow status without updating and modifying any simulation object and produce the input files without executing the workflow respectively.
:code:`machine` defines which machine is currently utilized to run the simulations. :code:`ws16` is a generic 16 core workstation.
However, most US-based and some Europe and Japan based supercomputers, and in some cases smaller institutional clusters are also available and kept up-to-date under this setting. 
With :code:`machine` setting, Nexus can become aware of the computer architecture and identify how to interact with its job scheduler to submit jobs and check the status of running simulations. 
:code:`get_dft_settings` is a user defined function which only uses converged kpoint grid, kinetic energy cutoff and pseudopotential files as a minimal setup to define SCF, 
NSCF and Quantum Espresso to QMCPACK conversion run (:code:`scf_shared`, :code:`nscf_shared`, :code:`conv_shared`) settings.  
:code:`get_qmc_settings` is also a user defined function which uses the :code:`PhysicalSystem` object, hybrid-representation parameters and pseudopotential files as a minimal setup to define the settings for QMC calculations. 
Comparing DFT and QMC settings, we can see that DFT settings are oblivious to the material used, while QMC settings require :code:`PhysicalSystem` object as an input parameter, because some parameters such as Jastrow radii are dependent on the system size in our workflow. 
Although these functions aim to be general to large classes of materials by construction, they are far from being truly general and make assumptions on what kind of materials that can be studied using them. 
For example, if the DFT electronic structure is metallic, then DFT settings would need to be modified since :code:`occupations  = 'fixed'` would be invalid and more care would be required for the smearing in the SCF calculation. 
Similarly, if the user wants to achieve certain final uncertainty in the DMC total energies, then the statistical accumulation parameters (e.g. number of blocks, walkers) can be optimized with respect to several factors such as system size, complexity and electron count to minimize computational cost, after benchmarking relevant systems. 


