Binding energy of bilayer hexagonal BN 
=================

.. _aim:

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

* Convergence tests on BN shows that a kinetic energy cutoff of XX eV and a kpoint grid of (8x8x1) is sufficient to achieve a resolution of 1 meV/atom on the DFT total energy. 

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
  
  \BN_tutorial
  ├── pseudos/ 
  │   ├── B.ccECP.upf 
  │   ├── B.ccECP.xml 
  │   ├── N.ccECP.upf 
  │   └── N.ccECP.xml 
  ├── README 
  ├── run_library.py 
  ├── run.py 
  ├── structures/ 
  │   ├── hBN_d_2500.xsf 
  │   ├── hBN_d_3000.xsf 
  │   ├── hBN_d_3250.xsf 
  │   ├── hBN_d_3500.xsf 
  │   ├── hBN_d_4000.xsf 
  │   ├── hBN_d_4500.xsf 
  │   ├── hBN_d_5000.xsf 
  │   └── hBN_mono.xsf 

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

  # material specific settings start
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
  # material specific settings end 

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
    
Breakdown of the Nexus scripts
----------------

The workflow in this example is managed by `run.py`, while most settings that remain invariant between different materials/structures are stored in `run_library.py`. 
This tutorial assumes a familiarity with Nexus workflow scripts, therefore please refer to `https://nexus-workflows.readthedocs.io/en/latest/ <https://nexus-workflows.readthedocs.io/en/latest/>`_ for a general basic understanding of Nexus suite. 
In ``
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

