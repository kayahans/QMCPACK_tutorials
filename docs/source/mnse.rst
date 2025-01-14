Wavefunction optimization of bulk MnSe
=====

.. _MnSe:

Aim
------------

The aim of this tutorial is to provide users with publication quality examples on calculating binding energies of layered materials. 

.. .. code-block:: console

..    (.venv) $ pip install lumache

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

* Convergence tests on bulk MnSe shows that a kinetic energy cutoff of 700 eV and a kpoint grid of (:math:`8\times8\times8`) [TODO] is sufficient to achieve a resolution of 1 meV/atom on the DFT total energy. 

* A single self-consistent field calculation of the primitive cell of the material calculated using converged DFT parameters.

* Non self-consistent field calculations of the primitive cell at different general :math:`3\times3` supercell (tiling) vectors, but fixed k-point grid of :math:`2\times2\times2`. 

* Jastrow optimization for each supercell size

    * Variance minimization using 2-body Jastrows

    * Mixed energy/reweighed variance minimization (0.95/0.05) using 2-body Jastrows

    * Mixed energy/reweighed variance minimization (0.95/0.05) using 3-body Jastrows


* Twist-averaged DMC calculations at :math:`2\times2\times2` grid. 

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
  │   ├── Mn.ccECP.upf 
  │   ├── Mn.ccECP.xml 
  │   ├── Se.ccECP.upf 
  │   └── Se.ccECP.xml 
  ├── README 
  ├── opt_u.py
  ├── run_library.py 
  ├── run_dmc.py 
  └── structures/ 
      └── MnSe.xsf 

Complete Nexus scripts
----------------

Wavefunction optimization script (run_u.py)
^^^^^^^^^^^^^^^^^^^^^^^^^^^
DMC script (run_dmc.py)
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Workflow library script (run_library.py)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Work through of the Nexus scripts
----------------

The workflows in this example are managed by :code:`opt_u.py`, :code:`run_dmc.py`, while DFT and QMC settings are generated using functions imported from :code:`run_library.py`. 
Therefore, all scripts need to be in the same directory to complete the workflow. 
If you plan to use modified versions of the scripts in your own work repeatedly, you can alternatively place :code:`run_library.py` 
in a directory defined under :code:`PYTHONPATH` environment variable to make it accessible to Python interpreter. 

This workflow can provide a good starting point to write a general workflow that can be used to calculate DMC ground state energies of bulk materials. 
The workflow in this example is based on the workflow implemented in :cite:`Saritas2017` which aims to calculate the formation energies of uncorrelated solids using DMC. 
Basic structure of the scripts, especially the division of the workflow and the settings to separate files is very similar to the workflow explained in the :ref:`hBN tutorial <hBN>`, therefore will not be covered here again. 






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

