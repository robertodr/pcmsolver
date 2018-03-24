Time-dependent solvers
======================

We will here describe the inheritance hierarchy for generating time-dependent solvers, in
order to use and extend it properly.  The runtime creation of solver objects
relies on the Factory Method pattern :cite:`Gamma1994,Alexandrescu2001`,
implemented through the generic Factory class.

.. image:: ../gfx/td_solver.png
   :scale: 70 %
   :align: center

ITDSolver
---------
.. doxygenclass:: pcm::ITDSolver
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

TDIEFSolver
-----------
.. doxygenclass:: pcm::td_solver::TDIEFSolver
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

TDCPCMSolver
------------
.. doxygenclass:: pcm::td_solver::TDCPCMSolver
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

TDSingleIEFSolver
-----------------
.. doxygenclass:: pcm::td_solver::TDSingleIEFSolver
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:
