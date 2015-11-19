Solvers
=======

We will here describe the inheritance hierarchy for generating solvers, in
order to use and extend it properly.  The runtime creation of solver objects
relies on the Factory Method pattern :cite:`Gamma1994,Alexandrescu2001`,
implemented through the generic Factory class.

.. image:: ../gfx/solver.png
   :scale: 70 %
   :align: center

PCMSolver
---------
.. doxygenclass:: PCMSolver
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

IEFSolver
---------
.. doxygenclass:: IEFSolver
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

CPCMSolver
----------
.. doxygenclass:: CPCMSolver
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

VPCMSolver
----------
.. doxygenclass:: VPCMSolver
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

VIEFSolver
----------
.. doxygenclass:: VIEFSolver
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:

VCPCMSolver
-----------
.. doxygenclass:: VCPCMSolver
   :project: PCMSolver
   :members:
   :protected-members:
   :private-members:
