.. vLab documentation master file, created by
   sphinx-quickstart on Sun Apr 17 21:45:36 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to vLab's documentation!
================================


.. include:: README.md
   :parser: myst_parser.sphinx_

Functionality
==================
.. toctree::
   :maxdepth: 2
   :caption: Lectures

   TrainingMaterial/process_modeling.md

.. toctree::
   :maxdepth: 2
   :caption: Tutorial

   Tutorial/tutorial

.. toctree::
   :maxdepth: 2
   :caption: Integrated Bioprocess

   IntegratedBioprocess/Bioreactor
   IntegratedBioprocess/Chromatography
   IntegratedBioprocess/HarvestTank
   IntegratedBioprocess/ODESolver
   IntegratedBioprocess/Plantwise
   IntegratedBioprocess/PlantwiseSimulator
   IntegratedBioprocess/Util

.. toctree::
   :maxdepth: 2
   :caption: Steady-State Glycosylation

   PerfusionSimulator/PerfusionSimulator

.. toctree::
   :maxdepth: 2
   :caption: Dynamic Glycosylation Simulator

   DynamicGlycosylationSimulator/DynamicGlycosylationSimulator

.. toctree::
   :maxdepth: 2
   :caption: Glycosylation Model Base Class

   GlycosylationModelBase/GlycosylationModelBase

.. toctree::
   :maxdepth: 2
   :caption: Raman Spectropecty

   RamanAnalytics/RamanAnalytics
   RamanSimulator/RamanSimulator



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
