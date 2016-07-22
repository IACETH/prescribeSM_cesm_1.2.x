Usage
=====

#. Run a reference simulation with CESM/ CLM, :doc:`output</SMforcing>` ``SOILLIQ`` and ``SOILICE`` daily data.
#. Process the output, e.g. calculate a climatology.
#. :doc:`Convert</SMforcing>` the :term:`4D` output file to :term:`3D`, if necessary.
#. :doc:`Add the necessary files and compile CESM.</installation>`
#. Create the :doc:`namelist</namelist>` file, deciding what :doc:`method</methods>` to use to prescribe SM.
#. Run CESM with prescribed SM