Installation
============

.. NOTE::
  This guide only explains how to add the extended functionality to CLM and not how CESM can be installed.


Obtaining the Source Code
-------------------------
The source code must first be obtained from github:
``git clone https://github.com/IACETH/prescribeSM_cesm_1.2.x.git``


Compilation
-----------
There is one additional source code file and a number of original CLM files were changed (see :ref:`source_files`). 
The files need to go in to the :file:`SourceMods` folder 
(:file:`$CASEROOT/SourceMods/src.clm/`) before CESM is compiled.

