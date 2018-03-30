# gmosredux

These scripts takes the raw data from the GOA and and reduced them to the point of 2D transformed spectra for later use in our PYPIT pipeline. Unfortunately still uses PYRAF and the Gemini IRAF software.

Requirements:
  - Anaconda installation of the astroconda gemini package (http://www.gemini.edu/node/11823)


NOTE:
When reducing the standards, highly recommend using the notebook in standards. There are many things to change for each standard and some interactivity. Also, the notebook better follows the gemini cookbook (http://ast.noao.edu/sites/default/files/GMOS_Cookbook/Processing/PyrafProcMOS.html) rather than Andy Caseys (very helpful) code. 

NOTE: 
Might as well just use the notebook for running the MOS redux too. Wont have to go back and rerun the whole thing when things go wrong. Horrible I know. IRAF sucks. 
