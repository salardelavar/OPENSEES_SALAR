@Echo Off
Cls
Echo. #####################################################################
Echo. #                       IN THE NAME OF ALLAH                        #
Echo. # ELASTIC AND INELASTIC ANLYSIS OF SDOF STRUCTURE RESPONSE SPECTRUM #
Echo. #                       OUTPUT OF ANALYSIS:                         #
Echo. #                     PEACK GROUND DISPLACEMENT (mm)                 #
Echo. #                     PEACK GROUND VELOCITY     (mm/s)               #
Echo. #                     PEACK GROUND ACCELARATION (mm/s2)              #
Echo. #####################################################################
Echo.
Echo.
Echo.
Echo.
Echo.
Echo.           ## ELASTIC RESPONSE SPECTRUM ANALYSIS ##
OpenSees.exe E_SPEC.tcl
ren Data E_SPEC
Echo.           ## INELASTIC RESPONSE SPECTRUM ANALYSIS - UNCONFINED CONCRETE SECTION ##
OpenSees.exe INE_SPEC_UNCONFINED.tcl
ren Data INE_SPEC_UNCONFINED
Echo.           ## INELASTIC RESPONSE SPECTRUM ANALYSIS - CONFINED CONCRETE SECTION ##
OpenSees.exe INE_SPEC_CONFINED.tcl
ren Data INE_SPEC_CONFINED
PAUSE
