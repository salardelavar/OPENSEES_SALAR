def GRAVITY_ANALYSIS_FUN(NstepGravity, MAX_TOLERANCE, MAX_ITERATIONS):
    # WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    import openseespy.opensees as ops
    import ANALYSIS_FUNCTION as S02
    # Gravity Analysis
    #NstepGravity = 10
    DGravity = 1/NstepGravity
    ops.integrator('LoadControl', DGravity) # determine the next time step for an analysis
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/loadControl.html
    ops.numberer('Plain') # renumber dof's to minimize band-width (optimization), if you want to
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/numberer.html
    ops.system('BandGeneral') # how to store and solve the system of equations in the analysis
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/system.html
    ops.constraints('Plain') # how it handles boundary conditions
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/constraints.html
    ops.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 0) # determine if convergence has been achieved at the end of an iteration step
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/test.html
    ops.algorithm('Newton') # use Newton's solution algorithm: updates tangent stiffness at every iteration
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/algorithm.html
    ops.analysis('Static') # define type of analysis static or transient
    # INFO LINK: https://openseespydoc.readthedocs.io/en/stable/src/analysis.html
    OK = ops.analyze(NstepGravity) # apply gravity
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/analyze.html
    #S02.ANALYSIS(OK, NstepGravity1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
    print('Graviy Analysis Done.')