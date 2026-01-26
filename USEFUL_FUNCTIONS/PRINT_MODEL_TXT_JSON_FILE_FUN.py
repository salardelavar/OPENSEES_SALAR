def PRINT_MODEL_TXT_JSON_FILE_FUN():
    import openseespy.opensees as ops
    # Compute and print modal properties
    ops.modalProperties("-print", "-file", "SALAR_ModalReport.txt", "-unorm")
    
    # Print log all messages and errors in a file
    ops.logFile("SALAR_logReport.txt", "-append", "-noEcho")
    
    # Print the Model
    ops.printModel()
    ops.printModel("-JSON", "-file", "SALAR_Print_Model.json")
    
    ops.printA("-file", "SALAR_Print_A.txt", "-ret") # Print the contents of a FullGeneral system that the integrator creates to the screen or a file
    ops.printB("-file", "SALAR_Print_B.txt", "-ret") # Print the right hand side of a FullGeneral system that the integrator creates to the screen or a file