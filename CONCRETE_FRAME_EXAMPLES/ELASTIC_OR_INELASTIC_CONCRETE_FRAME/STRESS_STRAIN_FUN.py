import openseespy.opensees as ops
# OUTPUT STRESS-STRAIN OF CONRETE SECTION
#%% -----------------------------------------------------------------------
def STRESS_STRAIN_CONCRETE(NAME_TAG, ELE_TAG, SECTION_NUM, COOR_FIBER_Y, COOR_FIBER_Z): # (-COOR_FIBER_Y, 0) ** (+COOR_FIBER_Y, 0)
    #ops.recorder('Element', '-file', f'{NAME_TAG}_CONCRETE.txt', '-time', '-ele', ELE_TAG, 'section', SECTION_NUM, 'fiber', -COOR_FIBER_Y, 0, 'stressStrainTangent')
    ops.recorder('Element', '-file', f'{NAME_TAG}_CONCRETE.txt', '-time', '-ele', ELE_TAG, 'section', SECTION_NUM, 'fiber', -COOR_FIBER_Y, 0, 'stressStrain')

#%% -----------------------------------------------------------------------
def STRESS_STRAIN_REBAR(NAME_TAG, ELE_TAG, SECTION_NUM, COOR_FIBER_Y, COOR_FIBER_Z):# (-DD+COVER, -0.5*B+COVER) ** (+DD-COVER, -0.5*B+COVER) 
    ops.recorder('Element', '-file', f'{NAME_TAG}_REBAR.txt', '-time', '-ele', ELE_TAG, 'section', SECTION_NUM, 'fiber', COOR_FIBER_Y, COOR_FIBER_Z, 'stressStrain')

#%% -----------------------------------------------------------------------    
def STRESS_STRAIN_PLOT(ELE_TAG):
    def OUTPUT_SECOND_COLUMN(X, COLUMN):
        import numpy as np
        # Time History
        filename = f"{X}.txt"
        data_collected = np.loadtxt(filename)
        X = data_collected[:, COLUMN]   
        return X 
    
    #%% STRESS-STRAIN OF ELEMENT ELE_TAG
    # CONCRETE
    strain_B_C = OUTPUT_SECOND_COLUMN(f'0{ELE_TAG}_STRESS_STRAIN_BOT_CONCRETE', 2) # Reading bottom concrete strain from Text file
    stress_B_C = OUTPUT_SECOND_COLUMN(f'0{ELE_TAG}_STRESS_STRAIN_BOT_CONCRETE', 1) # Reading bottom concrete stress from Text file
    strain_T_C = OUTPUT_SECOND_COLUMN(f'0{ELE_TAG}_STRESS_STRAIN_TOP_CONCRETE', 2) # Reading top concrete strain from Text file
    stress_T_C = OUTPUT_SECOND_COLUMN(f'0{ELE_TAG}_STRESS_STRAIN_TOP_CONCRETE', 1) # Reading top concrete stress from Text file
    # STEEL REBAR
    strain_B_R = OUTPUT_SECOND_COLUMN(f'0{ELE_TAG}_STRESS_STRAIN_BOT_REBAR', 2) # Reading bottom steel rebar strain from Text file
    stress_B_R = OUTPUT_SECOND_COLUMN(f'0{ELE_TAG}_STRESS_STRAIN_BOT_REBAR', 1) # Reading bottom steel rebar stress from Text file
    strain_T_R = OUTPUT_SECOND_COLUMN(f'0{ELE_TAG}_STRESS_STRAIN_TOP_REBAR', 2) # Reading top steel rebar strain from Text file
    stress_T_R = OUTPUT_SECOND_COLUMN(f'0{ELE_TAG}_STRESS_STRAIN_TOP_REBAR', 1) # Reading top steel rebar stress from Text file
    
    # PLOT
    import matplotlib.pyplot as plt
    plt.figure(100+ELE_TAG, figsize=(8, 6))
    plt.plot(strain_B_C, stress_B_C, color='blue', label='Bottom Fiber', linewidth=2)
    plt.plot(strain_T_C, stress_T_C, color='red', label='Top Fiber', linewidth=2)
    plt.xlabel('Strain (mm/mm)')
    plt.ylabel('Stress (N/mm^2)')
    plt.title(f'Stress-Strain Relation of Element {ELE_TAG} Concrete Top & Bottom Fibers')
    plt.grid()
    plt.legend()
    plt.show()
    
    plt.figure(200+ELE_TAG, figsize=(8, 6))
    plt.plot(strain_B_R, stress_B_R, color='blue', label='Bottom Fiber', linewidth=2)
    plt.plot(strain_T_R, stress_T_R, color='red', label='Top Fiber', linewidth=2)
    plt.xlabel('Strain (mm/mm)')
    plt.ylabel('Stress (N/mm^2)')
    plt.title(f'Stress-Strain Relation of Element {ELE_TAG} Steel Rebar Top & Bottom Fibers')
    plt.grid()
    plt.legend()
    plt.show()
#%% -----------------------------------------------------------------------    
    