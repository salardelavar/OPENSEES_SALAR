#          ######################################################################################################
#          #                                           IN THE NAME OF ALLAH                                     #
#          #                   STEEL-CONCRETE COMPOSITE PLATE GIRDERS BRIDGE SUPERSTRUCTURE                     #
#          #               RUNNING MOMENT-CURVATURE, PUSHOVER AND DYNAMIC ANALYSIS FOR CALCULATE                #
#          #                               OPTIMUM STRUCTURAL DUCTILIY DAMAGE INEX                              #
#          #----------------------------------------------------------------------------------------------------#
#          #                    THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)                      #
#          #                               EMAIL: salar.d.ghashghaei@gmail.com                                  #
#          ######################################################################################################
"""
Key Points:
[1]  Performs advanced nonlinear analyses on steel-concrete composite plate girder bridges.  
[2]  Implements fiber-section modeling for accurate material behavior representation.  
[3]  Conducts moment-curvature analysis to determine section-level ductility parameters.  
[4]  Executes displacement-controlled pushover analysis for capacity curve development.  
[5]  Performs nonlinear time history analysis with real ground motion records.  
[6]  Utilizes Newmark-beta integration for dynamic solution stability.  
[7]  Implements Rayleigh damping with mass/stiffness proportionality.  
[8]  Calculates structural period and mode shapes through eigenvalue analysis.  
[9]  Determines bilinear approximations for force-deformation relationships.  
[10] Computes ductility ratios and overstrength factors.  
[11] Quantifies structural damage via ductility-based damage indices.  
[12] Employs logarithmic decrement method for damping estimation.  
[13] Incorporates advanced concrete material models (Concrete02).  
[14] Uses hysteretic steel material models for cyclic behavior.  
[15] Implements displacement control with adaptive algorithms.  
[16] Performs large-deformation analysis capturing geometric nonlinearities.  
[17] Evaluates P-M interaction effects on structural capacity.  
[18] Implements rigorous convergence criteria (NormDispIncr).  
[19] Analyzes stiffness degradation and strength deterioration.  
[20] Integrates multiple analysis types for comprehensive seismic performance assessment.  

This implementation demonstrates state-of-the-art techniques in performance-based earthquake engineering,
 combining material-level modeling with system-level response evaluation to quantify structural damage
 under seismic loads. The methodology follows modern bridge design philosophy where ductility is explicitly
 designed and evaluated as a key performance metric.
"""
#import the os module
import os
import time as TI
import numpy as np
import openseespy.opensees as op
import opsvis as opsv
import matplotlib.pyplot as plt
import ANALYSIS_FUNCTION as S02
import PLOT_2D as S04
import SALAR_MATH as S05
#%%-------------------------------------------------
# Create a directory at specified path with name 'directory_path'
#import os
directory_path = 'C:\\OPENSEESPY_SALAR'

# Check if the directory already exists
if not os.path.exists(directory_path):
    os.mkdir(directory_path)
    print(f"Directory '{directory_path}' created successfully.")
else:
    print(f"Directory '{directory_path}' already exists. Skipping creation.")
#-------------------------------------------------------------
# Create folder name
FOLDER_NAME = 'COMPOSITE_BRIDGE_OPTIMIZATION'
dir = f"C:\\OPENSEESPY_SALAR\\{FOLDER_NAME}\\"
if not os.path.exists(dir):
    os.makedirs(dir)    
#-------------------------------------------------------------
# OUTPUT DATA ADDRESS:
SALAR_DIR = f'C://OPENSEESPY_SALAR//{FOLDER_NAME}//';
#-------------------------------------------------------------
## DELETE ALL FILES IN DIRECTORY 
def DELETE_FOLDER_CONTANTS(folder_path):
    import os
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")
    print("Deletion done")
   
FOLDER_PATH = f'C:\\OPENSEESPY_SALAR\\{FOLDER_NAME}'  # Specify the folder path
#DELETE_FOLDER_CONTANTS(FOLDER_PATH)  
#%%-------------------------------------------------
# Load the image
def PLOT_IMAGE(image):
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    image = mpimg.imread(image_path)

    # Display the image
    plt.figure(figsize=(30, 16))
    plt.imshow(image)
    plt.axis('off')  # Hide axes
    plt.show()
    
image_path = 'SINGLE_SPAN_SIMPLY_SUPPORTED_COMPOSITE_BRIDGE_00.png'    
PLOT_IMAGE(image_path)
#%%-------------------------------------------------
# BRIDGE IMAGES AND DATA
'https://www.steelconstruction.info/A1_Peterborough_to_Blyth'
# BOOK: Finite Element Analysis and Design of Steel and Steel–Concrete Composite Bridges
'https://www.sciencedirect.com/book/9780124172470/finite-element-analysis-and-design-of-steel-and-steel-concrete-composite-bridges'
# BOOK: Steel-concrete composite bridge design guide
'https://www.nzta.govt.nz/assets/resources/research/reports/525/docs/525.pdf' 
# PAPER: Seismic damage prediction by deterministic methods: Concepts and procedures
'https://onlinelibrary.wiley.com/doi/10.1002/eqe.4290160507'
#%%-------------------------------------------------
### -------------------------
###  DEFINE AND PLOT SECTION
### -------------------------
def Composite_Bridge_Section_Plot(Section_Tag, nFibZ, nFibY, SideWalk_z, SideWalk_y, SideWalk_C, Deck_z,
                                 Deck_y, Box_z, Box_y, tf01, bf, tf02, tw, DA01, DA02, DA03, PLOT, OUTPUT, I):

    import openseespy.opensees as op
    import opsvis as opsv
    import numpy as np
    
    Mat_Tag01 = 1 # Confined Concrete Section Tag
    Mat_Tag02 = 2 # Steel Plate Section Tag
    Mat_Tag03 = 3 # Steel Rebar Section Tag

    
    fc = -35                 # [N/mm^2] Nominal concrete compressive strength
    Ec = 4700 * np.sqrt(-fc) # [N/mm^2] Concrete Elastic Modulus

    # confined concrete
    Kfc = 1.3;			# ratio of confined to unconfined concrete strength
    fc1C = Kfc*fc;		# CONFINED concrete (mander model), maximum stress
    KfcB = 1.1;			# ratio of confined to unconfined concrete strength
    fc1CB = KfcB*fc;	# CONFINED concrete (mander model), maximum stress

    eps1C = 2*fc1C/Ec;	    # strain at maximum stress 
    fc2C = 0.2*fc1C;		# ultimate stress
    eps2C = 5*eps1C;		# strain at ultimate stress 
    # unconfined concrete
    fc1U = fc;			    # UNCONFINED concrete (todeschini parabolic model), maximum stress
    eps1U = -0.0025;	    # strain at maximum strength of unconfined concrete
    fc2U = 0.2*fc1U;		# ultimate stress
    eps2U = -0.012;			# strain at ultimate stress
    Lambda = 0.1;				# ratio between unloading slope at $eps2 and initial slope $Ec
    # tensile-strength properties
    ftC = -0.55*fc1C;		# tensile strength +tension
    ftU = -0.55*fc1U;		# tensile strength +tension
    Ets = ftU/0.002;		# tension softening stiffness
    
    op.uniaxialMaterial('Concrete02', Mat_Tag01, fc1C, eps1C, fc2C, eps2C, Lambda, ftC, Ets) # build core concrete (confined)
    
    # PLATE MATERIAL PROPERTIES:
    """
    FyP = 360			# Steel rebar yield stress
    CyP = 0.0018	    # Steel rebar yield strain
    EsP = FyP/CyP		# modulus of steel
    BsP = 0.01			# strain-hardening ratio 
    R0P = 18.0			# control the transition from elastic to plastic branches
    cR1P = 0.925		# control the transition from elastic to plastic branches
    cR2P = 0.15			# control the transition from elastic to plastic branches
    """    
    fy = 360          # [N/mm^2] Yield strength
    Es = 210e4        # [N/mm^2] Young's modulus 
    ey = fy/Es        # [mm/mm] Steel Yield Strain
    fu = 1.1818*fy    # [N/mm²] Steel Ultimate Strength
    esu = ey*75.2     # [mm/mm] Steel Ultimate Strain
    pinchX = 0.8   # Pinching factor in X direction
    pinchY = 0.5   # Pinching factor in Y direction
    damage1 = 0.0  # Damage due to ductility
    damage2 = 0.0  # Damage due to energy
    beta = 0.1     # Stiffness degradation parameter
    op.uniaxialMaterial('Hysteretic', Mat_Tag02, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -0.5*esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
    # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
    # REBAR MATERIAL PROPERTIES:
    """
    Fy = 4000			    # Steel rebar yield stress
    Cy = 0.02			    # Steel rebar yield strain
    Es = Fy/Cy				# modulus of steel
    Bs = 0.01				# strain-hardening ratio 
    R0 = 18.0				# control the transition from elastic to plastic branches
    cR1 = 0.925				# control the transition from elastic to plastic branches
    cR2 = 0.15				# control the transition from elastic to plastic branches
    """

    fy = 4000         # [N/mm^2] Yield strength
    Es = 210e4        # [N/mm^2] Young's modulus 
    ey = fy/Es        # [mm/mm] Steel Yield Strain
    fu = 1.1818*fy    # [N/mm²] Steel Ultimate Strength
    esu = ey*75.2     # [mm/mm] Steel Ultimate Strain
    pinchX = 0.8   # Pinching factor in X direction
    pinchY = 0.5   # Pinching factor in Y direction
    damage1 = 0.0  # Damage due to ductility
    damage2 = 0.0  # Damage due to energy
    beta = 0.1     # Stiffness degradation parameter
    op.uniaxialMaterial('Hysteretic', Mat_Tag03, fy, ey, fu, esu, 0.2*fu, 1.1*esu, -fy, -ey, -fu, -0.5*esu, -0.2*fu, -1.1*esu, pinchX, pinchY, damage1, damage2, beta)
    # INFO LINK: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material
    
    #op.uniaxialMaterial('Steel02', Mat_Tag02, FyP, EsP, BsP, R0P,cR1P,cR2P) # build plate material
    #op.uniaxialMaterial('Steel02', Mat_Tag03, Fy, Es, Bs, R0,cR1,cR2) # build reinforcement material
    


    #nFibZ = 5 # NUMBER OF FIBER SECTIONIN Z DIRECTION
    #nFibY = 10 # NUMBER OF FIBER SECTIONIN Y DIRECTION
    As01 = (np.pi * (0.5 * DA01)**2) # [mm^2] SideWalk Rebar Area
    As02 = (np.pi * (0.5 * DA02)**2) # [mm^2] Deck Rebar Area
    As03 = (np.pi * (0.5 * DA03)**2) # [mm^2] Box Rebar Area
    
    hw = H - tf01 - tf02# Top Web Plate Height Beam
    # SIDEWALK
    rey1 = SideWalk_y - cover
    rez1 = SideWalk_z - cover
    # DECK
    rey2 = Deck_y - cover
    rez2 = Deck_z - cover
    # BOX
    rey3 = Box_y - cover
    rez3 = Box_z - cover

    FIBER_SEC = [['section', 'Fiber', Section_Tag, '-GJ', 1.0e6],
                 # CONCRETE SECION SIDEWALK
                 ['patch', 'rect', Mat_Tag01, nFibY, nFibZ, SideWalk_y, -Deck_z+SideWalk_C, 0, -Deck_z-SideWalk_z],# Left Side Walk
                 ['patch', 'rect', Mat_Tag01, nFibY, nFibZ, SideWalk_y, Deck_z+SideWalk_z, 0, +Deck_z-SideWalk_C],# Right Side Walk
                 # CONCRETE SECION DECK
                 ['patch', 'rect', Mat_Tag01, nFibY, nFibZ, -Deck_y, -Deck_z, 0, +Deck_z],# Road Deck
                 ['patch', 'rect', Mat_Tag01, nFibY, nFibZ, -Deck_y, -Deck_z, 0, +Deck_z],# Road Deck
                 # CONCRETE SECION BOX
                 ['patch', 'rect', Mat_Tag01, nFibY, nFibZ, -Deck_y-Box_y, -Deck_z+Z1, -Deck_y, -Deck_z+Box_z+Z1],# Road Deck Box Right 01
                 ['patch', 'rect', Mat_Tag01, nFibY, nFibZ, -Deck_y-Box_y, -Deck_z+Box_z+Z1+Z2, -Deck_y, -Deck_z+Box_z+(2*Z1+Z2)],# Road Deck Box Right 02
                 ['patch', 'rect', Mat_Tag01, nFibY, nFibZ, -Deck_y-Box_y, -Box_z*0.5, -Deck_y, Box_z*0.5],# Road Deck Box Center
                 ['patch', 'rect', Mat_Tag01, nFibY, nFibZ, -Deck_y-Box_y, +Deck_z-Box_z-Z1, -Deck_y, +Deck_z-Z1],# Road Deck Box Left 01
                 ['patch', 'rect', Mat_Tag01, nFibY, nFibZ, -Deck_y-Box_y, +Deck_z-Box_z-(2*Z1+Z2), -Deck_y, +Deck_z-Box_z-(Z1+Z2)],# Road Deck Box Left 02
                 # STEEL TOP FLANGE PLATE
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01, -Deck_z+Z1, -Deck_y-Box_y, -Deck_z+Box_z+Z1],# Road Deck Top Flange Right 01
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01, -Deck_z+Box_z+(Z1+Z2), -Deck_y-Box_y, -Deck_z+Box_z+(2*Z1+Z2)],# Road Deck Top Flange Right 02
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01, -Box_z*0.5, -Deck_y-Box_y, Box_z*0.5],# Road Deck Top Flange Center
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01, +Deck_z-Box_z-Z1, -Deck_y-Box_y, +Deck_z-Z1],# Road Deck Top Flange Left 01
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01, +Deck_z-Box_z-(2*Z1+Z2), -Deck_y-Box_y, +Deck_z-Box_z-(Z1+Z2)],# Road Deck Top Flange Left 02
                 # STEEL BOTTOM FLANGE PLATE
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01-hw-tf02, -Deck_z+Z1, -Deck_y-Box_y-tf01-hw, -Deck_z+Box_z+Z1],# Road Deck Top Flange Right 01
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01-hw-tf02, -Deck_z+Box_z+(Z1+Z2), -Deck_y-Box_y-tf01-hw, -Deck_z+Box_z+(2*Z1+Z2)],# Road Deck Top Flange Right 02
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01-hw-tf02, -Box_z*0.5, -Deck_y-Box_y-tf01-hw, Box_z*0.5],# Road Deck Top Flange Center
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01-hw-tf02, +Deck_z-Box_z-Z1, -Deck_y-Box_y-tf01-hw, +Deck_z-Z1],# Road Deck Top Flange Left 01
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01-hw-tf02, +Deck_z-Box_z-(2*Z1+Z2), -Deck_y-Box_y-tf01-hw, +Deck_z-Box_z-(Z1+Z2)],# Road Deck Top Flange Left 02
                 # STEEL WEB PLATE
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01-hw, -Deck_z+0.5*Box_z+Z1-0.5*tw, -Deck_y-Box_y-tf01, -Deck_z+0.5*Box_z+Z1+0.5*tw],# Road Deck Web Right 01
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01-hw, -Deck_z+0.5*Box_z+(2*Z1+Z2)-0.5*tw, -Deck_y-Box_y-tf01, -Deck_z+0.5*Box_z+(2*Z1+Z2)+0.5*tw],# Road Deck Web Right 02
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01-hw, -Deck_z+0.5*Box_z+(3*Z1+2*Z2)-0.5*tw, -Deck_y-Box_y-tf01, -Deck_z+0.5*Box_z+(3*Z1+2*Z2)+0.5*tw],# Road Deck Web Right 03
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01-hw, -Deck_z+0.5*Box_z+(4*Z2)-0.5*tw, -Deck_y-Box_y-tf01, -Deck_z+0.5*Box_z+(4*Z2)+0.5*tw],# Road Deck Web Right 04
                 ['patch', 'rect', Mat_Tag02, nFibY, nFibZ, -Deck_y-Box_y-tf01-hw, -Deck_z+0.5*Box_z+(Z1+5*Z2)-0.5*tw, -Deck_y-Box_y-tf01, -Deck_z+0.5*Box_z+(Z1+5*Z2)+0.5*tw],# Road Deck Web Right 05
                 # STEEL REBAR
                 ['layer', 'straight', Mat_Tag03, 10, As01, rey1, Deck_z+rez1, rey1, Deck_z-SideWalk_C+cover],# TOP REBAR LAYERS - SIDEWALK RIGHT
                 ['layer', 'straight', Mat_Tag03, 10, As01, cover, Deck_z+rez1, cover, Deck_z-SideWalk_C+cover],# BOTTOM REBAR LAYERS - SIDEWALK RIGHT
                 ['layer', 'straight', Mat_Tag03, 10, As01, rey1, -(Deck_z+rez1), rey1, -(Deck_z-SideWalk_C+cover)],# TOP REBAR LAYERS - SIDEWALK LEFT
                 ['layer', 'straight', Mat_Tag03, 10, As01,  cover, -(Deck_z+rez1),  cover, -(Deck_z-SideWalk_C+cover)],# BOTTOM REBAR LAYERS - SIDEWALK LEFT
                 ['layer', 'straight', Mat_Tag03, 65, As02, -rey2, -rez2, -rey2, rez2],# BOTTOM REBAR LAYERS - DECK
                 ['layer', 'straight', Mat_Tag03, 65, As02, -Deck_y*0.5, -rez2, -Deck_y*0.5, rez2],# MIDDLE REBAR LAYERS - DECK
                 ['layer', 'straight', Mat_Tag03, 65, As02, -cover, -rez2, -cover, rez2],# TOP REBAR LAYERS - DECK
                 ['layer', 'straight', Mat_Tag03,  3, As03, -(Deck_y+rey3), -Deck_z+Z1+cover, -(Deck_y+rey3), -Deck_z+Box_z+Z1-cover],# BOTTOM REBAR LAYERS - BOX 01
                 ['layer', 'straight', Mat_Tag03,  3, As03, -(Deck_y+cover), -Deck_z+Z1+cover, -(Deck_y+cover), -Deck_z+Box_z+Z1-cover],# TOP REBAR LAYERS - BOX 01
                 ['layer', 'straight', Mat_Tag03,  3, As03, -(Deck_y+rey3), -Deck_z+Box_z+Z1+Z2+cover, -(Deck_y+rey3), -Deck_z+Box_z+(2*Z1+Z2)-cover],# BOTTOM REBAR LAYERS - BOX 02
                 ['layer', 'straight', Mat_Tag03,  3, As03, -(Deck_y+cover), -Deck_z+Box_z+Z1+Z2+cover, -(Deck_y+cover), -Deck_z+Box_z+(2*Z1+Z2)-cover],# TOP REBAR LAYERS - BOX 02
                 ['layer', 'straight', Mat_Tag03,  3, As03, -(Deck_y+rey3), -Box_z*0.5+cover, -(Deck_y+rey3), Box_z*0.5-cover],# BOTTOM REBAR LAYERS - BOX 03 - MIDDLE
                 ['layer', 'straight', Mat_Tag03,  3, As03, -(Deck_y+cover), -Box_z*0.5+cover, -(Deck_y+cover), Box_z*0.5-cover],# TOP REBAR LAYERS - BOX 03 - MIDDLE
                 ['layer', 'straight', Mat_Tag03,  3, As03, -(Deck_y+rey3), Deck_z-Z1-cover, -(Deck_y+rey3), Deck_z-Box_z-Z1+cover],# BOTTOM REBAR LAYERS - BOX 05
                 ['layer', 'straight', Mat_Tag03,  3, As03, -(Deck_y+cover), Deck_z-Z1-cover, -(Deck_y+cover), Deck_z-Box_z-Z1+cover],# TOP REBAR LAYERS - BOX 05
                 ['layer', 'straight', Mat_Tag03,  3, As03, -(Deck_y+rey3), Deck_z-Box_z-Z1-Z2-cover, -(Deck_y+rey3), Deck_z-Box_z-(2*Z1+Z2)+cover],# BOTTOM REBAR LAYERS - BOX 04
                 ['layer', 'straight', Mat_Tag03,  3, As03, -(Deck_y+cover), Deck_z-Box_z-Z1-Z2-cover, -(Deck_y+cover), Deck_z-Box_z-(2*Z1+Z2)+cover],# TOP REBAR LAYERS - BOX 04

                ]

    if PLOT == 1:
        matcolor = ['gold', 'lightgrey']
        plt.figure(1)
        opsv.plot_fiber_section(FIBER_SEC, matcolor=matcolor)
        # Set the x and y limits
        plt.ylim(-1000, 1000)
        plt.xlim(-1500,-6850) # RIGHT
        plt.show()
        
        plt.figure(2)
        opsv.plot_fiber_section(FIBER_SEC, matcolor=matcolor)
        plt.ylim(-2500, 300)
        plt.xlim(-1500,-6850) #  RIGHT
        plt.show()

        plt.figure(3)
        opsv.plot_fiber_section(FIBER_SEC, matcolor=matcolor)
        plt.ylim(-1000, 1000)
        plt.xlim(600,-600) #  MIDDLE
        plt.show()

        plt.figure(4)
        opsv.plot_fiber_section(FIBER_SEC, matcolor=matcolor)
        plt.ylim(-2500, 300)
        plt.xlim(6850,1500) # LEFT
        plt.show()
        
        plt.figure(5)
        opsv.plot_fiber_section(FIBER_SEC, matcolor=matcolor)
        plt.ylim(-1000, 1000)
        plt.xlim(6850,1500) # LEFT
        plt.show()

        plt.figure(6)
        opsv.plot_fiber_section(FIBER_SEC, matcolor=matcolor)
        plt.axis('equal')
        plt.show()
        
    # CONCRETE FIBER
    op.recorder('Element','-ele', 50,'-file',f'{SALAR_DIR}{OUTPUT}_fiberCon_StressStrain_01_{I}.txt','section', 5,'fiber', SideWalk_y, +Deck_z-SideWalk_C+SideWalk_z,'stressStrain')# concrete fiber 01
    op.recorder('Element','-ele', 50,'-file',f'{SALAR_DIR}{OUTPUT}_fiberCon_StressStrain_02_{I}.txt','section', 5,'fiber', 0, +Deck_z,'stressStrain')# concrete fiber 02
    op.recorder('Element','-ele', 50,'-file',f'{SALAR_DIR}{OUTPUT}_fiberCon_StressStrain_03_{I}.txt','section', 5,'fiber', -Deck_y, +Deck_z,'stressStrain')# concrete fiber 03
    # STEEL REBAR FIBER
    op.recorder('Element','-ele', 50,'-file',f'{SALAR_DIR}{OUTPUT}_fiberReb_StressStrain_01_{I}.txt','section', 5,'fiber', rey1, Deck_z+rez1,'stressStrain')# steel rebar fiber 01 - TOP REBAR LAYERS - SIDEWALK RIGHT
    op.recorder('Element','-ele', 50,'-file',f'{SALAR_DIR}{OUTPUT}_fiberReb_StressStrain_02_{I}.txt','section', 5,'fiber', cover, Deck_z+rez1,'stressStrain')# steel rebar fiber 02 - BOTTOM REBAR LAYERS - SIDEWALK RIGHT
    op.recorder('Element','-ele', 50,'-file',f'{SALAR_DIR}{OUTPUT}_fiberReb_StressStrain_03_{I}.txt','section', 5,'fiber', -cover, -rez2,'stressStrain')# steel rebar fiber 03 - TOP REBAR LAYERS - DECK
    op.recorder('Element','-ele', 50,'-file',f'{SALAR_DIR}{OUTPUT}_fiberReb_StressStrain_04_{I}.txt','section', 5,'fiber', -Deck_y*0.5, -rez2,'stressStrain')# steel rebar fiber 04 - MIDDLE REBAR LAYERS - DECK
    op.recorder('Element','-ele', 50,'-file',f'{SALAR_DIR}{OUTPUT}_fiberReb_StressStrain_05_{I}.txt','section', 5,'fiber', -rey2, -rez2,'stressStrain')# steel rebar fiber 05 - BOTTOM REBAR LAYERS - DECK
    # STEEL PLATE FIBER
    op.recorder('Element','-ele', 50,'-file',f'{SALAR_DIR}{OUTPUT}_fiberPlate_StressStrain_01_{I}.txt','section', 5,'fiber', -Deck_y-Box_y, +Deck_z-Box_z-(2*Z1+Z2),'stressStrain')# steel plate fiber 01 - TOP PLATE LAYERS
    op.recorder('Element','-ele', 50,'-file',f'{SALAR_DIR}{OUTPUT}_fiberPlate_StressStrain_02_{I}.txt','section', 5,'fiber', -Deck_y-Box_y-tf01-hw-tf02, +Deck_z-Box_z-(2*Z1+Z2),'stressStrain')# steel plate fiber 02 - BOTTOM PLATE LAYERS        
    
    return FIBER_SEC
#%%-------------------------------------------------
### -------------------------------
###    MOMENT-CURVATURE FUNCTION
### -------------------------------

def MC_ANALYSIS(Section_Tag, PX, PY, MZ, DR, numIncr, MAX_ITERATIONS, MAX_TOLERANCE, J):
    ####      Start Moment Curvature Analysis
    # Define model builder
    # --------------------
    op.wipe() # Reset model
    op.model('basic','-ndm',2,'-ndf',3)
    
    nFibZ = 5  # NUMBER OF FIBER SECTIONIN Z DIRECTION
    nFibY = 25 # NUMBER OF FIBER SECTIONIN Y DIRECTION
    FS = Composite_Bridge_Section_Plot(Section_Tag, nFibZ, nFibY, SideWalk_z, SideWalk_y, SideWalk_C,
                                      Deck_z, Deck_y, Box_z, Box_y, tf01, bf, tf02, tw, DA01, DA02, DA03,
                                      PLOT=0, OUTPUT='MC', I=J)
    
    opsv.fib_sec_list_to_cmds(FS)
    
    H = SideWalk_y + Deck_y + Box_y + tf01 + tw + tf02 # I Section Height
    # -----------------------------------------------------------------------------------------------------
    # Yield Concrete Strain
    Cy = 0.004
    # Set Yield Curvature
    Ky = Cy / (0.5 * H)
    #print('Ky', Ky)

    # set ultimate Curvature
    Ku = Ky * DR
    #print('Ku', Ku)


    # Define two nodes at (0,0)
    op.node(1, 0.0, 0.0)
    op.node(2, 0.0, 0.0)

    # Fix all degrees of freedom except axial and bending
    op.fix(1, 1, 1, 1)
    op.fix(2, 0, 1, 0)

    # Define element
    #                             tag ndI ndJ  secTag
    op.element('zeroLengthSection',  1,   1,   2,  Section_Tag)
    # Create recorder
    op.recorder('Node', '-file', f"{SALAR_DIR}CUR_{J}.txt",'-time', '-node', 2, '-dof', 3, 'disp')# Curvature Time History nodes 2
    op.recorder('Node', '-file', f"{SALAR_DIR}MOM_{J}.txt",'-time', '-node', 1, '-dof', 3, 'reaction')# Base Shear Time History nodes 1

    # Define constant axial load
    op.timeSeries('Constant', 1)
    op.pattern('Plain', 1, 1)
    op.load(2, PX, PY, MZ)

    # Define analysis parameters
    op.integrator('LoadControl', 0.001)
    op.system('SparseGeneral', '-piv')
    op.test('NormUnbalance', MAX_TOLERANCE, MAX_ITERATIONS)
    op.numberer('Plain')
    op.constraints('Plain')
    op.algorithm('Newton')
    op.analysis('Static')

    # Do one analysis for constant axial load
    op.analyze(1)

    # Define reference moment
    op.timeSeries('Linear', 2)
    op.pattern('Plain',2, 2)
    op.load(2, 0.0, 0.0, 1.0)

    # Compute curvature increment
    dK = Ku / numIncr

    # Use displacement control at node 2 for section analysis
    op.integrator('DisplacementControl', 2, 3, dK, 1, dK, dK)

    # Do the section analysis
    #op.analyze(numIncr)
        
    # Data storage
    FORCE_S, FORCE_A, MOMENT = [], [], []
    DISP_X, DISP_Y, ROT = [], [], []
    KA, KS, KI, STEP = [], [], [], []
        
    for step in range(numIncr):
        OK = op.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        # Record results
        op.reactions()
        S = op.nodeReaction(1, 1)  # SHEAR BASE REACTION
        A = op.nodeReaction(1, 2)  # AXIAL BASE REACTION
        M = op.nodeReaction(1, 3)  # MOMENT BASE REACTION
        #print(rot, M)
        disp_X = op.nodeDisp(2, 1) # LATERAL DISPLACEMENT IN X FOR NODE IN MIDDLE ON BRIDGE LENGTH
        disp_Y = op.nodeDisp(2, 2) # LATERAL DISPLACEMENT IN Y FOR NODE IN MIDDLE ON BRIDGE LENGTH
        rot = op.nodeDisp(2, 3)    # ROTATION IN Z FOR NODE IN MIDDLE ON BRIDGE LENGTH
        FORCE_S.append(S)
        FORCE_A.append(A)
        MOMENT.append(M)
        DISP_X.append(disp_X)
        DISP_Y.append(disp_Y)
        ROT.append(rot)
        KS.append(np.abs(S)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
        KA.append(np.abs(A)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
        KI.append(np.abs(M)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
        STEP.append(step)
        print(step+1, rot, M)
        
    print(f'{J+1} MC Done.')
    #op.wipe()

#%%-------------------------------------------------
### -----------------------
###    PUSHOVER FUNCTION
### -----------------------

def PUSHOVER_ANALYSIS(Length, SideWalk_z, SideWalk_y, SideWalk_C, Deck_z, Deck_y,
                      Box_z, Box_y, tf01, bf, tf02, tw, Z1, Z2, cover,
                      H, DA01, DA02, DA03, WEIGHT, ND, NUM_NODES, DMAX, DINCR, J, MAX_ITERATIONS, MAX_TOLERANCE):
    op.wipe()
    op.model('basic', '-ndm', 2, '-ndf', 3)  # frame 2D

    # DEFINE BRIDGE COMPOSITE SECTION
    Section_Tag = 1 # Composite Section Tag
    nFibZ = 5 # NUMBER OF FIBER SECTIONIN Z DIRECTION
    nFibY = 70 # NUMBER OF FIBER SECTIONIN Y DIRECTION
    FS = Composite_Bridge_Section_Plot(Section_Tag, nFibZ, nFibY, SideWalk_z, SideWalk_y, SideWalk_C,
                                      Deck_z, Deck_y, Box_z, Box_y, tf01, bf, tf02, tw, DA01, DA02, DA03,
                                      PLOT=0, OUTPUT='PUSH', I=J)
    opsv.fib_sec_list_to_cmds(FS)
    
    # nodal coordinates:
    LL = Length / (NUM_NODES-1) # Lnength of each Beam elements
    for i in range(0, NUM_NODES):
        op.node(i+1, LL*i, 0.0) # node#, X, Y
        #print(i+1, LL*i)
        
    # Constraints -- Boundary Conditions
    op.fix(1, 1, 1, 0) # node DX DY RZ
    op.fix(NUM_NODES, 1, 1, 0)
    
    op.geomTransf('Linear', 1)
    numIntgrPts = 5
    for i in range(0, NUM_NODES-1):
        op.element('nonlinearBeamColumn', i+1, i+1, i+2, numIntgrPts, Section_Tag, 1)
        #print(i+1, i+1, i+2)
        
    # OUTPUT DATA
    op.recorder('Node', '-file', f"{SALAR_DIR}DTH_PUSH_{J}.txt",'-time', '-node', ND, '-dof', 1,2,3, 'disp')# Displacement Time History Node 2
    op.recorder('Node', '-file', f"{SALAR_DIR}BTH_PUSH_01_{J}.txt",'-time', '-node', 1, '-dof', 1,2,3, 'reaction')# Base Shear Time History Node 1
    op.recorder('Node', '-file', f"{SALAR_DIR}BTH_PUSH_11_{J}.txt",'-time', '-node', NUM_NODES, '-dof', 1,2,3, 'reaction')# Base Shear Time History Node 11
    op.recorder('Element', '-file', f"{SALAR_DIR}DEF_PUSH_{J}.txt",'-time', '-ele', ND, 'section', 5, 'deformations')# Curvature Time History for ND
    for i in range(0, NUM_NODES-1):
        op.recorder('Element', '-file', f"{SALAR_DIR}LOCALFORCE_PUSH_{i}_{J}.txt",'-time', '-ele', i+1,'localForce')# Local Force Time History for Each Element 
    
    #defining gravity loads
    op.timeSeries('Linear', 1)
    op.pattern('Plain', 1, 1)
    for i in range(1, NUM_NODES+1):
        op.load(i, 0.0, -WEIGHT, 0.0)
    print('Model Built')    

    NstepGravity = 10
    DGravity = 1 / NstepGravity
    op.integrator('LoadControl', DGravity) # determine the next time step for an analysis
    op.numberer('Plain') # renumber dof's to minimize band-width (optimization), if you want to
    op.system('BandGeneral') # how to store and solve the system of equations in the analysis
    op.constraints('Plain') # how it handles boundary conditions
    op.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 0) # determine if convergence has been achieved at the end of an iteration step
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    op.algorithm('Newton') # use Newton's solution algorithm: updates tangent stiffness at every iteration
    op.analysis('Static') # define type of analysis static or transient
    op.analyze(NstepGravity) # apply gravity

    op.loadConst('-time', 0.0) #maintain constant gravity loads and reset time to zero
    Vload = 1#Weight

    op.timeSeries('Linear', 2)
    op.pattern('Plain', 200, 2)
    op.load(ND, 0.0, Vload, 0.0)

    op.wipeAnalysis()
    op.constraints('Plain')
    op.numberer('Plain')
    op.system('BandGeneral')
    op.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 2)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    op.algorithm('Newton')

    IDctrlDOF = 2   ## INCREMENTAL DISPLACEMENT NODE IN Y DIRECTION

    op.integrator('DisplacementControl', ND, IDctrlDOF, DINCR)
    op.analysis('Static')

    Nsteps =  int(np.abs(DMAX/ DINCR))
    
    #OK = op.analyze(Nsteps)
    #S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS

    # Data storage
    FORCE_S, FORCE_A, MOMENT = [], [], []
    DISP_X, DISP_Y, ROT = [], [], []
    KA, KS, KI, STEP = [], [], [], []
    
    for step in range(Nsteps):
        OK = op.analyze(1)
        S02.ANALYSIS(OK, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        # Record results
        op.reactions()
        S = op.nodeReaction(1, 2) + op.nodeReaction(NUM_NODES, 2) # SHEAR BASE REACTION
        A = op.nodeReaction(1, 1) + op.nodeReaction(NUM_NODES, 1) # AXIAL BASE REACTION
        M = op.nodeReaction(1, 3) + op.nodeReaction(NUM_NODES, 3) # MOMENT BASE REACTION
        #print(rot, M)
        disp_X = op.nodeDisp(ND, 1) # LATERAL DISPLACEMENT IN X FOR NODE IN MIDDLE ON BRIDGE LENGTH
        disp_Y = op.nodeDisp(ND, 2) # LATERAL DISPLACEMENT IN Y FOR NODE IN MIDDLE ON BRIDGE LENGTH
        rot = op.nodeDisp(ND, 3)    # ROTATION IN Z FOR NODE IN MIDDLE ON BRIDGE LENGTH
        FORCE_S.append(S)
        FORCE_A.append(A)
        MOMENT.append(M)
        DISP_X.append(disp_X)
        DISP_Y.append(disp_Y)
        ROT.append(rot)
        KS.append(np.abs(S)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
        KA.append(np.abs(A)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
        KI.append(np.abs(M)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
        STEP.append(step)
        print(step+1, disp_X, S)
    
    print(' \n\nPushover Done.')
    #op.wipe()
    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP  
    
#%%-------------------------------------------------
### ----------------------
###    DYNAMIC FUNCTION
### ----------------------
def DYNAMIC_ANALYSIS(Length, SideWalk_z, SideWalk_y, SideWalk_C, Deck_z, Deck_y,
                      Box_z, Box_y, tf01, bf, tf02, tw, Z1, Z2, cover,
                      H, DA01, DA02, DA03, MASS, WEIGHT, ND, NUM_NODES, DAMPING_RATIO, J, MAX_ITERATIONS, MAX_TOLERANCE, dt, Duration):
    op.wipe()
    op.model('basic', '-ndm', 2, '-ndf', 3)  # frame 2D
    
    # DEFINE BRIDGE COMPOSITE SECTION
    Section_Tag = 1 # Composite Section Tag
    nFibZ = 5 # NUMBER OF FIBER SECTIONIN Z DIRECTION
    nFibY = 70 # NUMBER OF FIBER SECTIONIN Y DIRECTION
    FS = Composite_Bridge_Section_Plot(Section_Tag, nFibZ, nFibY, SideWalk_z, SideWalk_y, SideWalk_C,
                                      Deck_z, Deck_y, Box_z, Box_y, tf01, bf, tf02, tw, DA01, DA02, DA03,
                                      PLOT=0, OUTPUT='DYN', I=J)
    opsv.fib_sec_list_to_cmds(FS)
    
    # Nodal coordinates:
    LL = Length / (NUM_NODES-1) # Lnength of each Beam elements
    for i in range(0, NUM_NODES):
        op.node(i+1, LL*i, 0.0) # node#, X, Y
        #print(i+1, LL*i)
        
    # Boundary Conditions
    op.fix(1, 1, 1, 0) # node DX DY RZ
    op.fix(NUM_NODES, 1, 1, 0)
    
    op.geomTransf('Linear', 1)
    numIntgrPts = 5
    for i in range(0, NUM_NODES-1):
        op.element('nonlinearBeamColumn', i+1, i+1, i+2, numIntgrPts, Section_Tag, 1)
        #print(i+1, i+1, i+2)

    # OUTPUT DATA
    op.recorder('EnvelopeNode','-file', f"{SALAR_DIR}MD_{J}.txt" ,'-time','-node',ND,'-dof',2,'disp');# max. displacements of free node ND
    op.recorder('EnvelopeNode','-file',f"{SALAR_DIR}MV_{J}.txt" ,'-time','-node',ND,'-dof',2,'vel');# max. vel of free node ND
    op.recorder('EnvelopeNode','-file', f"{SALAR_DIR}MA_{J}.txt" ,'-time','-node',ND,'-dof',2,'accel');# max. accel of free node ND	
    op.recorder('Node', '-file', f"{SALAR_DIR}DTH_DYN_{J}.txt",'-time', '-node', ND, '-dof', 1,2,3, 'disp')# Displacement Time History Node ND
    op.recorder('Node', '-file', f"{SALAR_DIR}VTH_DYN_{J}.txt",'-time', '-node', ND, '-dof', 1,2,3, 'vel')# Velocity Time History Node ND
    op.recorder('Node', '-file', f"{SALAR_DIR}ATH_DYN_{J}.txt",'-time', '-node', ND, '-dof', 1,2,3, 'accel')# Acceleration Time History Node ND
    op.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_01_{J}.txt",'-time', '-node', 1, '-dof', 1,2,3, 'reaction')# Base Shear Time History Node 1
    op.recorder('Node', '-file', f"{SALAR_DIR}BTH_DYN_11_{J}.txt",'-time', '-node', NUM_NODES, '-dof', 1,2,3, 'reaction')# Base Shear Time History Node 11
    op.recorder('Element', '-file', f"{SALAR_DIR}DEF_DYN_{J}.txt",'-time', '-ele', ND, 'section', 5, 'deformations')# Curvature Time History for ND
    for i in range(0, NUM_NODES-1):
        op.recorder('Element', '-file', f"{SALAR_DIR}LOCALFORCE_DYN_{i}_{J}.txt",'-time', '-ele', i+1,'localForce')# Local Force Time History for Each Element 
     
    # Define gravity loads
    op.timeSeries('Linear', 1)
    op.pattern('Plain', 1, 1)
    # Define Weight
    for j in range(1, NUM_NODES+1): 
        op.load(j, 0.0, -WEIGHT, 0.0)
        
    # Define Mass
    for j in range(1, NUM_NODES+1): 
        op.mass(j, 0.0, MASS, 0.0)    
    print('Model Built')    

    NstepGravity = 10
    DGravity = 1/NstepGravity
    op.integrator('LoadControl', DGravity) # determine the next time step for an analysis
    op.numberer('Plain') # renumber dof's to minimize band-width (optimization), if you want to
    op.system('BandGeneral') # how to store and solve the system of equations in the analysis
    op.constraints('Plain') # how it handles boundary conditions
    op.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 0) # determine if convergence has been achieved at the end of an iteration step
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    op.algorithm('Newton') # use Newton's solution algorithm: updates tangent stiffness at every iteration
    op.analysis('Static') # define type of analysis static or transient
    op.analyze(NstepGravity) # apply gravity

    op.loadConst('-time', 0.0) #maintain constant gravity loads and reset time to zero

    # Calculate Rayleigh damping factors
    #Lambda01 = op.eigen('-fullGenLapack', 2)  # eigenvalue mode 2
    Lambda01 = op.eigen('-genBandArpack', 2) # eigenvalue mode 2
    Omega01 = np.power(max(Lambda01), 0.5)
    Omega02 = np.power(min(Lambda01), 0.5)
    a0 = (2 * Omega01 * Omega02 * DAMPING_RATIO) / (Omega01 + Omega02) # c = a0 * m : Mass-proportional damping
    a1 = (DAMPING_RATIO * 2) / (Omega01 + Omega02)                     # c = a1 * k : Stiffness-proportional damping
    # Apply Rayleigh damping
    op.rayleigh(a0, a1, 0, 0)   # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    #op.rayleigh(0, 0, 2 * DR * Omega01, 0) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/reyleigh.html
    PERIOD_01 = (np.pi * 2) / Omega01 # Structure First Period
    PERIOD_02 = (np.pi * 2) / Omega02 # Structure Second Period
    #print('Structure First Period:  ', PERIOD_01)
    #print('Structure Second Period: ', PERIOD_02) 

    #applying Dynamic Ground motion analysis
    GMfact = 9810    # standard acceleration of gravity or standard acceleration
    SSF_X = 0.01     # Seismic Acceleration Scale Factor in X Direction
    SSF_Y = 0.01     # Seismic Acceleration Scale Factor in Y Direction
    iv0_X = 0.00005  # [mm/s] Initial velocity applied to the node  in X Direction
    iv0_Y = 0.00005  # [mm/s] Initial velocity applied to the node  in Y Direction
    st_iv0 = 0.0     # [s] Initial velocity applied starting time
    
    SEISMIC_TAG_01 = 100
    op.timeSeries('Path', SEISMIC_TAG_01, '-dt', dt, '-filePath', 'Ground_Acceleration_X.txt', '-factor', GMfact, '-startTime', st_iv0) # SEISMIC-X
    # Define load patterns
    # pattern UniformExcitation $patternTag $dof -accel $tsTag <-vel0 $vel0> <-fact $cFact>
    op.pattern('UniformExcitation', SEISMIC_TAG_01, 1, '-accel', SEISMIC_TAG_01, '-vel0', iv0_X, '-fact', SSF_X) # SEISMIC-X 
    SEISMIC_TAG_02 = 200
    op.timeSeries('Path', SEISMIC_TAG_02, '-dt', dt, '-filePath', 'Ground_Acceleration_Y.txt', '-factor', GMfact) # SEISMIC-Z
    op.pattern('UniformExcitation', SEISMIC_TAG_02, 2, '-accel', SEISMIC_TAG_02, '-vel0', iv0_Y, '-fact', SSF_Y)  # SEISMIC-Z 

    op.wipeAnalysis()
    op.constraints('Plain')
    op.numberer('Plain')
    op.system('BandGeneral')
    op.test('NormDispIncr', MAX_TOLERANCE, MAX_ITERATIONS, 2)
    # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/normDispIncr.html
    op.algorithm('Newton')

    NewmarkGamma = 0.5; NewmarkBeta = 0.25
    op.integrator('Newmark', NewmarkGamma, NewmarkBeta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/newmark.html
    #alpha = 1;gamma=1.5-alpha; beta=((2-alpha)**2)/4;
    #op.integrator('HHT', alpha, gamma, beta) # INFO LINK: https://openseespydoc.readthedocs.io/en/latest/src/hht.html
    op.analysis('Transient')

    # Perform Dynamic analysis
    #n_steps =  int(np.abs(Duration/ dt))# Analysis Steps
    
    #OK = op.analyze(n_steps, dt)
    #S02.ANALYSIS(OK, n_steps, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
  
    # Data storage
    FORCE_S, FORCE_A, MOMENT = [], [], []
    DISP_X, DISP_Y, ROT = [], [], []
    KA, KS, KI, STEP = [], [], [], []
    time = []
    displacement = []
    velocity_X, velocity_Y = [], []
    acceleration_X, acceleration_Y = [], []
    
    stable = 0
    current_time = 0.0
       
    while stable == 0 and current_time < Duration:
        stable = op.analyze(1, dt)
        S02.ANALYSIS(stable, 1, MAX_TOLERANCE, MAX_ITERATIONS) # CHECK THE ANALYSIS
        current_time = op.getTime()
        time.append(current_time)
        # Record results
        op.reactions()
        S = op.nodeReaction(1, 2) + op.nodeReaction(NUM_NODES, 2) # SHEAR BASE REACTION
        A = op.nodeReaction(1, 1) + op.nodeReaction(NUM_NODES, 1) # AXIAL BASE REACTION
        M = op.nodeReaction(1, 3) + op.nodeReaction(NUM_NODES, 3) # MOMENT BASE REACTION
        #print(rot, M)
        disp_X = op.nodeDisp(ND, 1) # LATERAL DISPLACEMENT IN X FOR NODE IN MIDDLE ON BRIDGE LENGTH
        disp_Y = op.nodeDisp(ND, 2) # LATERAL DISPLACEMENT IN Y FOR NODE IN MIDDLE ON BRIDGE LENGTH
        rot = op.nodeDisp(ND, 3)    # ROTATION IN Z FOR NODE IN MIDDLE ON BRIDGE LENGTH
        velocity_X.append(op.nodeVel(ND, 1))       # LATERAL VELOCITY IN X FOR NODE IN MIDDLE ON BRIDGE LENGTH
        acceleration_X.append(op.nodeAccel(ND, 1)) # LATERAL ACCELERATION IN X FOR NODE IN MIDDLE ON BRIDGE LENGTH
        velocity_Y.append(op.nodeVel(ND, 2))       # LATERAL VELOCITY IN Y FOR NODE IN MIDDLE ON BRIDGE LENGTH
        acceleration_Y.append(op.nodeAccel(ND, 2)) # LATERAL ACCELERATION IN Y FOR NODE IN MIDDLE ON BRIDGE LENGTH
        FORCE_S.append(S)
        FORCE_A.append(A)
        MOMENT.append(M)
        DISP_X.append(disp_X)
        DISP_Y.append(disp_Y)
        ROT.append(rot)
        KS.append(np.abs(S)/np.abs(disp_X)) # LATERAL STIFFNESS IN X
        KA.append(np.abs(A)/np.abs(disp_Y)) # LATERAL STIFFNESS IN Y
        KI.append(np.abs(M)/np.abs(rot))    # ROTATIONAL STIFFNESS IN Z
        print(current_time, disp_X, S)
        
    # Calculating Damping Ratio Using Logarithmic Decrement Analysis 
    displacement = np.array(DISP_Y)
    peaks = np.array([displacement[i] for i in range(1, len(displacement)-1) if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]])
    # Natural logarithm
    delta = np.log(peaks[:-1] / peaks[1:])    
        
    print(' \n\nDynamic Done.')
    #op.wipe()
    return FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, time, velocity_X, velocity_Y, acceleration_X, acceleration_Y, PERIOD_01, PERIOD_02, delta
#%%-------------------------------------------------
# COMPOSITE BRIDGE SECTION
Length = 30000    # [mm] Bridge Length
SideWalk_z = 1000 # [mm] SideWalk Width
SideWalk_y = 150  # [mm] SideWalk Height
SideWalk_C = 250  # [mm] SideWalk to deck Connection
Deck_z = 11500/2  # [mm] Bridge Deck Width Concrete Top Flange
Deck_y = 300      # [mm] Bridge Thickness Concrete Top Flange
Box_z = 500       # [mm] Width Concrete Box Under Deck
Box_y = 200       # [mm] Thickness Concrete Box Under Deck
tf01 = 30         # [mm] Top Flange Plate Thickness Beam
bf = Box_z        # [mm] Top Flange Plate Thickness Beam (BOX WIDTH, TOP AND BOTTTOM PLATE WIDTH IS THE SAME AN EQUAL)
tf02 = 30         # [mm] Bottom Flange Plate Thickness Beam
tw = 16           # [mm] Top Web Plate Thickness Beam 01

Z1 = 500          # [mm]
Z2 = 2000         # [mm]
cover = 50        # [mm] Concrete Cover
H = 1800          # [mm] I Section Beam Height

DA01 = 18         # [mm] Rebar Diameter for Sidewalk
DA02 = 25         # [mm] Rebar Diameter for Deck
DA03 = 16         # [mm] Rebar Diameter for Box
#%%-------------------------------------------------
#import matplotlib.pyplot as plt
#import numpy as np
#from mpl_toolkits.mplot3d import Axes3D

# Define beam length (example value)
Length = 30000  # mm

# Number of nodes
num_nodes = 11

# Generate node positions
x = np.linspace(0, Length, num_nodes)
y = np.zeros(num_nodes)
z = np.zeros(num_nodes)

def MATPLOT_2D_BEAM(x, y, num_nodes, title):
    fig, ax = plt.subplots(figsize=(15, 3))
    
    # Beam line
    ax.plot([x[0], x[-1]], [0, 0], 'k-', lw=3, label='Beam')
    
    # Nodes (only plot markers, no labels to reduce clutter)
    ax.plot(x, y, 'bo', markersize=5, label='Nodes')
    
    # Support points (triangles)
    ax.plot([x[0], x[-1]], [0, 0], 'r^', markersize=12, label='Supports')
    
    # Label only first and last nodes to avoid clutter
    ax.text(x[0], 0.02*Length, 'Node 1', ha='center', fontsize=9)
    ax.text(x[-1], 0.02*Length, f'Node {num_nodes}', ha='center', fontsize=9)
    
    # Add distance markers every 10% of length
    for pos in np.linspace(0, Length, 6):
        ax.plot([pos, pos], [-0.01*Length, 0.01*Length], 'k-', lw=0.5)
        ax.text(pos, -0.03*Length, f'{pos:.0f} mm', ha='center', fontsize=8)

    ax.set_title(title, fontsize=14)
    ax.set_xlabel('Length (mm)', fontsize=12)
    ax.set_ylim(-0.1*Length, 0.1*Length)
    ax.legend(loc='upper right')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

def MATPLOT_3D_BEAM(x, y, z, num_nodes, title):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Beam line
    ax.plot(x, y, z, 'k-', lw=3, label='Beam')
    
    # Nodes
    ax.scatter(x, y, z, s=20, c='b', label='Nodes')
    
    # Support points
    ax.scatter([x[0], x[-1]], [0, 0], [0, 0], s=100, c='r', marker='^', label='Supports')
    
    # Label key nodes
    ax.text(x[0], 0, 0, 'Node 1', fontsize=9)
    ax.text(x[-1], 0, 0, f'Node {num_nodes}', fontsize=9)
    
    # Set axis properties
    ax.set_title(title, fontsize=14)
    ax.set_xlabel('Length (mm)', fontsize=12)
    ax.set_ylabel('Y Position', fontsize=12)
    ax.set_zlabel('Z Position', fontsize=12)
    
    # Set equal aspect ratio
    max_dim = max(Length, Length*0.1)
    ax.set_box_aspect([max_dim, max_dim/10, max_dim/10])
    
    # Add grid and legend
    ax.grid(True)
    ax.legend()
    
    # Set viewing angle
    ax.view_init(elev=30, azim=-45)
    plt.tight_layout()
    plt.show()

# Generate plots
MATPLOT_2D_BEAM(x, y, num_nodes, f'2D Simply Supported Beam with {num_nodes} Nodes')
MATPLOT_3D_BEAM(x, y, z, num_nodes, f'3D Simply Supported Beam with {num_nodes} Nodes')
#%%-------------------------------------------------
Section_Tag = 0 # Composite Section Tag
nFibZ = 1
nFibY = 1
PLOT = 1
OUTPUT = 'PUSH'
I = 0
d = Composite_Bridge_Section_Plot(Section_Tag, nFibZ, nFibY, SideWalk_z, SideWalk_y, SideWalk_C, Deck_z,
                                 Deck_y, Box_z, Box_y, tf01, bf, tf02, tw, DA01, DA02, DA03, PLOT, OUTPUT, I)
#%%-------------------------------------------------
# RUN MOMENT-CURVATURE ANLYSIS
PX = 0.0                   # [N] Axial Section Force
PY, MZ = 0.0, 0.0
DR = 26.5                  # Ductility ratio for moment curvature
numIncr = 200              # Number of analysis increments
MAX_ITERATIONS = 20000     # Convergence iteration for test
MAX_TOLERANCE = 1.0e-8     # Convergence tolerance for test
MC_ANALYSIS(Section_Tag, PX, PY, MZ, DR, numIncr, MAX_ITERATIONS, MAX_TOLERANCE, J=0)
#%%-------------------------------------------------
CUR = S05.OUTPUT_SECOND_COLUMN(FOLDER_NAME,'CUR', 1, 0, 2)
MOM = S05.OUTPUT_SECOND_COLUMN(FOLDER_NAME,'MOM', 1, 0, 2)
xxc, yyc, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = S05.BILNEAR_CURVE(CUR, -MOM, 2)
xxc = np.abs(xxc); yyc = np.abs(yyc); # ABSOLUTE VALUE
XLABEL = 'Curvature'
YLABEL = 'Moment'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = 'Moment and Curvature Analysis'
COLOR = 'black'
S05.PLOT_2D(CUR, -MOM, xxc, yyc, xxc, yyc, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 

print('###         SECTION PARAMETERS BASED ON ANALYSIS      ###')
print('=========================================================')
print(f' Section Elastic Stiffness :      {Elastic_ST:.2f}')
print(f' Section Plastic Stiffness :      {Plastic_ST:.2f}')
print(f' Section Tangent Stiffness :      {Tangent_ST:.2f}')
print(f' Section Ductility Ratio :        {Ductility_Rito:.2f}')
print(f' Section Over Strength Factor:    {Over_Strength_Factor:.2f}')
print(f' Section Yield Curvature:         {xxc[1]:.6f}')
print(f' Section Ultimate Curvature:      {xxc[2]:.6f}')
#%%-------------------------------------------------
### -----------------------------------------------------------
###   Composite Bridge Superstructure Ductility Damage Index 
### -----------------------------------------------------------

# Define parameters (units: mm, N)
# ------------------------------------------
g = 9810                # [mm/s^2] Acceleration due to Gravity
DMAX = -70              # [mm] Max. Pushover Incremental Displacement
DINCR = -0.01           # [mm] Pushover Increment 
ND = 6                  # NODE NUMBER FOR INCREMENTAL DISPLACEMENT
NUM_NODES = 11          # NUMBER OF NODES IN BEAM ELEMENT

DAMPING_RATIO = 0.05	# damping ratio
dt = 0.01               # [s] Time increment
Duration = 15.0         # [s] Total time duration 

MAX_ITERATIONS = 20000  # Convergence iteration for test
MAX_TOLERANCE = 1.0e-8  # Convergence tolerance for test

st = 5 # Sleep Time 

starttime = TI.process_time()

t = TI.localtime()
current_time = TI.strftime("%H:%M:%S", t)
print(f"Current time (HH:MM:SS): {current_time}\n\n")

Massef = 95000.0               # [kg] Total Mass of Structure
MASS_NODE = Massef / NUM_NODES # Mass of each Node (Truss has 62 Nodes)
WEIGHT_NODE = MASS_NODE * g    # Weight of each Node
#%%------------------------------------------------- 
#%% RUN NONLINEAR STATIC ANALYSIS:     
DATA = PUSHOVER_ANALYSIS(Length, SideWalk_z, SideWalk_y, SideWalk_C, Deck_z, Deck_y,
                      Box_z, Box_y, tf01, bf, tf02, tw, Z1, Z2, cover,
                      H, DA01, DA02, DA03, WEIGHT_NODE, ND, NUM_NODES, DMAX, DINCR, 0, MAX_ITERATIONS, MAX_TOLERANCE)
FORCE_Sp, FORCE_Ap, MOMENTp, DISP_Xp, DISP_Yp, ROTp, KAp, KSp, KIp, STEPp = DATA

# %% Plot 2D Frame Shapes for Nonlinear Static Analysis
S04.PLOT_2D_FRAME(deformed_scale=10)  # Adjust scale factor as needed 
#%%-------------------------------------------------  
#%% RUN NONLINEAR DYNAMIC ANALYSIS:  
DATA = DYNAMIC_ANALYSIS(Length, SideWalk_z, SideWalk_y, SideWalk_C, Deck_z, Deck_y,
                         Box_z, Box_y, tf01, bf, tf02, tw, Z1, Z2, cover,
                         H, DA01, DA02, DA03, MASS_NODE, WEIGHT_NODE, ND, NUM_NODES, DAMPING_RATIO, 0, MAX_ITERATIONS, MAX_TOLERANCE, dt, Duration)
FORCE_Sd, FORCE_Ad, MOMENTd, DISP_Xd, DISP_Yd, ROTd, KAd, KSd, KId, timed, velocity_Xd, velocity_Yd, acceleration_Xd, acceleration_Yd, PERIOD_01d, PERIOD_02d, deltad = DATA

print('Structure First Period:  ', PERIOD_01d, ' (s)')
print('Structure Second Period: ', PERIOD_02d, ' (s)') 

# %% Plot 2D Frame Shapes for Nonlinear Dynamic Analysis
S04.PLOT_2D_FRAME(deformed_scale=100)  # Adjust scale factor as needed 
#%%-------------------------------------------------  
TI.sleep(st);# Sleep time
#%% STRUCTURE DUCTILITY DAMAGE INDEX
dispD = S05.OUTPUT_SECOND_COLUMN(FOLDER_NAME,'DTH_DYN', 2, 0, 2) # Reading Disp from Text file - DYNAMIC - NODE 16
dispP = S05.OUTPUT_SECOND_COLUMN(FOLDER_NAME,'DTH_PUSH', 2, 0, 2) # Reading Disp from Text file - PUSHOVER
base01 = S05.OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_PUSH_01', 2, 0, 2) # Reading base shear from Text file - PUSHOVER - NODE 1
base02 = S05.OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_PUSH_11', 2, 0, 2) # Reading base shear from Text file - PUSHOVER - NODE 11
dispP = np.abs(dispP); baseP = np.abs(base01 + base02); # ABSOLUTE VALUE
xx, yy, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = S05.BILNEAR_CURVE(dispP, baseP, 30)
xx = np.abs(xx); yy = np.abs(yy) # ABSOLUTE VALUE
demand_disp = np.max(np.abs(dispD))# DIPLACEMENT - DYNAMIC
# WHEN SOFTENING HAPPENED THIS FORMULA HAS ERROR BEACAUSE OF LARGE DISPLACEMENT
DI = (demand_disp - xx[1]) / (xx[2] - xx[1])  # DUCTIITY DAMGE INDEX
#DI = (demand_disp - xx[1]) / (DMAX - xx[1])

t = TI.localtime()
current_time = TI.strftime("%H:%M:%S", t)
print(f"Current time (HH:MM:SS): {current_time}\n\n")

totaltime = TI.process_time() - starttime
print(f'\nTotal time (s): {totaltime:.4f} \n\n') 
#%%------------------------------------------------------------------------------
plt.figure(1, figsize=(12, 8))
plt.plot(MOMENTp, FORCE_Ap, color='black')
plt.plot(MOMENTd, FORCE_Ad, color='purple')
#plt.scatter(MOMENTp, FORCE_Ap, color='black', linewidth=2)
#plt.scatter(MOMENTd, FORCE_Ad, color='purple', linewidth=2)
plt.title('P-M Interaction')
plt.ylabel('Axial Force [N]')
plt.xlabel('Bending Moment [N.mm]')
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(2, figsize=(12, 8))
plt.plot(DISP_Xp, FORCE_Sp, color='green', linewidth=2)
plt.plot(DISP_Xd, FORCE_Sd, color='lime', linewidth=2)
#plt.scatter(DISP_Xp, FORCE_Sp, color='green', linewidth=2)
#plt.scatter(DISP_Xd, FORCE_Sd, color='lime', linewidth=2)
plt.title('SHEAR FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Shear Force [N]')
plt.xlabel('Displacement  in X [mm]')
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(3, figsize=(12, 8))
plt.plot(DISP_Yp, FORCE_Ap, color='purple', linewidth=2)
plt.plot(DISP_Yd, FORCE_Ad, color='green', linewidth=2)
#plt.scatter(DISP_Yp, FORCE_Ap, color='purple', linewidth=2)
#plt.scatter(DISP_Yd, FORCE_Ad, color='green', linewidth=2)
plt.title('AXIAL FORCE-DISPLACEMENT DIAGRAM')
plt.ylabel('Axial Force [N]')
plt.xlabel('Displacement in Y [mm]')
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(4, figsize=(12, 8))
plt.plot(ROTp, MOMENTp, color='red', linewidth=2)
plt.plot(ROTd, MOMENTd, color='pink', linewidth=2)
#plt.scatter(ROTp, MOMENTp, color='red', linewidth=2)
#plt.scatter(ROTd, MOMENTd, color='pink', linewidth=2)
plt.title('MOMENT-ROTATION DIAGRAM')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Rotation [rad]')
plt.legend(['PUSHOVER', 'DYNAMIC'])

plt.grid()
plt.show()

plt.figure(5, figsize=(12, 8))
#plt.plot(KIp, KSp, color='black', linewidth=2)
#plt.plot(KId, KSd, color='grey', linewidth=2)
plt.scatter(KIp, KSp, color='black', linewidth=2)
plt.scatter(KId, KSd, color='grey', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in X Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(6, figsize=(12, 8))
#plt.plot(KIp, KAp, color='black', linewidth=2)
#plt.plot(KId, KAd, color='grey', linewidth=2)
plt.scatter(KIp, KAp, color='black', linewidth=2)
plt.scatter(KId, KAd, color='grey', linewidth=2)
plt.title('ROTATIONAL STIFFNESS-LATERAL STIFFNESS DIAGRAM')
plt.ylabel('Rotational Stiffness [N.mm/Rad]')
plt.xlabel('Lateral Stiffness in Y Dir. [N/mm]')
plt.semilogx()
plt.semilogy()
plt.legend(['PUSHOVER', 'DYNAMIC'])
plt.grid()
plt.show()

plt.figure(7, figsize=(12, 8))
plt.plot(timed, FORCE_Ad, color='brown', linewidth=2)
#plt.scatter(timed, FORCE_Ad, color='brown', linewidth=2)
plt.title('Axial Force During the Analysis')
plt.ylabel('Axial Force [N]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(8, figsize=(12, 8))
plt.plot(timed, FORCE_Sd, color='purple', linewidth=2)
#plt.scatter(timed, FORCE_Sd, color='purple', linewidth=2)
plt.title('Shear Force During the Analysis')
plt.ylabel('Shear Force [N]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(9, figsize=(12, 8))
plt.plot(timed, MOMENTd, color='green', linewidth=2)
#plt.scatter(timed, MOMENTd, color='green', linewidth=2)
plt.title('Moment During the Analysis')
plt.ylabel('Moment [kN.mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(10, figsize=(12, 8))
plt.plot(timed, DISP_Xd, color='brown', linewidth=2)
#plt.scatter(timed, DISP_Xd, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - X [mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(11, figsize=(12, 8))
plt.plot(timed, DISP_Yd, color='blue', linewidth=2)
#plt.scatter(timed, DISP_Xd, color='brown', linewidth=2)
plt.title('Displacement During the Analysis')
plt.ylabel('Displacement - Y [mm]')
plt.xlabel('Times')
plt.grid()
plt.show()

plt.figure(12, figsize=(12, 8))
plt.plot(timed, ROTd, color='black', linewidth=2)
#plt.scatter(timed, ROTd, color='black', linewidth=2)
plt.title('Rotation During the Analysis')
plt.ylabel('Rotation [rad]')
plt.xlabel('Times')
plt.grid()
plt.show()
#%%------------------------------------------------------------------------------  
# EXACT SOLUTION:
from scipy.optimize import fsolve

# Define the equation for natural logarithm of this ratio, called the logarithmic decrement, we denote by δ
def EQUATION(x, delta):
    if np.any(x == 0):  # Avoid division by zero
        return np.inf  
        
    # Calculate the value of the equation
    A = x**2 - 1 + ((2 * np.pi * x) / np.mean(delta)) ** 2
    #print(f"x: {x}, A: {A}")  # Debugging output
    # Return the difference (for root finding)
    return A
      
# Initial guess for root(s)
x0 = 1  # Intial Guess for Damping Ratio
# Solve for x
solution = fsolve(EQUATION, x0, args=(deltad))
print(f"Exact Damping Ratio: {solution[0]:.8e}")
#%%------------------------------------------------------------------------------
# Compute the Cumulative Maximum Absolute Value of Last Analysis Data
def MAX_ABS(X):
    import numpy as np
    X = np.asarray(X)  # Convert input to a numpy array for faster operations
    X_MAX = np.zeros_like(X)  # Initialize an array to store cumulative max values
    X_MAX[0] = np.abs(X[0])  # Set the first value

    # Compute cumulative maximum absolute values
    for i in range(1, len(X)):
        X_MAX[i] = max(X_MAX[i-1], np.abs(X[i]))
    
    return X_MAX  

DISP_ZX = MAX_ABS(DISP_Xd)  
DISP_ZY = MAX_ABS(DISP_Yd) 
VELO_Z = MAX_ABS(velocity_Yd) 
ACCE_Z = MAX_ABS(acceleration_Yd) 
BASE_Z = MAX_ABS(FORCE_Sd) 

plt.figure(1, figsize=(8, 6))
plt.plot(timed, DISP_Xd, color='blue', linewidth=2)
plt.plot(timed, DISP_ZX, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in X [mm]')
plt.title(f'Time vs Displacement  - MAX. ABS: {DISP_ZY[-1]}')
plt.grid()
plt.show()

plt.figure(2, figsize=(8, 6))
plt.plot(timed, DISP_Yd, color='blue', linewidth=2)
plt.plot(timed, DISP_ZY, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Displacement in Y [mm]')
plt.title(f'Time vs Displacement - MAX. ABS: {DISP_ZX[-1]}')
plt.grid()
plt.show()

plt.figure(3, figsize=(8, 6))
plt.plot(timed, velocity_Yd, color='blue', linewidth=2)
plt.plot(timed, VELO_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Velocity in Y [mm/s]')
plt.title(f'Time vs Velocity - MAX. ABS: {VELO_Z[-1]}')
plt.grid()
plt.show()

plt.figure(4, figsize=(8, 6))
plt.plot(timed, acceleration_Yd, color='blue', linewidth=2)
plt.plot(timed, ACCE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Acceleration in Y [mm/s^2]')
plt.title(f'Time vs Acceleration - MAX. ABS: {ACCE_Z[-1]}')
plt.grid()
plt.show()

plt.figure(5, figsize=(8, 6))
plt.plot(timed, FORCE_Sd, color='blue', linewidth=2)
plt.plot(timed, BASE_Z, color='red', linewidth=2)
plt.xlabel('Time [s]')
plt.ylabel('Base-reaction [N]')
plt.title(f'Time vs Base-reaction - MAX. ABS: {BASE_Z[-1]}')
plt.grid()
plt.show()

#%%------------------------------------------------------------------------------
import BILINEAR_CURVE as BC

# --------------------------------------
#  Plot Base Shear-Displacement Analysis 
# --------------------------------------
XX = np.abs(DISP_Yp); YY = np.abs(FORCE_Sp); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in X [mm]'
YLABEL = 'Base-Axial Reaction [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Base Shear-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
BC.PLOT_2D(np.abs(DISP_Yp), np.abs(FORCE_Sp), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
#print(f'\t\t Ductility Ratio: {Y[2]/Y[1]:.4f}')

# Calculate Over Strength Coefficient (Ω0)
Omega_0 = Y[2] / Y[1]
# Calculate Displacement Ductility Ratio (μ)
mu = X[2] / X[1]
# Calculate Ductility Coefficient (Rμ)
#R_mu = 1
#R_mu = (2 * mu - 1) ** 0.5
R_mu = mu
# Calculate Structural Behavior Coefficient (R)
R = Omega_0 * R_mu
print(f'Over Strength Coefficient (Ω0):      {Omega_0:.4f}')
print(f'Displacement Ductility Ratio (μ):    {mu:.4f}')
print(f'Ductility Coefficient (Rμ):          {R_mu:.4f}')
print(f'Structural Behavior Coefficient (R): {R:.4f}')
Dd = np.max(np.abs(DISP_Yd))
DIx = (Dd - X[1]) /(X[2] - X[1])
print(f'Structural Ductility Damage Index in X Direction: {DIx:.4f}')
# ---------------------------------------
#  Plot Base Axial-Displacement Analysis
# ---------------------------------------
XX = np.abs(DISP_Xp); YY = np.abs(FORCE_Ap); # ABSOLUTE VALUE
SLOPE_NODE = 10

DATA = BC.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

XLABEL = 'Displacement in Y [mm]'
YLABEL = 'Base-Axial Reaction [N]'
LEGEND01 = 'Curve'
LEGEND02 = 'Bilinear Fitted'
LEGEND03 = 'Undefined'
TITLE = f'Base Axial-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}'
COLOR = 'black'
BC.PLOT_2D(np.abs(DISP_Xp), np.abs(FORCE_Ap), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
#print(f'\t\t Ductility Ratio: {YY[2]/YY[1]:.4f}')

# Calculate Over Strength Coefficient (Ω0)
Omega_0 = Y[2] / Y[1]
# Calculate Displacement Ductility Ratio (μ)
mu = X[2] / X[1]
# Calculate Ductility Coefficient (Rμ)
#R_mu = 1
#R_mu = (2 * mu - 1) ** 0.5
R_mu = mu
# Calculate Structural Behavior Coefficient (R)
R = Omega_0 * R_mu
print(f'Over Strength Coefficient (Ω0):      {Omega_0:.4f}')
print(f'Displacement Ductility Ratio (μ):    {mu:.4f}')
print(f'Ductility Coefficient (Rμ):          {R_mu:.4f}')
print(f'Structural Behavior Coefficient (R): {R:.4f}')
Dd = np.max(np.abs(DISP_Xd))
DIy = (Dd - X[1]) /(X[2] - X[1])
print(f'Structural Ductility Damage Index in Y Direction: {DIy:.4f}')
#%%-------------------------------------------------        
print('###      STRUCTURAL PARAMETERS BASED ON ANALYSIS      ###')
print('=========================================================')
print(f' Structure Elastic Stiffness :      {Elastic_ST:.2f}')
print(f' Structure Plastic Stiffness :      {Plastic_ST:.2f}')
print(f' Structure Tangent Stiffness :      {Tangent_ST:.2f}')
print(f' Structure Ductility Ratio :        {Ductility_Rito:.2f}')
print(f' Structure Over Strength Factor:    {Over_Strength_Factor:.2f}')
print(f' Structure Yield Displacement:      {xx[1]:.2f}')
print(f' Structure Ultimate Displacement:   {xx[2]:.2f}')
print(f' Structure Demand Displacement:     {demand_disp:.2f}')
print(f' Structure Ductility Damage index:  {100* DI:.2f} %')
dispD = S05.OUTPUT_SECOND_COLUMN(FOLDER_NAME,'DTH_DYN', 2, 0, 2) # Reading Disp from Text file - DYNAMIC - NODE 16
dispP = S05.OUTPUT_SECOND_COLUMN(FOLDER_NAME,'DTH_PUSH', 2, 0, 2) # Reading Disp from Text file - PUSHOVER
base01 = S05.OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_PUSH_01', 2, 0, 2) # Reading base shear from Text file - PUSHOVER - NODE 1
base02 = S05.OUTPUT_SECOND_COLUMN(FOLDER_NAME,'BTH_PUSH_11', 2, 0, 2) # Reading base shear from Text file - PUSHOVER - NODE 11
dispP = np.abs(dispP); baseP = np.abs(base01 + base02); # ABSOLUTE VALUE
xx, yy, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = S05.BILNEAR_CURVE(dispP, baseP, 30)
xx = np.abs(xx); yy = np.abs(yy) # ABSOLUTE VALUE
demand_disp = np.max(np.abs(dispD))# DIPLACEMENT - DYNAMIC
XLABEL = 'DISPLACEMENT DOF (17)'
YLABEL = 'BASE SHEAR DOF (2) + DOF (32)'
TITLE = f'DISPLACEMENT BASE-SHEAR CURVE FOR DYNAMIC AND PUSHOVER ANALYSIS  - DUCTILITY DAMAGE INDEX: {100* DI:.2f} %'
LEGEND01 = 'PUSHOVER'
LEGEND02 = 'PUSHOVER BILINEAR FITTED'
LEGEND03 = 'DYNAMIC'
COLOR = 'blue'
S05.PLOT_2D(dispP, baseP, xx, yy, xx, yy, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR, Z=2)
#%%-------------------------------------------------
# EXPORT DATA TO EXCEL
import pandas as pd
DATA_TOTAL = { 
    'DISP_X_P': DISP_Xp, # P: NONLINEAR STATIC ANALYSIS
    'DISP_X_D': DISP_Xd, # D: NONLINEAR DYNAMIC ANALYSIS
    'DISP_Y_P': DISP_Yp,
    'DISP_Y_D': DISP_Yd,
    'ROTATION_P': ROTp,
    'ROTATION_D': ROTd,
    'VELO_D': velocity_Xd,
    'ACCEL_D': acceleration_Xd,
    'AXIAL_FORCE_P': FORCE_Ap,
    'AXIAL_FORCE_D': FORCE_Ad,
    'SHEAR_FORCE_P': FORCE_Sp,
    'SHEAR_FORCE_D': FORCE_Sd,
    'MOMENT_P': MOMENTp,
    'MOMENT_D': MOMENTd,
    'AXIAL_RIGIDITY_P': np.abs(FORCE_Ap),
    'AXIAL_RIGIDITY_D': np.abs(FORCE_Ad),
    'ROTATIONAL_ST_P': KIp,
    'ROTATIONAL_ST_D': KId,
    'LATERAL_ST_Y_P': KAp,
    'LATERAL_ST_Y_D': KAd,
    'LATERAL_ST_X_P': KSp,
    'LATERAL_ST_X_D': KSd,
}
# Convert to DataFrame
results_df = pd.DataFrame(DATA_TOTAL)
# Export the DataFrame to an Excel file
results_df.to_excel('SINGLE_SPAN_SIMPLY_SUPPORTED_COMPOSITE_BRIDGE_RESULTS.xlsx', index=False) 
#%%------------------------------------------------------------------------------
# Print out the state of nodes
op.printModel("node",1, ND, NUM_NODES)
# Print out the state of element 1 , 2 and 3
op.printModel("ele", 1, 2 , 3)
# Print the Model
#printModel()
op.printModel("-JSON", "-file", "SINGLE_SPAN_SIMPLY_SUPPORTED_COMPOSITE_BRIDGE.json")
#%%-------------------------------------------------------------------------------
