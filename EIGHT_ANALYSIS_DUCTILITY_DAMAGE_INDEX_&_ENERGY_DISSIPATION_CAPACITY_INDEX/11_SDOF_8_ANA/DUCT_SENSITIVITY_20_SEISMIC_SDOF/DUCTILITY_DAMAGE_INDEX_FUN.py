def PLOT_2D(X, Y, Xfit, Yfit, X2, Y2, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR, Z):
    import matplotlib.pyplot as plt
    plt.figure(figsize=(12, 8))
    if Z == 1:
        # Plot 1 line
        plt.plot(X, Y,color=COLOR)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        plt.title(TITLE)
        plt.grid(True)
        plt.show()
    if Z == 2:
        # Plot 2 lines
        plt.plot(X, Y, Xfit, Yfit, 'r--', linewidth=3)
        plt.title(TITLE)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        plt.legend([LEGEND01, LEGEND02], loc='lower right')
        plt.grid(True)
        plt.show()
    if Z == 3:
        # Plot 3 lines
        plt.plot(X, Y, Xfit, Yfit, 'r--', X2, Y2, 'g-*', linewidth=3)
        plt.title(TITLE)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        plt.legend([LEGEND01, LEGEND02, LEGEND03], loc='lower right')
        plt.grid(True)
        plt.show() 
        
def DUCTILITY_DAMAGE_INDEX_FUN(disp_PUSH, reaction_PUSH, SLOPE_NODE, disp_SEI, PLOT):
        import BILINEAR_CURVE as S07
        import numpy as np
        
        XX = np.abs(disp_PUSH); YY = np.abs(reaction_PUSH); # ABSOLUTE VALUE
        #SLOPE_NODE = 10
        DATA = S07.BILNEAR_CURVE(XX, YY, SLOPE_NODE)
        X, Y, Elastic_ST, Plastic_ST, Tangent_ST, Ductility_Rito, Over_Strength_Factor = DATA

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
        print(f'Structural Behavior Coefficient (R): {R:.4f}')
        Dd = np.max(np.abs(disp_SEI))
        DIx = 100*(Dd - X[1]) /(X[2] - X[1])
        print(f'Structural Ductility Damage Index:   {DIx:.4f} (%)')
        
        if PLOT == True:
            XLABEL = 'Displacement [m]'
            YLABEL = 'Base Reaction [N]'
            LEGEND01 = 'Curve'
            LEGEND02 = 'Bilinear Fitted'
            LEGEND03 = 'Undefined'
            TITLE = f'BaseShear-Displacement Analysis - Ductility Ratio: {X[2]/X[1]:.4f} - Over Strength Factor: {Y[2]/Y[1]:.4f}\n Ductility Damage Index: {DIx:.2f} [%] - Structural Behavior Coefficient (R): {R:.4f}'
            COLOR = 'black'
            PLOT_2D(np.abs(disp_PUSH), np.abs(reaction_PUSH), X, Y, X, Y, XLABEL, YLABEL, TITLE, LEGEND01, LEGEND02, LEGEND03, COLOR='black', Z=2) 
        #print(f'\t\t Ductility Ratio: {Y[2]/Y[1]:.4f}')
        
        return DIx, Omega_0, mu, R_mu, R