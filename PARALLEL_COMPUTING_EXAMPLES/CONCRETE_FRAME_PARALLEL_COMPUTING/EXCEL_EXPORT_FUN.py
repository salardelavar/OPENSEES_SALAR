#%% EXPORT DATA TO EXCEL FOR PUSHOVER ANALYSIS
def EXPORT_DATA_STATIC(FORCE_S, FORCE_A, MOMENT, DISP_X, DISP_Y, ROT, KA, KS, KI, STEP):
    import pandas as pd
    # Create a dictionary of the data
    dfP = pd.DataFrame({
        'STEP': STEP,
        'FORCE_S': FORCE_S,
        'FORCE_A': FORCE_A,
        'MOMENT': MOMENT,
        'DISP_X': DISP_X,
        'DISP_Y': DISP_Y,
        'ROT': ROT,
        'KA': KA,
        'KS': KS,
        'KI': KI
        })
        
    # Write to Excel
    dfP.to_excel('PUSHOVER_ANALYSIS_RESULTS.xlsx', index=False)  # index=False avoids writing row numbers

#%% EXPORT DATA TO EXCEL FOR DYNAMIC ANALYSIS
def EXPORT_DATA_DYNAMIC(time, DISP_X, DISP_Y, velocity_X, velocity_Y, acceleration_X, acceleration_Y, FORCE_S, FORCE_A, MOMENT, ROT, delta, KA, KS, KI, PERIOD_01, PERIOD_02):
    import pandas as pd
    # Create a DataFrame from time series data
    dfD = pd.DataFrame({
        'time': time,
        'FORCE_S': FORCE_S,
        'FORCE_A': FORCE_A,
        'MOMENT': MOMENT,
        'DISP_X': DISP_X,
        'DISP_Y': DISP_Y,
        'ROT': ROT,
        'KA': KA,
        'KS': KS,
        'KI': KI,
        'velocity_X': velocity_X,
        'velocity_Y': velocity_Y,
        'acceleration_X': acceleration_X,
        'acceleration_Y': acceleration_Y
        })
            
    # Write time series data to an Excel file
    with pd.ExcelWriter('DYNAMIC_ANALYSIS_RESULTS.xlsx', engine='openpyxl') as writer:
        dfD.to_excel(writer, sheet_name='Time Series Data', index=False)
            
        # Save scalar values in a new sheet
        scalar_df = pd.DataFrame({
            'PERIOD_01': [PERIOD_01],
            'PERIOD_02': [PERIOD_02],
            'delta': [delta]
            })
        scalar_df.to_excel(writer, sheet_name='Summary', index=False)
