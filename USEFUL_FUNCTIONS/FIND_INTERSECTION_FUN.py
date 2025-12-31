# Find intersection points between Response Spectrum and Pushover curve
def FIND_INTERSECTION(x1, y1, x2, y2, DATA01_LABEL, DATA02_LABEL, X_LABEL, Y_LABEL):
    """
    Find intersection points between two curves
    x1, y1: points of first curve (Response Spectrum)
    x2, y2: points of second curve (Pushover)
    """
    import numpy as np
    intersections = []
    
    # Convert to arrays for computation
    x1 = np.array(x1)
    y1 = np.array(y1)
    x2 = np.array(x2)
    y2 = np.array(y2)
    
    # Interpolate curves to continuous functions
    from scipy import interpolate
    
    # Interpolate first curve (Response Spectrum)
    f1 = interpolate.interp1d(x1, y1, kind='cubic', bounds_error=False, fill_value='extrapolate')
    
    # Interpolate second curve (Pushover) - may need sorting
    # Sort Pushover points by T
    sorted_indices = np.argsort(x2)
    x2_sorted = x2[sorted_indices]
    y2_sorted = y2[sorted_indices]
    f2 = interpolate.interp1d(x2_sorted, y2_sorted, kind='linear', bounds_error=False, fill_value='extrapolate')
    
    # Find intersection points in common range
    x_min = max(np.min(x1), np.min(x2))
    x_max = min(np.max(x1), np.max(x2))
    
    if x_min < x_max:
        # Create function to find roots
        def diff_func(x):
            return f1(x) - f2(x)
        
        # Sample to find sign changes
        x_samples = np.linspace(x_min, x_max, 1000)
        diff_samples = diff_func(x_samples)
        
        # Find intervals where function changes sign
        sign_changes = np.where(np.diff(np.sign(diff_samples)))[0]
        
        for idx in sign_changes:
            # Use bisection method to find precise intersection
            a = x_samples[idx]
            b = x_samples[idx + 1]
            
            try:
                # Use brentq to find root
                from scipy.optimize import brentq
                intersection_x = brentq(diff_func, a, b)
                intersection_y = f1(intersection_x)
                
                # Ensure point is in reasonable range
                if (x_min <= intersection_x <= x_max and 
                    intersection_y >= 0 and 
                    not np.isnan(intersection_x) and 
                    not np.isnan(intersection_y)):
                    intersections.append((intersection_x, intersection_y))
            except:
                continue
    
    # Create a separate plot for better visualization of intersection points
    if intersections:
        print(f"\n{'='*60}")
        print(f"Intersection Points between {DATA01_LABEL} and {DATA02_LABEL}:")
        print('='*60)
        
        intersection_x = []
        intersection_y = []
        
        for i, (x, y) in enumerate(intersections):
            print(f"Intersection Point {i+1}:")
            print(f"  - {X_LABEL}: {x:.6f}")
            print(f"  - {Y_LABEL}: {y:.6f}")
            print('-'*40)
            
            intersection_x.append(x)
            intersection_y.append(y)  
            
        import matplotlib.pyplot as plt
        plt.figure(figsize=(12, 10))
        # Plot curves
        plt.plot(x1, y1, 'b-', linewidth=3, label=DATA01_LABEL)
        plt.plot(x2, y2, 'r--', linewidth=3, label=DATA02_LABEL, alpha=0.7)
        
        # Intersection points
        plt.plot(intersection_x, intersection_y, 'm*', markersize=20, 
                 markeredgecolor='black', markeredgewidth=2, 
                 label=f'Intersection Points: X: {intersection_x} - Y: {intersection_y}')
        
        # Guide lines from intersection points
        for x, y in intersections:
            plt.axvline(x=x, color='gray', linestyle=':', alpha=0.5)
            plt.axhline(y=y, color='gray', linestyle=':', alpha=0.5)
            plt.plot([x, x], [0, y], 'g:', alpha=0.5)
            plt.plot([0, x], [y, y], 'g:', alpha=0.5)
        
        plt.xlabel(X_LABEL, fontsize=12)
        plt.ylabel(Y_LABEL, fontsize=12)
        plt.title(f'Intersection Points between {DATA01_LABEL} and {DATA02_LABEL}', fontsize=14)
        plt.legend(loc='best')
        plt.grid(True, alpha=0.3)
        #plt.xlim(0, max(intersection_x) * 1.3)
        #plt.ylim(0, max(intersection_y) * 1.3)
        plt.tight_layout()
        plt.show()
    else:
        print("No intersection points found between {DATA01_LABEL} and {DATA01_LABEL}")
         
        
    return intersections