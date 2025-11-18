def SECTION_ANALYSIS():
    import numpy as np
    import matplotlib.pyplot as plt
    
    # Define rectangles: [x_centroid, y_centroid, width (b), height (h)]
    # Replace these with actual dimensions based on your section
    
    # x_centroid, y_centroid, width, height
    rects = np.array([
        [-55, 250, 10, 580],   # Right Web Plate
        [250, 520, 500, 20],   # Top Plate
        [250, 500, 600, 20],   # Top flange
        [450, 250, 20, 480],   # Right web
        [50, 250, 20, 480],    # Left web
        [250, 0, 600, 20],     # Bottom flange
        [250, -20, 500, 20],   # Bottom Plate
        [555, 250, 10, 580],   # Right Web Plate
    ])
    
    # Compute area of each rectangle
    areas = rects[:, 2] * rects[:, 3]
    
    # Composite centroid
    x_c = np.sum(areas * rects[:, 0]) / np.sum(areas)
    y_c = np.sum(areas * rects[:, 1]) / np.sum(areas)
    
    # Moment of inertia about centroidal axes
    Ix_total = 0
    Iy_total = 0
    A_total = 0
    
    for i in range(len(rects)):
        x, y, b, h = rects[i]
        A = b * h
        Ix_c = (b * h**3) / 12
        Iy_c = (h * b**3) / 12
        A_total += A 
        Ix_total += Ix_c + A * (y - y_c)**2
        Iy_total += Iy_c + A * (x - x_c)**2
    
    print(f"Composite Centroid: x̄ = {x_c:.2f}, ȳ = {y_c:.2f}")
    print(f"Moment of Inertia about x-axis (Ix): {Ix_total:.2f}")
    print(f"Moment of Inertia about y-axis (Iy): {Iy_total:.2f}")
    
    # Optional: Plot the section
    fig, ax = plt.subplots()
    for x, y, b, h in rects:
        ax.add_patch(plt.Rectangle((x - b/2, y - h/2), b, h, edgecolor='black', facecolor='lightgray'))
    ax.plot(x_c, y_c, 'ro', label='Centroid')
    ax.set_aspect('equal')
    ax.set_title('Steel Section')
    ax.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()
    
    return x_c, y_c, A_total, Ix_total, Iy_total
