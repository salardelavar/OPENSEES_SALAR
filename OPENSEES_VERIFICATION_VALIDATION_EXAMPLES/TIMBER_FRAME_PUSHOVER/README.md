# NONLINEAR PUSHOVER ANALYSIS OF A ONE-BAY, ONE-STORY TIMBER FRAME WITH SEMI-RIGID AND RIGID CONNECTIONS IN OPENSEES

# تحلیل پوش‌آور غیرخطی یک قاب چوبی یک‌دهانه و یک‌طبقه با اتصالات نیمه‌صلب و صلب در نرم‌افزار اوپنسیس

This Python script develops a nonlinear finite-element model of a one-bay, one-story (1B1S) timber frame with
 hysteretic beam-to-column and column-to-foundation connections using OpenSeesPy. The frame has a single 4,100
 mm bay and a story height of 2,740 mm, with elastic beam and column members connected through zero-length elements
 representing nonlinear rotational behavior. The connection response is defined using the `Pinching4` uniaxial material
 model, which captures pinching effects, stiffness and strength degradation, and energy-based deterioration under cyclic
 loading. Alternatively, the beam-to-column connections can be assigned rigid elastic behavior by changing `BCJType` to 
 `"Rigid"`. The model uses elastic material properties for the timber structural members and highly stiff elastic materials
 to represent shear constraints at the connections. After defining the nodal coordinates, elements, boundary conditions, and
 connection properties, the script establishes nodal recorders to capture roof displacement and base reactions.
 A displacement-controlled nonlinear static pushover analysis is then performed by incrementally increasing the horizontal
 displacement of the roof control node up to 260 mm, with a displacement increment of 1 mm per analysis step.
 The resulting roof displacement and accumulated base shear are extracted from the recorder files and used to generate the
 global pushover capacity curve. Finally, the `PLOT_2D_FRAME` function is called to visualize the undeformed and deformed
 frame configurations. Overall, this model provides a computational framework for investigating the nonlinear lateral response,
 hysteretic connection behavior, and global seismic capacity of a one-bay, one-story timber frame with semi-rigid or rigid connection assumptions.
 
![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/OPENSEES_VERIFICATION_VALIDATION_EXAMPLES/TIMBER_FRAME_PUSHOVER/COVER.png)
