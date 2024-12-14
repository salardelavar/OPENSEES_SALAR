#          #####################################################################################
#          #                                  IN THE NAME OF ALLAH                             #
#          #          FREE VIBRATION ANAYSIS OF INELASTIC SPRING WITH PULSE LOADING            #
#          #-----------------------------------------------------------------------------------#
#          #              THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)           #
#          #                       EMAIL: salar.d.ghashghaei@gmail.com                         #
#          #####################################################################################

# Model parameters
set k 1250.0  ;# Stiffness of the spring
set E 210000.0  ;# Young's modulus in kPa
set A 0.25  ;# Cross-sectional area
set L 1000.0  ;# Length of the spring element
set k [expr $E * $A / $L]  ;# Truss Axial Stiffness
set m 50.0  ;# Mass
set u0 0.001  ;# Initial displacement
set damping_ratio 0.02  ;# Damping ratio

# Analysis parameters
set duration 50.0 ;# [s]
set dt 0.01       ;# Time step

# ReinforcingSteel material properties
set fy 400.0  ;# Yield strength
set Es 210000.0  ;# Elastic modulus
set fu 600.0  ;# Ultimate strength
set Esh 20000.0 ;# Hardening modulus
set esh 0.01  ;# Strain at start of hardening
set esu 0.1   ;# Ultimate strain

# Pulse loading parameters
set START_TIME 2.0
set END_TIME 6.0
set PERIOD 1.2
set LOAD_FACTOR 2.0
set PULSEWIDTH 2.5 ;# Fraction of the period

# Procedure to perform the analysis
proc perform_analysis {damping output_file} {
    wipe
    model BasicBuilder -ndm 2 -ndf 3

    # Define nodes
    node 1 0.0 0.0
    node 2 $L 0.0

    # Fixities
    fix 1 1 1 1
    fix 2 0 1 1

    # Mass
    mass 2 $m 0.0 0.0

    # Material
    uniaxialMaterial ReinforcingSteel 1 $fy $Es $fu $Esh $esh $esu

    # Truss element
    element Truss 1 1 2 $A 1

    # Apply initial static displacement
    timeSeries Linear 1
    pattern Plain 1 1 {
        load 2 1.0 0.0 0.0
    }
    constraints Transformation
    numberer RCM
    system BandGeneral
    algorithm Linear
    test NormDispIncr 1.0e-8 10
    integrator DisplacementControl 2 1 $u0
    analysis Static
    analyze 1
    loadConst -time 0.0

    # Remove static load
    wipeAnalysis
    remove loadPattern 1
    system UmfPack

    # Define triangle loading 
    #timeSeries Triangle 2 $START_TIME $END_TIME #PERIOD -factor $LOAD_FACTOR
    # Define dynamic pulse loading
    timeSeries Rectangular 2 0.0 2.0 -factor $LOAD_FACTOR

    pattern Plain 2 2 {
        load 2 1.0 0.0 0.0
    }

    # Dynamic analysis setup
    constraints Transformation
    numberer RCM
    system UmfPack
    test NormDispIncr 1.0e-8 10
    integrator Newmark 0.5 0.25
    algorithm Newton
    analysis Transient

    # Apply Rayleigh damping if required
    if {$damping} {
        set omega1 [expr sqrt($k / $m)]
        set omega2 [expr 2.0 * $omega1]
        set a0 [expr $damping_ratio * (2.0 * $omega1 * $omega2) / ($omega1 + $omega2)]
        set a1 [expr $damping_ratio * 2.0 / ($omega1 + $omega2)]
        rayleigh $a0 $a1 0.0 0.0
    }

    # Open the output file
    set fid [open $output_file "w"]
    puts $fid "Time\tDisplacement\tVelocity\tAcceleration\tSpringForce"

    # Perform transient analysis and store results
    set current_time 0.0
    while {$current_time < $duration} {
        set stable [analyze 1 $dt]
        if {$stable != 0} {
            break
        }
        set current_time [getTime]
        set disp [nodeDisp 2 1]
        set vel [nodeVel 2 1]
        set accel [nodeAccel 2 1]
        set spring_force [expr -1.0 * $k * $disp]

        # Write results to the file
        puts $fid "$current_time\t$disp\t$vel\t$accel\t$spring_force"
    }

    # Close the output file
    close $fid

    puts "Results written to $output_file"
}

# Run analyses
perform_analysis 0 "undamped_results.txt"
perform_analysis 1 "damped_results.txt"

puts "Analysis Complete"
