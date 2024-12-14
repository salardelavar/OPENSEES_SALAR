###########################################################################################################
#                                                 IN THE NAME OF ALLAH                                    #
#---------------------------------------------------------------------------------------------------------#
#                                   THIS PROGRAM IS WRITTEN BY SALAR DELAVAR GHASHGHAEI                   #
#                                          EMAIL: SALAR.D.GHASHGHAEI@GMAIL.COM                            #
###########################################################################################################
#
# Title:
# Dynamic Response Analysis of a Single-Degree-of-Freedom (SDOF) System Subjected to Explosion Loading Using OpenSees
#
# Target:
# This code simulates the dynamic response of a single-degree-of-freedom (SDOF) structural system subjected to explosion-induced loading.
# It achieves the following:
#
# 1. Define Explosion Loading:
#    - Implements the Friedlander equation to model pressure-time history from an explosive event.
#    - Plots the explosion pressure profile over time.
#
# 2. Structural Model:
#    - Models the SDOF system in OpenSeesPy with a linear elastic spring (stiffness `k`), mass `m`, and damping ratio.
#    - Applies the explosion-induced time-dependent loading to the mass node.
#
# 3. Dynamic Analysis:
#    - Simulates the time history response using the Newmark method for transient analysis.
#    - Tracks system responses, including displacement, velocity, acceleration, and base reactions.
#
# 4. Visualization:
#    - Generates plots for displacement, velocity, acceleration, and base reactions to evaluate the impact of the explosion loading on the structure.
#
# This simulation is useful for structural engineers studying the effects of blast loads on structures, aiding in the design and assessment of resilient systems.

# Define parameters for the Friedlander equation
set P0 1.0e5               ;# Peak pressure (Pa)
set t0 0.1                ;# Positive phase duration (s)
set A 1.3                 ;# Wave decay coefficient
set dt 0.001              ;# Time step (s)
set impact_duration 2.0   ;# Total duration of explosion impact (s)
set duration 10.0         ;# Total simulation duration (s)

# Define structure parameters
set k 1.0e6               ;# Stiffness of the structure (N/m)
set m 1000.0              ;# Mass of the structure (kg)
set damping_ratio 0.05    ;# Damping ratio

# Generate Friedlander explosion pressure time series
set time_series ""
set time 0
while {$time <= $impact_duration} {
    if {$time < 0} {
        set pressure 0
    } else {
        set pressure [expr $P0 * (1 - $time / $t0) * exp(-$A * $time / $t0)]
    }
    set time_series "$time_series $pressure"
    set time [expr $time + $dt]
}

# Initialize OpenSees model
wipe
model Basic -ndm 1 -ndf 1

# Define nodes
node 1 0.0    ;# Fixed base
node 2 0.0    ;# Mass node

# Define boundary conditions
fix 1 1

# Define mass
mass 2 $m

# Calculate natural frequency and damping coefficient
set wn [expr sqrt($k / $m)]
set CS [expr 2 * $wn * $m * $damping_ratio]

# Define materials
uniaxialMaterial Elastic 1 $k
uniaxialMaterial Elastic 2 0.0 $CS

# Define element
element zeroLength 1 1 2 -mat 1 2 -dir 1

# Define time series for explosion loading
timeSeries Path 1 -dt $dt -values $time_series
pattern Plain 1 1 {
    load 2 1.0
}

# Define analysis parameters
constraints Plain
numberer Plain
system BandGeneral
test NormDispIncr 1.0e-6 10
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient

# Perform dynamic analysis and write results to text files
set time_list ""
set disp_list ""
set vel_list ""
set accel_list ""
set base_reaction_list ""

set current_time 0.0
while {$current_time < $duration} {
    if {[analyze 1 $dt] != 0} {
        puts "Analysis failed at time $current_time"
        break
    }

    set current_time [getTime]
    lappend time_list $current_time
    lappend disp_list [nodeDisp 2 1]
    lappend vel_list [nodeVel 2 1]
    lappend accel_list [nodeAccel 2 1]
    lappend base_reaction_list [reaction 1 1]
}

# Write data to text files
set out_time [open "time.txt" "w"]
set out_disp [open "displacement.txt" "w"]
set out_vel [open "velocity.txt" "w"]
set out_accel [open "acceleration.txt" "w"]
set out_reaction [open "base_reaction.txt" "w"]

foreach t $time_list d $disp_list v $vel_list a $accel_list r $base_reaction_list {
    puts $out_time $t
    puts $out_disp $d
    puts $out_vel $v
    puts $out_accel $a
    puts $out_reaction $r
}

close $out_time
close $out_disp
close $out_vel
close $out_accel
close $out_reaction

puts "Analysis completed and data written to files."
