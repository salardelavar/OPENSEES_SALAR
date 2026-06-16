# ############################################################################
#     >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<       #
#                       NEWTON-RAPHSON ALGORITHM METHOD                      #
#     THIS TCL SCRIPT TRANSLATED FROM PYTHON BY SALAR DELAVAR GHASHGHAEI     #
# ############################################################################

# Initial Guess for X1 and X2
set X01 100.0
set X02 100.0

# Finite difference derivative Convergence Tolerance
set ESP 1e-3

# Convergence Tolerance
set TOLERANCE 1e-6

# Convergence Residual
set RESIDUAL 100.0

# Initial Iteration
set IT 0

# Max. Iteration
set ITMAX 100000

# Target Value
set DEMAND 0.0

# Analysis Duration start
set starttime [clock milliseconds]

# ---------------------------------------------------------------------------
# FIND THE OPTIMUM VALUE (NEWTON-RAPHSON SOLVER FOR OPTIMAL X1 and X2)
# ---------------------------------------------------------------------------
proc FUNCTION {X1 X2} {
    set A01 [expr {$X1**2 + $X2**2 - 25}]   ;# FUNCTION 01
    set A02 [expr {$X1*$X2 - 12}]            ;# FUNCTION 02
    return [list $A01 $A02]
}

while {$RESIDUAL > $TOLERANCE} {
    # X -------------------
    lassign [FUNCTION $X01 $X02] SUPPLY01 SUPPLY02
    set F01 [expr {$SUPPLY01 - $DEMAND}]
    set F02 [expr {$SUPPLY02 - $DEMAND}]
    puts "F:     $F01 $F02"
    
    # XMIN -------------------
    # Evaluate at Xmin and Fmin
    set Xmin01 [expr {$X01 - $ESP}]
    set Xmin02 [expr {$X02 - $ESP}]
    puts "Xmin:  $Xmin01 $Xmin02"
    lassign [FUNCTION $Xmin01 $Xmin02] SUPPLYmin01 SUPPLYmin02
    set Fmin01 [expr {$SUPPLYmin01 - $DEMAND}]
    set Fmin02 [expr {$SUPPLYmin02 - $DEMAND}]
    puts "Fmin:  $Fmin01 $Fmin02"
    
    # XMAX -------------------
    # Evaluate at Xmax and Fmax
    set Xmax01 [expr {$X01 + $ESP}]
    set Xmax02 [expr {$X02 + $ESP}]
    puts "Xmax:  $Xmax01 $Xmax02"
    lassign [FUNCTION $Xmax01 $Xmax02] SUPPLYmax01 SUPPLYmax02
    set Fmax01 [expr {$SUPPLYmax01 - $DEMAND}]
    set Fmax02 [expr {$SUPPLYmax02 - $DEMAND}]
    puts "Fmax:  $Fmax01 $Fmax02"
    
    # DF -------------------
    # Calculate the Finite difference derivative of F1
    set DF01 [expr {($Fmax01 - $Fmin01) / (2.0 * $ESP)}]
    # Calculate the Finite difference derivative of F2
    set DF02 [expr {($Fmax02 - $Fmin02) / (2.0 * $ESP)}]
    puts "DF:    $DF01 $DF02"
    
    # DX -------------------
    set DX01 [expr {$F01 / $DF01}]   ;# Calculate dx
    set DX02 [expr {$F02 / $DF02}]   ;# Calculate dx
    puts "DX:    $DX01 $DX02"
    
    # RESIDUAL -------------------
    # Calculate residual
    if {[expr {abs($DX01)}] > [expr {abs($DX02)}]} {
        set RESIDUAL [expr {abs($DX01)}]
    } else {
        set RESIDUAL [expr {abs($DX02)}]
    }
    
    incr IT
    puts "IT: $IT - RESIDUAL: $RESIDUAL - X01: $X01 - X02: $X02\n"
    
    # Update X1 and X2
    set X01 [expr {$X01 - $DX01}]
    set X02 [expr {$X02 - $DX02}]
    
    # CONTROLLING -------------------
    if {$IT == $ITMAX} {
        puts "\t\t Iteration reached to Max. Iteration"
        puts "\t\t Change ESP and TOLERANCE for better Convergence"
        break
    }
    
    if {$RESIDUAL < $TOLERANCE} {
        puts [format "\t\t Optimum X1:                      %.4f" $X01]
        puts [format "\t\t Optimum X2:                      %.4f" $X02]
        puts "\t\t Iteration Counts:                $IT"
        puts [format "\t\t Convergence Residual:            %.10e" $RESIDUAL]
    }
}

set endtime [clock milliseconds]
set totaltime [expr {($endtime - $starttime) / 1000.0}]
puts [format "\nTotal time (s): %.4f \n\n" $totaltime]