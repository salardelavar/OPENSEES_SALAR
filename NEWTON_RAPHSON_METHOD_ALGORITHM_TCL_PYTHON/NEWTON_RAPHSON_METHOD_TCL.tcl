# ############################################################################
#     >> IN THE NAME OF ALLAH, THE MOST GRACIOUS, THE MOST MERCIFUL <<       #
#                       NEWTON-RAPHSON ALGORITHM METHOD                      #
#       THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)        #
# ############################################################################

set X 100.0          ;# Initial guess
set ESP 1e-3         ;# Finite difference derivative convergence tolerance
set TOLERANCE 1e-6   ;# Convergence tolerance
set RESIDUAL 100.0   ;# Convergence residual (initialised large)
set IT 0             ;# Initial iteration count
set ITMAX 100000     ;# Max. iterations
set DEMAND 60.0      ;# Target Value

# Analysis durations:
set starttime [clock microseconds]

# ---------------------------------------------------------------------------
# FIND THE OPTIMUM VALUE (NEWTON-RAPHSON SOLVER FOR OPTIMAL X)
# ---------------------------------------------------------------------------
proc FUNCTION {x} {
    return [expr {$x**3 - 4.0 * $x**2 + $x - 250.0}]
}

while {$RESIDUAL > $TOLERANCE} {
    # X -------------------
    set answer [FUNCTION $X]
    set SUPPLY $answer
    set F [expr {$SUPPLY - $DEMAND}]
    puts "F:    $F"

    # XMIN -------------------
    set Xmin [expr {$X - $ESP}]
    puts "Xmin:    $Xmin"
    set answer [FUNCTION $Xmin]
    set SUPPLYmin $answer
    set Fmin [expr {$SUPPLYmin - $DEMAND}]
    puts "Fmin: $Fmin"

    # XMAX -------------------
    set Xmax [expr {$X + $ESP}]
    puts "Xmax:    $Xmax"
    set answer [FUNCTION $Xmax]
    set SUPPLYmax $answer
    set Fmax [expr {$SUPPLYmax - $DEMAND}]
    puts "Fmax: $Fmax"

    # DF -------------------
    set DF [expr {($Fmax - $Fmin) / (2.0 * $ESP)}]   ;# Finite difference derivative of F
    puts "DF:   $DF"

    # DX -------------------
    set DX [expr {$F / $DF}]                         ;# Newton step
    puts "DX:   $DX"

    # RESIDUAL -------------------
    set RESIDUAL [expr {abs($DX)}]
    puts "IT: [expr {$IT + 1}] - RESIDUAL: $RESIDUAL - X: $X\n"

    # Update X and iteration counter
    set X [expr {$X - $DX}]
    incr IT

    # CONTROLLING -------------------
    if {$IT == $ITMAX} {
        puts "\t\t Iteration reached to Max. Iteration"
        puts "\t\t Change ESP and TOLERANCE for better Convergence"
        set X [expr {-$DX}]
        break
    }
    if {$RESIDUAL < $TOLERANCE} {
        puts [format "\t\t Optimum X:  %.4f" $X]
        puts "\t\t Iteration Counts:                $IT"
        puts [format "\t\t Convergence Residual:            %.10e" $RESIDUAL]
    }
}

set endtime [clock microseconds]
set totaltime [expr {($endtime - $starttime) / 1e6}]
puts [format "\nTotal time (s): %.4f" $totaltime]