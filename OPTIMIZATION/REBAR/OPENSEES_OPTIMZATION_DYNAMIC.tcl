#  #########################################################################
#  #                           IN THE NAME OF ALLAH                        #
#  #      STEEL REBAR DIAMETER AND CONCRETE CONFINED SECTION COLUMN        #
#  #            OPTIMIZATION ANALYSIS WITH NONLINEAR DYNAMIC               #
#  #                         NEWTON-RAPHSON METHOD                         #
#  #-----------------------------------------------------------------------#
#  #        THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)     #
#  #                   EMAIL: salar.d.ghashghaei@gmail.com                 #
#  #########################################################################

proc DYNAMIC_ANALYSIS {HCol BCol LCol coverCol Weight numBarsCol BD fc} {
    wipe
    model Basic -ndm 2 -ndf 3
    set PCol $Weight
    set g 9810
    set Mass [expr $PCol / $g]

    node 1 0.0 0.0
    node 2 0.0 $LCol
    fix 1 1 1 1
    mass 2 $Mass 1e-9 0.0

    set barAreaCol [expr (3.1415 * pow($BD, 2)) / 4]

    set ColSecTag 1
    set IDconcCore 1
    set IDconcCover 2
    set IDreinf 3

    set Ec [expr 4700 * sqrt(-$fc)]
    set Kfc 1.3
    set fc1C [expr $Kfc * $fc]
    set eps1C [expr 2 * $fc1C / $Ec]
    set fc2C [expr 0.2 * $fc1C]
    set eps2C [expr 5 * $eps1C]
    set fc1U $fc
    set eps1U -0.0025
    set fc2U [expr 0.2 * $fc1U]
    set eps2U -0.012
    set Lambda 0.1
    set ftC [expr -0.55 * $fc1C]
    set ftU [expr -0.55 * $fc1U]
    set Ets [expr $ftU / 0.002]

    set Fy 4000
    set Cy 0.02
    set Es [expr $Fy / $Cy]
    set Bs 0.01
    set R0 18.0
    set cR1 0.925
    set cR2 0.15

    uniaxialMaterial Concrete02 $IDconcCore $fc1C $eps1C $fc2C $eps2C $Lambda $ftC $Ets
    uniaxialMaterial Concrete02 $IDconcCover $fc1U $eps1U $fc2U $eps2U $Lambda $ftU $Ets
    uniaxialMaterial Steel02 $IDreinf $Fy $Es $Bs $R0 $cR1 $cR2

    set coverY [expr $HCol / 2.0]
    set coverZ [expr $BCol / 2.0]
    set coreY [expr $coverY - $coverCol]
    set coreZ [expr $coverZ - $coverCol]
    set coreY02 [expr $coreY - 50]
    set coreZ02 $coreZ

    set nfCoreY 16
    set nfCoreZ 4
    set nfCoverY 16
    set nfCoverZ 4

    section Fiber $ColSecTag {
        patch quad $IDconcCore $nfCoreZ $nfCoreY -$coreY $coreZ -$coreY -$coreZ $coreY -$coreZ $coreY $coreZ
        patch quad $IDconcCover $nfCoverZ $nfCoverY -$coverY $coverZ -$coreY $coreZ $coreY $coreZ $coverY $coverZ
        patch quad $IDconcCover $nfCoverZ $nfCoverY -$coreY -$coreZ -$coverY -$coverZ $coverY -$coverZ $coreY -$coreZ
        patch quad $IDconcCover $nfCoverZ $nfCoverY -$coverY $coverZ -$coverY -$coverZ -$coreY -$coreZ -$coreY $coreZ
        patch quad $IDconcCover $nfCoverZ $nfCoverY $coreY $coreZ $coreY -$coreZ $coverY -$coverZ $coverY $coverZ
        layer straight $IDreinf $numBarsCol $barAreaCol $coreY $coreZ $coreY -$coreZ
        layer straight $IDreinf 2 $barAreaCol $coreY02 $coreZ02 $coreY02 -$coreZ02
        layer straight $IDreinf 2 $barAreaCol -$coreY02 $coreZ02 -$coreY02 -$coreZ02
        layer straight $IDreinf $numBarsCol $barAreaCol -$coreY $coreZ -$coreY -$coreZ
    }

    set ColTransfTag 1
    geomTransf Linear $ColTransfTag
    set numIntgrPts 5
    set eleTag 1
    element nonlinearBeamColumn $eleTag 1 2 $numIntgrPts $ColSecTag $ColTransfTag

    recorder EnvelopeNode -file "MD.txt" -time -node 2 -dof 1 disp
    recorder EnvelopeNode -file "MV.txt" -time -node 2 -dof 1 disp
    recorder EnvelopeNode -file "MA.txt" -time -node 2 -dof 1 disp
    recorder Node -file "DTH.txt" -time -node 2 -dof 1 2 3 disp
    recorder Node -file "VTH.txt" -time -node 2 -dof 1 2 3 vel
    recorder Node -file "ATH.txt" -time -node 2 -dof 1 2 3 accel
    recorder Node -file "BTH.txt" -time -node 1 -dof 1 2 3 reaction
    recorder Element -file "FCol.txt" -time -ele 1 globalForce
    recorder Element -file "ForceColSec1.txt" -time -ele 1 section 1 force
    recorder Element -file "DCol.txt" -time -ele 1 deformations

    timeSeries Linear 1
    pattern Plain 1 1 {
        load 2 0.0 -$PCol 0.0
    }

    set Tol 1e-8
    set NstepGravity 10
    set DGravity [expr 1.0 / $NstepGravity]
    integrator LoadControl $DGravity
    numberer Plain
    system BandGeneral
    constraints Plain
    test NormDispIncr $Tol 6
    algorithm Newton
    analysis Static
    analyze $NstepGravity

    loadConst -time 0.0

    set GMdirection 1
    set GMfile "BM68elc.acc"
    set GMfact 1.0

    set Lambda01 [eigen -fullGenLapack 1]
    set Lambda02 [eigen -genBandArpack 1]
    set Omega [expr pow([expr min($Lambda01, $Lambda02)], 0.5)]
    set betaKcomm [expr 2 * (0.02 / $Omega)]

    set xDamp 0.02
    set alphaM 0.0
    set betaKcurr 0.0
    set betaKinit 0.0

    rayleigh $alphaM $betaKcurr $betaKinit $betaKcomm

    set IDloadTag 400
    set dt 0.01
    set GMfatt 1.0
    set maxNumIter 10
    timeSeries Path 2 -dt $dt -filePath $GMfile -factor $GMfact
    pattern UniformExcitation $IDloadTag $GMdirection -accel 2

    wipeAnalysis
    constraints Transformation
    numberer Plain
    system BandGeneral
    test EnergyIncr $Tol $maxNumIter
    algorithm ModifiedNewton

    set NewmarkGamma 0.5
    set NewmarkBeta 0.25
    integrator Newmark $NewmarkGamma $NewmarkBeta
    analysis Transient

    set DtAnalysis 0.01
    set TmaxAnalysis 10.0
    set Nsteps [expr int($TmaxAnalysis / $DtAnalysis)]

    set ok [analyze $Nsteps $DtAnalysis]

    set tCurrent [getTime]

    set test {1 NormDispIncr 2 RelativeEnergyIncr 4 RelativeNormUnbalance 5 RelativeNormDispIncr 6 NormUnbalance}
    set algorithm {1 KrylovNewton 2 SecantNewton 4 RaphsonNewton 5 PeriodicNewton 6 BFGS 7 Broyden 8 NewtonLineSearch}

    foreach i $test {
        foreach j $algorithm {
            if {$ok != 0} {
                if {$j < 4} {
                    algorithm [lindex $algorithm $j] -initial
                } else {
                    algorithm [lindex $algorithm $j]
                }
                test [lindex $test $i] $Tol 1000
                set ok [analyze $Nsteps $DtAnalysis]
                puts "[lindex $test $i] [lindex $algorithm $j] $ok"
                if {$ok == 0} {
                    break
                }
            } else {
                continue
            }
        }
    }

    puts "Dyn Done."
    wipe
}

# -------------------------

proc MAXABS_FUN {DATA_FILE} {
    # Read and process displacement data
    set NameFiles $DATA_FILE
    set filename "${NameFiles}.txt"
    
    # Load the data from the file
    set fileId [open $filename r]
    set data [read $fileId]
    close $fileId
    
    # Split the data into lines and then into fields
    set lines [split $data "\n"]
    set D {}
    foreach line $lines {
        if {[string length $line] > 0} {
            lappend D [split $line]
        }
    }
    
    # Extract the second column and compute the max absolute value
    set maxAbs 0
    foreach row $D {
        set value [lindex $row 1]
        set absValue [expr {abs($value)}]
        if {$absValue > $maxAbs} {
            set maxAbs $absValue
        }
    }
    
    return $maxAbs
}
# -------------------------

# --------------------------------------------------------------
#                COLUMN REBAR DIAMETER OPTIMIZATION
# --------------------------------------------------------------

# define section geometry
set LCol 30000.0 ;# [mm] column length
set HCol 300 ;# [mm] Column Depth
set BCol 300 ;# [mm] Column Width
set coverCol 50.0 ;# [mm] Column cover to reinforcing steel NA.
set Weight 1000000.0 ;# [N] superstructure weight
set numBarsCol 4 ;# number of longitudinal-reinforcement bars in column. (symmetric top & bot)
set BD 20 ;# [mm] Rebar Diamater
set fc -25.0 ;# [N/mm^2] Concrete Compressive Strength (+Tension, -Compression)
set X $BD ;# Initial Guess

set ESP 1e-3 ;# Finite difference derivative Convergence Tolerance
set TOLERANCE 1e-6 ;# Convergence Tolerance
set RESIDUAL 100 ;# Convergence Residual 
set IT 0 ;# Initial Iteration
set ITMAX 100000 ;# Max. Iteration
set TARGET_PGD 20 ;# [mm] Target Demand Max. Abs. Displacement Time History

set DATA_FILE "C:/OPENSEESPY_SALAR/OPENSEESPY_DATA/MD" ;# MAX DISPLACEMENT
# monitor cpu time
set starttime [clock seconds]

# FIND THE OPTIMUM VALUE 
while {$RESIDUAL > $TOLERANCE} {
    DYNAMIC_ANALYSIS $HCol $BCol $LCol $coverCol $Weight $numBarsCol $X $fc
    after 10000 ;# Sleep for 10 seconds
    set F [expr [MAXABS_FUN $DATA_FILE] - $TARGET_PGD]
    puts "Current Max. Abs. displacement: [MAXABS_FUN $DATA_FILE]"
    
    # Evaluate at Xmain and Fmin
    set Xmin [expr $X - $ESP]
    DYNAMIC_ANALYSIS $HCol $BCol $LCol $coverCol $Weight $numBarsCol $Xmin $fc
    after 10000 ;# Sleep for 10 seconds
    set Fmin [expr [MAXABS_FUN $DATA_FILE] - $TARGET_PGD]
    
    # Evaluate at Xmax and Fmax
    set Xmax [expr $X + $ESP]
    DYNAMIC_ANALYSIS $HCol $BCol $LCol $coverCol $Weight $numBarsCol $Xmax $fc
    after 10000 ;# Sleep for 10 seconds
    set Fmax [expr [MAXABS_FUN $DATA_FILE] - $TARGET_PGD]
    
    set DF [expr ($Fmax - $Fmin) / (2 * $ESP)] ;# Calculate the Finite difference derivative of F
    set DX [expr $F / $DF] ;# Calculate dx
    set RESIDUAL [expr abs($DX)] ;# Calculate residual
    puts "RESIDUAL: $RESIDUAL"
    set X [expr $X - $DX] ;# update X
    incr IT ;# update iteration
    
    if {$IT == $ITMAX} {
        puts "\t\t Iteration reached to Max. Iteration"
        puts "\t\t Change ESP and TOLERANCE for better Convergence"
        set X -1
        break
    }
    
    if {$RESIDUAL < $TOLERANCE} {
        puts "\t\t Optimum Rebar Diamater: [format %.4f $X]"
        puts "\t\t Iteration Counts: $IT"
        puts "\t\t Convergence Residual [format %.10e $RESIDUAL]"
    }
}

set totaltime [expr [clock seconds] - $starttime]
puts "\nTotal time (s): [format %.4f $totaltime] \n\n"