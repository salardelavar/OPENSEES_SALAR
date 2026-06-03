proc COMPOSITE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_TWO {secTag STEEL_TYPE fc Kfc CONCRETE_DENSITY} {
    # THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
    # ------------------------------------------------------------------
    # 1. Material properties (units: mm, N)
    # ------------------------------------------------------------------
    set fcU   [expr { -$fc }]
    set Ec    [expr { 4700.0 * sqrt($fcU) }]
    set ec0U  [expr { 2.0 * $fcU / $Ec }]
    set fcUU  [expr { 0.2 * $fcU }]
    set ecuU  [expr { 5.0 * $ec0U }]
    set LambdaU 0.1

    set fcC   [expr { $Kfc * $fcU }]
    set EcC   [expr { 4700.0 * sqrt($fcC) }]
    set ec0C  [expr { 2.0 * $fcC / $EcC }]
    set fcUC  [expr { 0.65 * $fcC }]
    set ecuC  [expr { 15.0 * $ec0C }]
    set LambdaC 0.1

    set ftC   [expr { 0.7 * sqrt($fcC) }]
    set ftU   [expr { 0.7 * sqrt($fcU) }]
    set EtsC  [expr { $ftC / abs($ec0C) }]
    set EtsU  [expr { $ftU / abs($ec0U) }]

    # Reinforcing steel (rebars)
    set fy    400.0
    set Es    200000.0
    set ey    [expr { $fy / $Es }]
    set fu    [expr { 1.1818 * $fy }]
    set esu   0.09
    set Esh   [expr { ($fu - $fy) / ($esu - $ey) }]
    set Bs    [expr { $Esh / $Es }]

    # Steel I‑section (plates and built‑up shape)
    set fyI   240.0
    set EsI   200000.0
    set eyI   [expr { $fyI / $EsI }]
    set fuI   [expr { 1.1818 * $fyI }]
    set esuI  0.25
    set EshI  [expr { ($fuI - $fyI) / ($esuI - $eyI) }]
    set BsI   [expr { $EshI / $EsI }]

    # Material tags
    set coreTag    [expr { $secTag + 100 }]
    set coverTag   [expr { $secTag + 200 }]
    set steelTag   [expr { $secTag + 300 }]
    set steelITag  [expr { $secTag + 400 }]

    # ------------------------------------------------------------------
    # 2. Define uniaxial materials
    # ------------------------------------------------------------------
    if {$STEEL_TYPE == "ELASTIC"} {
        uniaxialMaterial Elastic $steelTag   $Es
        uniaxialMaterial Elastic $steelITag  $EsI
    } elseif {$STEEL_TYPE == "INELASTIC"} {
        set pinchX  0.8
        set pinchY  0.5
        set damage1 0.0
        set damage2 0.0
        set beta    0.1
        uniaxialMaterial Hysteretic $steelTag \
            $fy $ey $fu $esu [expr {0.2*$fu}] [expr {1.1*$esu}] \
            [expr {-$fy}] [expr {-$ey}] [expr {-$fu}] [expr {-$esu}] \
            [expr {-0.2*$fu}] [expr {-1.1*$esu}] $pinchX $pinchY $damage1 $damage2 $beta
        uniaxialMaterial Hysteretic $steelITag \
            $fyI $eyI $fuI $esuI [expr {0.2*$fuI}] [expr {1.1*$esuI}] \
            [expr {-$fyI}] [expr {-$eyI}] [expr {-$fuI}] [expr {-$esuI}] \
            [expr {-0.2*$fuI}] [expr {-1.1*$esuI}] $pinchX $pinchY $damage1 $damage2 $beta
    } else {
        error "STEEL_TYPE must be ELASTIC or INELASTIC"
    }

    uniaxialMaterial Concrete02 $coreTag  $fcC  $ec0C  $fcUC  $ecuC  $LambdaC  $ftC  $EtsC
    uniaxialMaterial Concrete02 $coverTag $fcU  $ec0U  $fcUU  $ecuU  $LambdaU  $ftU  $EtsU   ;# not used, but defined

    # ------------------------------------------------------------------
    # 3. Create fiber section
    # ------------------------------------------------------------------
    section Fiber $secTag

    # Section geometry parameters
    set PT   10.0
    set CD  500.0
    set CW  700.0
    set N     7
    set stepY [expr { $CD / $N }]
    set stepX [expr { $CW / $N }]
    set NFY  10
    set NFX  10

    # Variables for area and height computation
    set totalArea 0.0
    set minBottom 1e20
    set maxTop   -1e20

    # Helper to add a rectangular patch and update area/bounds
    proc addPatch {matTag numY numX xLeft yBot xRight yTop} {
        upvar totalArea totalArea minBottom minBottom maxTop maxTop
        patch rect $matTag $numY $numX $xLeft $yBot $xRight $yTop
        set width  [expr { $xRight - $xLeft }]
        set depth  [expr { $yTop   - $yBot  }]
        set area   [expr { $width * $depth }]
        set totalArea [expr { $totalArea + $area }]
        if {$yBot < $minBottom} { set minBottom $yBot }
        if {$yTop > $maxTop}    { set maxTop    $yTop }
    }

    # 3a. Concrete grid patches (all confined concrete)
    for {set i 0} {$i < $N} {incr i} {
        for {set j 0} {$j < $N} {incr j} {
            set centerY [expr { ($i + 0.5) * $stepY }]
            set centerX [expr { ($j + 0.5) * $stepX }]
            set yBot    [expr { $centerY - $stepY/2.0 }]
            set yTop    [expr { $centerY + $stepY/2.0 }]
            set xLeft   [expr { $centerX - $stepX/2.0 }]
            set xRight  [expr { $centerX + $stepX/2.0 }]
            addPatch $coreTag $NFY $NFX $xLeft $yBot $xRight $yTop
        }
    }

    # 3b. Four outer steel plates (using steelITag)
    # Top plate
    set topYcenter [expr { $CD + $PT/2.0 }]
    addPatch $steelITag $NFY $NFX \
        [expr { -$PT }] [expr { $topYcenter - $PT/2.0 }] \
        [expr { $CW + $PT }] [expr { $topYcenter + $PT/2.0 }]
    # Bottom plate
    set botYcenter [expr { -$PT/2.0 }]
    addPatch $steelITag $NFY $NFX \
        [expr { -$PT }] [expr { $botYcenter - $PT/2.0 }] \
        [expr { $CW + $PT }] [expr { $botYcenter + $PT/2.0 }]
    # Right plate
    set rightXcenter [expr { $CW + $PT/2.0 }]
    addPatch $steelITag $NFY $NFX \
        [expr { $rightXcenter - $PT/2.0 }] 0.0 \
        [expr { $rightXcenter + $PT/2.0 }] $CD
    # Left plate
    set leftXcenter [expr { -$PT/2.0 }]
    addPatch $steelITag $NFY $NFX \
        [expr { $leftXcenter - $PT/2.0 }] 0.0 \
        [expr { $leftXcenter + $PT/2.0 }] $CD

    # 3c. Built‑up steel shape (cross‑section) – all using steelITag
    # Horizontal bar (10x150) at y = 125
    addPatch $steelITag $NFY $NFX \
        [expr {0.5*$CW - 75.0}] [expr {0.5*$CD - 125.0 - 5.0}] \
        [expr {0.5*$CW + 75.0}] [expr {0.5*$CD - 125.0 + 5.0}]
    # Vertical bar (240x10) at centre
    addPatch $steelITag $NFY $NFX \
        [expr {0.5*$CW - 5.0}] [expr {0.5*$CD - 120.0}] \
        [expr {0.5*$CW + 5.0}] [expr {0.5*$CD + 120.0}]
    # Horizontal bar (10x150) at y = 375
    addPatch $steelITag $NFY $NFX \
        [expr {0.5*$CW - 75.0}] [expr {0.5*$CD + 125.0 - 5.0}] \
        [expr {0.5*$CW + 75.0}] [expr {0.5*$CD + 125.0 + 5.0}]
    # Vertical bar (10x105) at x = 292.5
    addPatch $steelITag $NFY $NFX \
        [expr {0.5*$CW - 57.5 - 5.0}] [expr {0.5*$CD - 5.0}] \
        [expr {0.5*$CW - 57.5 + 5.0}] [expr {0.5*$CD + 5.0}]
    # Vertical bar (10x105) at x = 407.5
    addPatch $steelITag $NFY $NFX \
        [expr {0.5*$CW + 57.5 - 5.0}] [expr {0.5*$CD - 5.0}] \
        [expr {0.5*$CW + 57.5 + 5.0}] [expr {0.5*$CD + 5.0}]
    # Horizontal bar (150x10) at y = 250, x = 235
    addPatch $steelITag $NFY $NFX \
        [expr {0.5*$CW - 115.0 - 5.0}] [expr {0.5*$CD - 75.0}] \
        [expr {0.5*$CW - 115.0 + 5.0}] [expr {0.5*$CD + 75.0}]
    # Horizontal bar (150x10) at y = 250, x = 465
    addPatch $steelITag $NFY $NFX \
        [expr {0.5*$CW + 115.0 - 5.0}] [expr {0.5*$CD - 75.0}] \
        [expr {0.5*$CW + 115.0 + 5.0}] [expr {0.5*$CD + 75.0}]

    # 3d. Rebars (circular steel fibers)
    set rebars {
        {25.0  50.0   50.0}
        {18.0  50.0  350.0}
        {25.0  50.0  650.0}
        {25.0 450.0   50.0}
        {18.0 450.0  350.0}
        {25.0 450.0  650.0}
        {16.0 250.0   50.0}
        {16.0 250.0  650.0}
        {16.0 125.0   50.0}
        {16.0 125.0  650.0}
        {16.0 375.0   50.0}
        {16.0 375.0  650.0}
        {18.0  50.0  175.0}
        {18.0  50.0  525.0}
        {18.0 450.0  175.0}
        {18.0 450.0  525.0}
    }
    foreach bar $rebars {
        lassign $bar dia y x
        set area [expr { 3.141592653589793 * $dia * $dia / 4.0 }]
        fiber $x $y $area $steelTag
    }

    # ------------------------------------------------------------------
    # 4. Compute section height and mass per unit length
    # ------------------------------------------------------------------
    set sectionHeight [expr { $maxTop - $minBottom }]
    # Mass per mm in N·s²/mm  (CONCRETE_DENSITY in kg/m³, area in mm²)
    set eleMass [expr { $CONCRETE_DENSITY * $totalArea * 1.0e-12 }]

    puts "Section height = $sectionHeight mm"
    puts "Section area   = $totalArea mm²"
    puts "Element mass per mm = $eleMass N·s²/mm"

    return [list $sectionHeight $eleMass]
}

# Define a section with tag 1, inelastic steel, concrete strength -30 MPa,
# confinement factor 1.2, density 2400 kg/m³.
#set result [COMPOSITE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_TWO 1 "INELASTIC" -30.0 1.2 2400.0]
#set H [lindex $result 0]
#set M [lindex $result 1]