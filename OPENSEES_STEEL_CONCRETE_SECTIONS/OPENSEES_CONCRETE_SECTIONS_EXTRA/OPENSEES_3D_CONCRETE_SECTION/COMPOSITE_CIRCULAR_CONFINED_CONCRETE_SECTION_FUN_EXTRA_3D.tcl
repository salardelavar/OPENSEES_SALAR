proc COMPOSITE_CIRCULAR_CONFINED_CONCRETE_SECTION_FUN_EXTRA_3D {secTag RI RO COVER fc Kfc THICKNESS STEEL_TYPE STEEL_DENSITY CONCRETE_DENSITY} {
    # ----------------------------------------------------------------------
    # Create a circular composite confined concrete section with external steel pipe.
    # Materials and fiber section are built in OpenSees.
    # Returns a list: {SECTION_HEIGHT ELE_MASS}
    #   SECTION_HEIGHT : total outer diameter (mm) = 2*(RO + THICKNESS)
    #   ELE_MASS       : mass per unit length (kg/mm)
    # THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI) 
    # ----------------------------------------------------------------------

    # ----- Concrete material parameters (units: mm, N) --------------------
    set fcU [expr {-$fc}]                ;# compressive strength (negative)
    set Ec  [expr {4700.0 * sqrt(-$fcU)}]
    set ec0U [expr {2.0 * $fcU / $Ec}]   ;# negative
    set fcUU [expr {0.2 * $fcU}]
    set ecuU [expr {5.0 * $ec0U}]
    set ftU  [expr {0.7 * sqrt(-$fcU)}]
    set EtsU [expr {$ftU / abs($ec0U)}]
    set LambdaU 0.1

    # Confined concrete
    set fcC [expr {$Kfc * $fcU}]
    set Ec  [expr {4700.0 * sqrt(-$fcC)}]
    set ec0C [expr {2.0 * $fcC / $Ec}]
    set fcUC [expr {0.65 * $fcC}]
    set ecuC [expr {15.0 * $ec0C}]
    set ftC  [expr {0.7 * sqrt(-$fcC)}]
    set EtsC [expr {$ftC / abs($ec0C)}]
    set LambdaC 0.1

    # ----- Steel material parameters ---------------------------------------
    # Reinforcing steel (rebars)
    set fy   400.0
    set Es   200000.0
    set ey   [expr {$fy / $Es}]
    set fu   [expr {1.1818 * $fy}]
    set esu  0.09
    set Esh  [expr {($fu - $fy) / ($esu - $ey)}]
    set Bs   [expr {$Esh / $Es}]

    # Steel for I‑sections and outer tube
    set fyI  240.0
    set EsI  200000.0
    set eyI  [expr {$fyI / $EsI}]
    set fuI  [expr {1.1818 * $fyI}]
    set esuI 0.25
    set EshI [expr {($fuI - $fyI) / ($esuI - $eyI)}]
    set BsI  [expr {$EshI / $EsI}]

    # Material tags
    set coreTag   [expr {$secTag + 100}]
    set coverTag  [expr {$secTag + 200}]
    set steelTag  [expr {$secTag + 300}]   ;# rebars
    set steelITag [expr {$secTag + 400}]   ;# tube and I‑profiles

    # ----- Define steel materials ------------------------------------------
    if {$STEEL_TYPE == "ELASTIC"} {
        uniaxialMaterial Elastic $steelTag  $Es
        uniaxialMaterial Elastic $steelITag $EsI
    } elseif {$STEEL_TYPE == "INELASTIC"} {
        set pinchX 0.8
        set pinchY 0.5
        set damage1 0.0
        set damage2 0.0
        set beta    0.1
        uniaxialMaterial Hysteretic $steelTag \
            $fy $ey $fu $esu [expr {0.2*$fu}] [expr {1.1*$esu}] \
            -$fy -$ey -$fu -$esu [expr {-0.2*$fu}] [expr {-1.1*$esu}] \
            $pinchX $pinchY $damage1 $damage2 $beta
        uniaxialMaterial Hysteretic $steelITag \
            $fyI $eyI $fuI $esuI [expr {0.2*$fuI}] [expr {1.1*$esuI}] \
            -$fyI -$eyI -$fuI -$esuI [expr {-0.2*$fuI}] [expr {-1.1*$esuI}] \
            $pinchX $pinchY $damage1 $damage2 $beta
    } else {
        puts "WARNING: Unknown STEEL_TYPE '$STEEL_TYPE' – using ELASTIC"
        uniaxialMaterial Elastic $steelTag  $Es
        uniaxialMaterial Elastic $steelITag $EsI
    }

    # ----- Define concrete materials (Concrete02) -------------------------
    uniaxialMaterial Concrete02 $coreTag  $fcC $ec0C $fcUC $ecuC $LambdaC $ftC $EtsC
    uniaxialMaterial Concrete02 $coverTag $fcU $ec0U $fcUU $ecuU $LambdaU $ftU $EtsU

    # ----- Fiber section ---------------------------------------------------
    section Fiber $secTag -GJ 1.0e7

    # Subdivision numbers
    set nfCoreR 8
    set nfCoreT 8
    set nfCoverR 4
    set nfCoverT 8

    # Concrete core (from RI to rc = RO - COVER)
    set rc [expr {$RO - $COVER}]
    patch circ $coreTag  $nfCoreT $nfCoreR 0.0 0.0 $RI $rc 0.0 360.0

    # Concrete cover (from rc to RO)
    patch circ $coverTag $nfCoverT $nfCoverR 0.0 0.0 $rc $RO 0.0 360.0

    # External steel pipe (from RO to RO+THICKNESS)
    set pipeOuter [expr {$RO + $THICKNESS}]
    patch circ $steelITag $nfCoverT $nfCoverR 0.0 0.0 $RO $pipeOuter 0.0 360.0

    # ----- Internal steel I‑profiles (rectangular patches) ----------------
    # Each entry: depth, width, centerY, centerX, nfY, nfX, matTag, colour (colour ignored)
    # Data from Python mat_layers:
    set CD $RO   ;# section depth (diameter)
    # 1: top horizontal bar
    set y1 [expr {0.5 * $CD}]
    set xMin1 [expr {0.0 - 150.0/2.0}]
    set xMax1 [expr {0.0 + 150.0/2.0}]
    set yMin1 [expr {$y1 - 10.0/2.0}]
    set yMax1 [expr {$y1 + 10.0/2.0}]
    patch rect $steelITag 10 10 $xMin1 $yMin1 $xMax1 $yMax1

    # 2: central vertical bar
    set xMin2 [expr {0.0 - 10.0/2.0}]
    set xMax2 [expr {0.0 + 10.0/2.0}]
    set yMin2 [expr {0.0 - 240.0/2.0}]
    set yMax2 [expr {0.0 + 240.0/2.0}]
    patch rect $steelITag 10 10 $xMin2 $yMin2 $xMax2 $yMax2

    # 3: bottom horizontal bar
    set y3 [expr {-0.5 * $CD}]
    set xMin3 [expr {0.0 - 150.0/2.0}]
    set xMax3 [expr {0.0 + 150.0/2.0}]
    set yMin3 [expr {$y3 - 10.0/2.0}]
    set yMax3 [expr {$y3 + 10.0/2.0}]
    patch rect $steelITag 10 10 $xMin3 $yMin3 $xMax3 $yMax3

    # 4: left horizontal bar (at y=0, x=-57.5)
    set x4 -57.5
    set xMin4 [expr {$x4 - 105.0/2.0}]
    set xMax4 [expr {$x4 + 105.0/2.0}]
    set yMin4 [expr {0.0 - 10.0/2.0}]
    set yMax4 [expr {0.0 + 10.0/2.0}]
    patch rect $steelITag 10 10 $xMin4 $yMin4 $xMax4 $yMax4

    # 5: right horizontal bar (at y=0, x=57.5)
    set x5 57.5
    set xMin5 [expr {$x5 - 105.0/2.0}]
    set xMax5 [expr {$x5 + 105.0/2.0}]
    set yMin5 [expr {0.0 - 10.0/2.0}]
    set yMax5 [expr {0.0 + 10.0/2.0}]
    patch rect $steelITag 10 10 $xMin5 $yMin5 $xMax5 $yMax5

    # 6: left vertical bar (at x=-115)
    set x6 -115.0
    set xMin6 [expr {$x6 - 10.0/2.0}]
    set xMax6 [expr {$x6 + 10.0/2.0}]
    set yMin6 [expr {0.0 - 150.0/2.0}]
    set yMax6 [expr {0.0 + 150.0/2.0}]
    patch rect $steelITag 10 10 $xMin6 $yMin6 $xMax6 $yMax6

    # 7: right vertical bar (at x=115)
    set x7 115.0
    set xMin7 [expr {$x7 - 10.0/2.0}]
    set xMax7 [expr {$x7 + 10.0/2.0}]
    set yMin7 [expr {0.0 - 150.0/2.0}]
    set yMax7 [expr {0.0 + 150.0/2.0}]
    patch rect $steelITag 10 10 $xMin7 $yMin7 $xMax7 $yMax7

    # ----- Longitudinal rebars (fibers) -----------------------------------
    # List of {diameter y x}
    set rebars {
        {25.0  200.0000  0.0000}
        {16.0  184.7759  76.5367}
        {25.0    0.0000 200.0000}
        {16.0 -184.7759  76.5367}
        {25.0 -200.0000   0.0000}
        {16.0 -184.7759 -76.5367}
        {25.0    0.0000 -200.0000}
        {16.0  184.7759 -76.5367}
        {16.0  141.4214 -141.4214}
        {16.0   76.5367 -184.7759}
        {16.0  -76.5367 -184.7759}
        {16.0 -141.4214 -141.4214}
        {16.0  -76.5367  184.7759}
        {16.0  141.4214  141.4214}
        {16.0 -141.4214  141.4214}
        {16.0   76.5367  184.7759}
    }

    set pi 3.141592653589793
    foreach bar $rebars {
        lassign $bar dia y x
        set area [expr {$pi * $dia * $dia / 4.0}]
        fiber $x $y $area $steelTag
    }

    # ----- Compute section height and mass per unit length ----------------
    set SECTION_HEIGHT [expr {2.0 * ($RO + $THICKNESS)}]

    # Concrete area = π*(RO^2 - RI^2)   (core + cover)
    set CONCRETE_AREA [expr {$pi * ($RO*$RO - $RI*$RI)}]
    # Steel pipe area = π*((RO+THICKNESS)^2 - RO^2)
    set STEEL_PIPE_AREA [expr {$pi * (($RO+$THICKNESS)*($RO+$THICKNESS) - $RO*$RO)}]

    # Mass per unit length (kg/mm)
    set ELE_MASS [expr {$CONCRETE_AREA * $CONCRETE_DENSITY + $STEEL_PIPE_AREA * $STEEL_DENSITY}]

    puts "Section height = $SECTION_HEIGHT mm"
    puts "Mass per unit length = $ELE_MASS kg/mm"

    # Return both values as a Tcl list
    return [list $SECTION_HEIGHT $ELE_MASS]
}

# Define section tag and parameters
#set secTag 1
#set RI 0.0
#set RO 200.0
#set COVER 40.0
#set fc 30.0
#set Kfc 1.3
#set THICKNESS 10.0
#set STEEL_TYPE "INELASTIC"
#set STEEL_DENSITY 7.85e-6      ;# kg/mm³ (7850 kg/m³)
#set CONCRETE_DENSITY 2.4e-6    ;# kg/mm³ (2400 kg/m³)

#set result [COMPOSITE_CIRCULAR_CONFINED_CONCRETE_SECTION_FUN_EXTRA_3D $secTag $RI $RO $COVER $fc $Kfc $THICKNESS $STEEL_TYPE $STEEL_DENSITY $CONCRETE_DENSITY]
#lassign $result Diameter MassPerLength
#puts "Diameter = $Diameter mm, Mass = $MassPerLength kg/mm"