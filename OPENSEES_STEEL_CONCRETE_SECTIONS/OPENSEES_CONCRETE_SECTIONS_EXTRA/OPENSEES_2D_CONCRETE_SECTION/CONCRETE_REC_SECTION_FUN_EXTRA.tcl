# CONCRETE_REC_SECTION_FUN
# Create a rectangular reinforced concrete fiber section (OpenSees) with unconfined concrete.
# Returns section height (mm) and mass per unit length (kg/mm).
#
# Arguments:
#   secTag            - section identifier (integer)
#   STEEL_TYPE        - "ELASTIC" or "INELASTIC"
#   fc                - unconfined concrete compressive strength (MPa, positive value)
#   Kfc               - ratio of confined to unconfined concrete strength (not used for patch, but kept for compatibility)
#   Bsec, Hsec        - width and height of rectangle (mm)
#   cover             - concrete cover for rebars (mm)
#   nFib              - number of fibers per direction in the rectangular patch
#   CONCRETE_DENSITY  - density in kg/mm³ (e.g., 2.5e-9)
#
# NOTE: The original Python uses only coverTag (unconfined) for the whole section.
#       The coreTag is defined but not used; it is kept here for completeness.
# THIS PROGRAM WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
proc CONCRETE_REC_SECTION_FUN {secTag STEEL_TYPE fc Kfc Bsec Hsec cover nFib CONCRETE_DENSITY} {
    # ------------------------------------------------------------------
    #  Material properties (units: mm, N)
    # ------------------------------------------------------------------
    # Cover concrete (unconfined)
    set fc_abs [expr {abs($fc)}]
    set fpc_cover [expr {-$fc_abs}]                     ; # negative compressive strength
    set Ec_cover  [expr {4700.0 * sqrt($fc_abs)}]
    set ec0_cover [expr {2.0 * $fpc_cover / $Ec_cover}] ; # negative strain
    set fpcu_cover [expr {0.2 * $fpc_cover}]            ; # negative
    set ecu_cover  [expr {5.0 * $ec0_cover}]            ; # negative
    set lambda_cover 0.1
    set ft_cover   [expr {0.7 * sqrt($fc_abs)}]         ; # positive tensile strength
    set Ets_cover  [expr {$ft_cover / abs($ec0_cover)}] ; # tension softening stiffness

    # Core concrete (confined) – defined but not used in patch (kept for compatibility)
    set fpc_core [expr {-$Kfc * $fc_abs}]               ; # negative, stronger
    set Ec_core  [expr {4700.0 * sqrt($Kfc * $fc_abs)}]
    set ec0_core [expr {2.0 * $fpc_core / $Ec_core}]
    set fpcu_core [expr {0.65 * $fpc_core}]
    set ecu_core  [expr {15.0 * $ec0_core}]
    set lambda_core 0.1
    set ft_core   [expr {0.7 * sqrt($Kfc * $fc_abs)}]
    set Ets_core  [expr {$ft_core / abs($ec0_core)}]

    # Steel for reinforcing bars
    set fy   400.0        ; # N/mm²
    set Es   200000.0
    set ey   [expr {$fy / $Es}]
    set fu   [expr {1.1818 * $fy}]
    set esu  0.09
    set Esh  [expr {($fu - $fy) / ($esu - $ey)}]
    set Bs   [expr {$Esh / $Es}]

    # Material tags
    set coreTag   [expr {$secTag + 100}]
    set coverTag  [expr {$secTag + 200}]
    set steelTag  [expr {$secTag + 300}]

    # ------------------------------------------------------------------
    #  Create uniaxial materials
    # ------------------------------------------------------------------
    # Concrete
    uniaxialMaterial Concrete02 $coreTag  $fpc_core $ec0_core $fpcu_core $ecu_core $lambda_core $ft_core $Ets_core
    uniaxialMaterial Concrete02 $coverTag $fpc_cover $ec0_cover $fpcu_cover $ecu_cover $lambda_cover $ft_cover $Ets_cover

    # Steel rebars
    if {$STEEL_TYPE eq "ELASTIC"} {
        uniaxialMaterial Elastic $steelTag $Es
    } elseif {$STEEL_TYPE eq "INELASTIC"} {
        set pinchX   0.8
        set pinchY   0.5
        set damage1  0.0
        set damage2  0.0
        set beta     0.1
        uniaxialMaterial Hysteretic $steelTag \
                $fy $ey $fu $esu [expr {0.2*$fu}] [expr {1.1*$esu}] \
               -$fy -$ey -$fu -$esu [expr {-0.2*$fu}] [expr {-1.1*$esu}] \
                $pinchX $pinchY $damage1 $damage2 $beta
    } else {
        puts "ERROR: STEEL_TYPE must be 'ELASTIC' or 'INELASTIC'"
        return
    }

    # ------------------------------------------------------------------
    #  Section area
    # ------------------------------------------------------------------
    set AREA [expr {$Bsec * $Hsec}]
    puts "Total Section Area = $AREA mm²"

    # ------------------------------------------------------------------
    #  Fiber section definition
    # ------------------------------------------------------------------
    section Fiber $secTag

    # Whole rectangle using cover concrete (unconfined)
    set x_left [expr {-$Bsec / 2.0}]
    set x_right [expr {$Bsec / 2.0}]
    set y_bot [expr {-$Hsec / 2.0}]
    set y_top [expr {$Hsec / 2.0}]
    patch rect $coverTag $nFib $nFib $x_left $y_bot $x_right $y_top

    # ------------------------------------------------------------------
    #  Reinforcing bars (fibers)
    # ------------------------------------------------------------------
    # rebars list: diameter (mm), y-coord (mm), x-coord (mm)
    # coordinates are relative to section centroid
    set rebars {
        {25.0 [expr {$Hsec/2.0 - $cover}] [expr {-$Bsec/3.0}]}
        {25.0 [expr {$Hsec/2.0 - $cover}] [expr {$Bsec/3.0}]}
        {18.0 [expr {$Hsec/2.0 - $cover}] 0.0}
        {18.0 0.0 [expr {-$Bsec/3.0}]}
        {18.0 0.0 [expr {$Bsec/3.0}]}
        {25.0 [expr {-$Hsec/2.0 + $cover}] [expr {-$Bsec/3.0}]}
        {25.0 [expr {-$Hsec/2.0 + $cover}] [expr {$Bsec/3.0}]}
        {18.0 [expr {-$Hsec/2.0 + $cover}] 0.0}
    }

    foreach rebar $rebars {
        lassign $rebar dia y x
        set area [expr {3.141592653589793 * $dia * $dia / 4.0}]
        fiber $x $y $area $steelTag
    }

    # ------------------------------------------------------------------
    #  Mass per unit length (kg/mm) – follows original Python
    # ------------------------------------------------------------------
    set ele_mass [expr {$CONCRETE_DENSITY * $AREA}]

    return [list $Hsec $ele_mass]
}

# Define section properties
#set secTag 1
#set STEEL_TYPE "INELASTIC"
#set fc 30.0               ; # MPa (positive)
#set Kfc 1.24
#set Bsec 400.0            ; # mm
#set Hsec 600.0            ; # mm
#set cover 50.0            ; # mm
#set nFib 10               ; # fibers per direction
#set density 2.5e-9        ; # kg/mm³ (concrete density)

#set result [CONCRETE_REC_SECTION_FUN $secTag $STEEL_TYPE $fc $Kfc $Bsec $Hsec $cover $nFib $density]
#set height [lindex $result 0]
#set mass   [lindex $result 1]
#puts "Height = $height mm, Mass = $mass kg/mm"