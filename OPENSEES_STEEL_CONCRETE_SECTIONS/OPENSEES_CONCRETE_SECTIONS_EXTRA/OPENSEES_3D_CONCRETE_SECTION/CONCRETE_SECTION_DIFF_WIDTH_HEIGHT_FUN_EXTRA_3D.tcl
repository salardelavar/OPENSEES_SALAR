# CONCRETE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D
# Create a reinforced concrete fiber section (OpenSees) with multiple concrete layers
# and discrete rebars. Returns section height (mm) and mass per unit length (kg/mm).
#
# Arguments:
#   secTag            - section identifier (integer)
#   STEEL_TYPE        - "ELASTIC" or "INELASTIC"
#   fc                - unconfined concrete compressive strength (MPa, positive magnitude)
#   Kfc               - ratio of confined to unconfined concrete strength (not used in patches)
#   CONCRETE_DENSITY  - density in kg/m³ (used directly: density * area -> kg/mm)
#
# NOTE: The mass calculation follows the original Python code (only concrete area, no rebar mass)
#       and is physically incorrect (should be multiplied by 1e-9), but kept for compatibility.
# THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)

proc CONCRETE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D {secTag STEEL_TYPE fc Kfc CONCRETE_DENSITY} {
    # ------------------------------------------------------------------
    #  Material properties
    # ------------------------------------------------------------------
    # Concrete (units: mm, N)
    set fc_abs [expr {abs($fc)}]                     ; # positive magnitude
    # Cover (unconfined) – all concrete layers use this material
    set fpc_cover [expr {-$fc_abs}]                  ; # negative compressive strength
    set Ec_cover  [expr {4700.0 * sqrt($fc_abs)}]
    set ec0_cover [expr {2.0 * $fpc_cover / $Ec_cover}]   ; # negative strain
    set fpcu_cover [expr {0.2 * $fpc_cover}]              ; # negative
    set ecu_cover  [expr {5.0 * $ec0_cover}]              ; # negative
    set lambda_cover 0.1
    set ft_cover   [expr {0.7 * sqrt($fc_abs)}]           ; # positive tensile strength
    set Ets_cover  [expr {$ft_cover / abs($ec0_cover)}]   ; # tension softening stiffness

    # Core (confined) – defined but NOT used in patches (kept for compatibility)
    set fpc_core [expr {-$Kfc * $fc_abs}]
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
    set coreTag   [expr {$secTag + 100}]   ; # confined (not used in patches)
    set coverTag  [expr {$secTag + 200}]   ; # unconfined (used for all concrete layers)
    set steelTag  [expr {$secTag + 300}]   ; # rebar steel

    # ------------------------------------------------------------------
    #  Create uniaxial materials
    # ------------------------------------------------------------------
    # Concrete (cover only, since no core patches)
    uniaxialMaterial Concrete02 $coverTag $fpc_cover $ec0_cover $fpcu_cover $ecu_cover $lambda_cover $ft_cover $Ets_cover

    # Steel for rebars
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
    #  Fiber section definition
    # ------------------------------------------------------------------
    section Fiber $secTag -GJ 1.0e7

    # ------------------------------------------------------------------
    #  Concrete layers (20 layers, all using coverTag)
    #  Format: depth, width, center_y, center_x, nfy, nfx
    # ------------------------------------------------------------------
    set concrete_layers {
        {50.0 700.0  25.0 350.0 10 10}
        {50.0 680.0  75.0 350.0 10 10}
        {50.0 660.0 125.0 350.0 10 10}
        {50.0 640.0 175.0 350.0 10 10}
        {50.0 620.0 225.0 350.0 10 10}
        {50.0 600.0 275.0 350.0 10 10}
        {50.0 580.0 325.0 350.0 10 10}
        {50.0 560.0 375.0 350.0 10 10}
        {50.0 540.0 425.0 350.0 10 10}
        {50.0 520.0 475.0 350.0 10 10}
        {50.0 500.0 525.0 350.0 10 10}
        {50.0 520.0 575.0 350.0 10 10}
        {50.0 540.0 625.0 350.0 10 10}
        {50.0 560.0 675.0 350.0 10 10}
        {50.0 580.0 725.0 350.0 10 10}
        {50.0 620.0 775.0 350.0 10 10}
        {50.0 640.0 825.0 350.0 10 10}
        {50.0 660.0 875.0 350.0 10 10}
        {50.0 680.0 925.0 350.0 10 10}
        {50.0 700.0 975.0 350.0 10 10}
    }

    set min_bottom 1e9
    set max_top   -1e9
    set total_area 0.0

    foreach layer $concrete_layers {
        lassign $layer depth width cy cx nfy nfx
        set x_left [expr {$cx - $width/2.0}]
        set x_right [expr {$cx + $width/2.0}]
        set y_bot  [expr {$cy - $depth/2.0}]
        set y_top  [expr {$cy + $depth/2.0}]

        # Update extents
        if {$y_bot < $min_bottom} { set min_bottom $y_bot }
        if {$y_top > $max_top}    { set max_top    $y_top }

        # Accumulate area (concrete only)
        set total_area [expr {$total_area + $depth * $width}]

        # Create rectangular patch
        patch rect $coverTag $nfy $nfx $x_left $y_bot $x_right $y_top
    }

    set section_height [expr {$max_top - $min_bottom}]
    puts "Section height = $section_height mm"
    puts "Section area   = $total_area mm²"

    # ------------------------------------------------------------------
    #  Reinforcing bars (fibers)
    # ------------------------------------------------------------------
    # rebar list: diameter (mm), y-coord (mm), x-coord (mm)
    set rebars {
        {25.0  50.0  50.0}
        {25.0  50.0 350.0}
        {25.0  50.0 650.0}
        {18.0 100.0  70.0}
        {18.0 100.0 350.0}
        {18.0 100.0 630.0}
        {18.0 150.0  90.0}
        {18.0 150.0 350.0}
        {18.0 150.0 610.0}
        {18.0 200.0 110.0}
        {18.0 200.0 350.0}
        {18.0 200.0 590.0}
        {25.0 950.0  50.0}
        {25.0 950.0 350.0}
        {25.0 950.0 650.0}
        {18.0 900.0  70.0}
        {18.0 900.0 350.0}
        {18.0 900.0 630.0}
        {18.0 850.0  90.0}
        {18.0 850.0 350.0}
        {18.0 850.0 610.0}
        {18.0 800.0 110.0}
        {18.0 800.0 350.0}
        {18.0 800.0 590.0}
        {16.0 500.0 150.0}
        {16.0 500.0 350.0}
        {16.0 500.0 550.0}
    }

    foreach rebar $rebars {
        lassign $rebar dia y x
        set area [expr {3.141592653589793 * $dia * $dia / 4.0}]
        fiber $x $y $area $steelTag
    }

    # ------------------------------------------------------------------
    #  Mass per unit length (kg/mm) – following original Python formula
    #  (only concrete area, rebar mass neglected)
    # ------------------------------------------------------------------
    set ele_mass [expr {$CONCRETE_DENSITY * $total_area}]

    return [list $section_height $ele_mass]
}

#set secTag 1
#set steelType "INELASTIC"
#set fc 30.0          ; # concrete compressive strength (MPa, positive)
#set Kfc 1.3          ; # confinement factor (unused in patches, kept for compatibility)
#set density 2400.0   ; # concrete density (kg/m³)

#set result [CONCRETE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D $secTag $steelType $fc $Kfc $density]
#set height [lindex $result 0]
#set mass   [lindex $result 1]
#puts "Height = $height mm, Mass = $mass kg/mm"