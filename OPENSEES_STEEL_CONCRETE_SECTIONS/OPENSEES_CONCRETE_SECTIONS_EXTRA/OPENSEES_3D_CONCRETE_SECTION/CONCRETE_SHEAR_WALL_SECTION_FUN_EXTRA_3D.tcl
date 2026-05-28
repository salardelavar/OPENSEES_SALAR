# CONCRETE_SHEAR_WALL_SECTION_FUN_EXTRA_EXTRA_3D
# Create a rectangular reinforced concrete shear wall fiber section (OpenSees)
# using Concrete02 for concrete and Hysteretic or Elastic steel for rebars.
#
# Arguments:
#   secTag            - section identifier (integer)
#   STEEL_TYPE        - "ELASTIC" or "INELASTIC"
#   fc                - unconfined concrete compressive strength (MPa, positive)
#   Kfc               - ratio of confined to unconfined concrete strength (not used for cover-only section, but kept for compatibility)
#   Bsec, Hsec        - width and height of the rectangle (mm)
#   cover             - concrete cover for rebars (mm)
#   nFib              - number of fiber divisions in the X direction for the patch
#   CONCRETE_DENSITY  - concrete density in kg/mm³ (e.g., 2.5e-9)
#
# Returns:
#   A Tcl list containing: [section_height_mm mass_per_unit_length_kg_per_mm]
# THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)

proc CONCRETE_SHEAR_WALL_SECTION_FUN_EXTRA_EXTRA_3D {secTag STEEL_TYPE fc Kfc Bsec Hsec cover nFib CONCRETE_DENSITY} {
    # ------------------------------------------------------------------
    #  Material properties (units: mm, N)
    # ------------------------------------------------------------------
    # Concrete cover (unconfined)
    set fc_abs [expr {abs($fc)}]                ; # positive magnitude
    set fpc_cover [expr {-$fc_abs}]             ; # negative compressive strength
    set Ec_cover  [expr {4700.0 * sqrt($fc_abs)}]
    set ec0_cover [expr {2.0 * $fpc_cover / $Ec_cover}]
    set fpcu_cover [expr {0.2 * $fpc_cover}]
    set ecu_cover  [expr {5.0 * $ec0_cover}]
    set lambda_cover 0.1
    set ft_cover   [expr {0.7 * sqrt($fc_abs)}] ; # tensile strength (positive)
    set Ets_cover  [expr {$ft_cover / abs($ec0_cover)}]

    # Concrete core (confined) – defined but not used (kept for compatibility)
    set fpc_core [expr {-$Kfc * $fc_abs}]
    set Ec_core  [expr {4700.0 * sqrt($Kfc * $fc_abs)}]
    set ec0_core [expr {2.0 * $fpc_core / $Ec_core}]
    set fpcu_core [expr {0.65 * $fpc_core}]
    set ecu_core  [expr {15.0 * $ec0_core}]
    set lambda_core 0.1
    set ft_core   [expr {0.7 * sqrt($Kfc * $fc_abs)}]
    set Ets_core  [expr {$ft_core / abs($ec0_core)}]

    # Reinforcing steel
    set fy   400.0        ; # N/mm²
    set Es   200000.0
    set ey   [expr {$fy / $Es}]
    set fu   [expr {1.1818 * $fy}]
    set esu  0.09
    set Esh  [expr {($fu - $fy) / ($esu - $ey)}]
    set Bs   [expr {$Esh / $Es}]

    # Material tags
    set coreTag  [expr {$secTag + 100}]
    set coverTag [expr {$secTag + 200}]
    set steelTag [expr {$secTag + 300}]

    # ------------------------------------------------------------------
    #  Create uniaxial materials
    # ------------------------------------------------------------------
    # Concrete (only cover is used in the section)
    uniaxialMaterial Concrete02 $coverTag $fpc_cover $ec0_cover $fpcu_cover $ecu_cover $lambda_cover $ft_cover $Ets_cover
    # Core concrete is defined but not used (original Python defines both but only coverTag is patched)
    uniaxialMaterial Concrete02 $coreTag $fpc_core $ec0_core $fpcu_core $ecu_core $lambda_core $ft_core $Ets_core

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
    #  Section area (gross concrete area)
    # ------------------------------------------------------------------
    set AREA [expr {$Bsec * $Hsec}]
    puts "Total Section Area = $AREA mm²"

    # ------------------------------------------------------------------
    #  Fiber section definition
    # ------------------------------------------------------------------
    section Fiber $secTag -GJ 1.0e7

    # Number of fiber divisions in Y direction (hardcoded to 10, as in original)
    set nFibZ 10

    # Add concrete cover patch covering the whole rectangle
    set x_left [expr {-$Bsec / 2.0}]
    set x_right [expr {$Bsec / 2.0}]
    set y_bot  [expr {-$Hsec / 2.0}]
    set y_top  [expr {$Hsec / 2.0}]
    patch rect $coverTag $nFibZ $nFib $x_left $y_bot $x_right $y_top

    # ------------------------------------------------------------------
    #  Reinforcing bars (fibers)
    # ------------------------------------------------------------------
    set DIST 220.0   ; # spacing between successive bar positions (mm)
    # Helper: add a bar at given coordinates with given diameter
    proc add_bar {x y dia steelTag} {
        set area [expr {3.141592653589793 * $dia * $dia / 4.0}]
        fiber $x $y $area $steelTag
    }

    # Generate bars on the right side (x = Bsec/2 - cover)
    set x_right [expr {$Bsec/2.0 - $cover}]
    # y positions: from -Hsec/2+cover to Hsec/2-cover in steps of DIST, plus zero
    for {set i 0} {$i <= 12} {incr i} {
        set y_neg [expr {-$Hsec/2.0 + $cover + $i * $DIST}]
        set y_pos [expr {$Hsec/2.0 - $cover - $i * $DIST}]
        if {$y_neg <= $y_pos} {
            add_bar $x_right $y_neg 16.0 $steelTag
            if {$y_neg != $y_pos} {
                add_bar $x_right $y_pos 16.0 $steelTag
            }
        }
    }
    # add the bar at y = 0
    add_bar $x_right 0.0 16.0 $steelTag

    # Generate bars on the center line (x = 0)
    set x_center 0.0
    for {set i 0} {$i <= 12} {incr i} {
        set y_neg [expr {-$Hsec/2.0 + $cover + $i * $DIST}]
        set y_pos [expr {$Hsec/2.0 - $cover - $i * $DIST}]
        if {$y_neg <= $y_pos} {
            add_bar $x_center $y_neg 16.0 $steelTag
            if {$y_neg != $y_pos} {
                add_bar $x_center $y_pos 16.0 $steelTag
            }
        }
    }
    add_bar $x_center 0.0 16.0 $steelTag

    # Generate bars on the left side (x = -Bsec/2 + cover)
    set x_left [expr {-$Bsec/2.0 + $cover}]
    for {set i 0} {$i <= 12} {incr i} {
        set y_neg [expr {-$Hsec/2.0 + $cover + $i * $DIST}]
        set y_pos [expr {$Hsec/2.0 - $cover - $i * $DIST}]
        if {$y_neg <= $y_pos} {
            add_bar $x_left $y_neg 16.0 $steelTag
            if {$y_neg != $y_pos} {
                add_bar $x_left $y_pos 16.0 $steelTag
            }
        }
    }
    add_bar $x_left 0.0 16.0 $steelTag

    # ------------------------------------------------------------------
    #  Mass per unit length (kg/mm) – follows original formula
    # ------------------------------------------------------------------
    set ele_mass [expr {$CONCRETE_DENSITY * $AREA}]

    return [list $Hsec $ele_mass]
}

#set secTag 1
#set steelType "INELASTIC"
#set fc 30.0              ; # MPa
#set Kfc 1.3              ; # confinement factor (not used for cover-only section)
#set Bsec 1000.0          ; # mm
#set Hsec 2000.0          ; # mm
#set cover 40.0           ; # mm
#set nFib 10              ; # number of fibers in X direction per patch
#set density 2.5e-9       ; # kg/mm³

#set result [CONCRETE_SHEAR_WALL_SECTION_FUN_EXTRA_3D $secTag $steelType $fc $Kfc $Bsec $Hsec $cover $nFib $density]
#set height [lindex $result 0]
#set mass   [lindex $result 1]
#puts "Height = $height mm, Mass = $mass kg/mm"