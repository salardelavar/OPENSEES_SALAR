# CONCRETE_DUMBBELL_SHEAR_WALL_DIFF_FLANGE_FUN_EXTRA_3D
# Create a concrete dumbbell-shaped shear wall fiber section (OpenSees) with different flange widths.
# Returns section height (mm) and mass per unit length (kg/mm).
#
# Arguments:
#   secTag            - section identifier (integer)
#   STEEL_TYPE        - "ELASTIC" or "INELASTIC"
#   fc                - unconfined concrete compressive strength (MPa, positive value)
#   Kfc               - ratio of confined to unconfined concrete strength (not used in patches, but material defined)
#   bfT               - top flange width (mm)
#   tfT               - top flange thickness (mm)
#   bfB               - bottom flange width (mm)
#   tfB               - bottom flange thickness (mm)
#   h                 - total section height (mm)
#   tw                - web thickness (mm)
#   nFib              - number of fiber divisions per patch (in each direction)
#   CONCRETE_DENSITY  - density in kg/m³ (used directly: density * area -> kg/mm)
#
# NOTE: The mass calculation follows the original Python code and is physically
#       incorrect (should be multiplied by 1e-9), but kept for compatibility.
# THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)

proc CONCRETE_DUMBBELL_SHEAR_WALL_DIFF_FLANGE_FUN_EXTRA_3D {secTag STEEL_TYPE fc Kfc bfT tfT bfB tfB h tw nFib CONCRETE_DENSITY} {
    # ------------------------------------------------------------------
    #  Material properties (units: mm, N)
    # ------------------------------------------------------------------
    set fc_abs $fc   ; # positive magnitude

    # Cover concrete (unconfined) – used for all concrete patches
    set fpc_cover [expr { -$fc_abs }]                     ; # negative compressive strength
    set Ec_cover  [expr {4700.0 * sqrt($fc_abs)}]
    set ec0_cover [expr {2.0 * $fpc_cover / $Ec_cover}]   ; # negative strain
    set fpcu_cover [expr {0.2 * $fpc_cover}]              ; # negative
    set ecu_cover  [expr {5.0 * $ec0_cover}]              ; # negative (more compressive)
    set lambda_cover 0.1
    set ft_cover   [expr {0.7 * sqrt($fc_abs)}]           ; # positive tensile strength
    set Ets_cover  [expr {$ft_cover / 0.002}]             ; # tension softening stiffness

    # Core concrete (confined) – defined but not used in patches (kept for consistency)
    set fpc_core [expr { $Kfc * $fpc_cover }]             ; # negative, stronger
    set Ec_core  [expr {4700.0 * sqrt( -$fpc_core )}]     ; # = 4700 * sqrt(Kfc*fc_abs)
    set ec0_core [expr {2.0 * $fpc_core / $Ec_core}]
    set fpcu_core [expr {0.65 * $fpc_core}]
    set ecu_core  [expr {15.0 * $ec0_core}]
    set lambda_core 0.1
    set ft_core   [expr {0.7 * sqrt( -$fpc_core )}]
    set Ets_core  [expr {$ft_core / 0.002}]

    # Steel for reinforcing bars
    set fy   400.0
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
    # Concrete02 for cover (unconfined) – note: core material also created but not used in patches
    uniaxialMaterial Concrete02 $coverTag $fpc_cover $ec0_cover $fpcu_cover $ecu_cover $lambda_cover $ft_cover $Ets_cover
    uniaxialMaterial Concrete02 $coreTag   $fpc_core   $ec0_core   $fpcu_core   $ecu_core   $lambda_core $ft_core $Ets_core

    # Steel reinforcement
    if {$STEEL_TYPE eq "ELASTIC"} {
        uniaxialMaterial Steel01 $steelTag $fy $Es $Bs
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
    #  Section area and neutral axis
    # ------------------------------------------------------------------
    set hw [expr {$h - $tfT - $tfB}]
    set area_top_flange    [expr {$bfT * $tfT}]
    set area_bottom_flange [expr {$bfB * $tfB}]
    set area_web           [expr {$tw * $hw}]
    set total_area         [expr {$area_top_flange + $area_bottom_flange + $area_web}]
    puts "Total Section Area = $total_area mm²"

    set ybar [expr {($area_top_flange * ($tfB + $hw + 0.5*$tfT) +
                     $area_web        * ($tfB + 0.5*$hw) +
                     $area_bottom_flange * $tfB) / $total_area}]
    puts "Neutral axis in the Y‑direction = $ybar mm"

    # ------------------------------------------------------------------
    #  Fiber section definition
    # ------------------------------------------------------------------
    section Fiber $secTag -GJ 1.0e7

    # Helper to add a rectangular patch
    proc add_rect {matTag y_bot y_top x_left x_right nFib} {
        patch rect $matTag $nFib $nFib $x_left $y_bot $x_right $y_top
    }

    # Concrete patches (all using coverTag – unconfined concrete)
    # Top flange
    add_rect $coverTag [expr {$h - $tfT}] $h [expr {-$bfT/2.0}] [expr {$bfT/2.0}] $nFib
    # Bottom flange
    add_rect $coverTag 0.0 $tfB [expr {-$bfB/2.0}] [expr {$bfB/2.0}] $nFib
    # Web
    add_rect $coverTag $tfB [expr {$tfB + $hw}] [expr {-$tw/2.0}] [expr {$tw/2.0}] $nFib

    # ------------------------------------------------------------------
    #  Reinforcing bars (fibers)
    # ------------------------------------------------------------------
    set cover 50.0
    set RD01 25.0
    set RD02 16.0
    set RD03 20.0

    set rebars {
        [list $RD01 [expr {$h - $cover}] [expr {-$bfT/2.0 + $cover}]]
        [list $RD01 [expr {$h - $cover}] [expr { $bfT/2.0 - $cover}]]
        [list $RD01 [expr {$h - $cover}] 0.0]
        [list $RD01 [expr {$h - $tfT + $cover}] [expr {-$bfT/2.0 + $cover}]]
        [list $RD01 [expr {$h - $tfT + $cover}] [expr { $bfT/2.0 - $cover}]]
        [list $RD01 [expr {$h - $tfT + $cover}] 0.0]
        [list $RD01 [expr {$h - $tfT/2.0}] [expr {-$bfT/2.0 + $cover}]]
        [list $RD01 [expr {$h - $tfT/2.0}] [expr { $bfT/2.0 - $cover}]]

        [list $RD02 [expr {$tfB + $hw*0.1}] [expr {-$tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.1}] [expr { $tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.2}] [expr {-$tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.2}] [expr { $tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.3}] [expr {-$tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.3}] [expr { $tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.4}] [expr {-$tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.4}] [expr { $tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.5}] [expr {-$tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.5}] [expr { $tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.6}] [expr {-$tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.6}] [expr { $tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.7}] [expr {-$tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.7}] [expr { $tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.8}] [expr {-$tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.8}] [expr { $tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.9}] [expr {-$tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.9}] [expr { $tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*1.0}] [expr {-$tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*1.0}] [expr { $tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.0}] [expr {-$tw/4.0}]]
        [list $RD02 [expr {$tfB + $hw*0.0}] [expr { $tw/4.0}]]

        [list $RD03 $cover [expr {-$bfB/2.0 + $cover}]]
        [list $RD03 $cover [expr { $bfB/2.0 - $cover}]]
        [list $RD03 [expr {$tfB - $cover}] [expr {-$bfB/2.0 + $cover}]]
        [list $RD03 [expr {$tfB - $cover}] [expr { $bfB/2.0 - $cover}]]
    }

    foreach rebar $rebars {
        lassign $rebar dia y x
        set area [expr {3.141592653589793 * $dia * $dia / 4.0}]
        fiber $x $y $area $steelTag
    }

    # ------------------------------------------------------------------
    #  Mass per unit length (kg/mm) – following original Python formula
    # ------------------------------------------------------------------
    set ele_mass [expr {$CONCRETE_DENSITY * $total_area}]

    return [list $h $ele_mass]
}

#set secTag 1
#set steelType "INELASTIC"
#set fc 30.0          ; # positive unconfined concrete strength (MPa)
#set Kfc 1.3
#set bfT 600.0
#set tfT 150.0
#set bfB 400.0
#set tfB 120.0
#set h 2000.0
#set tw 200.0
#set nFib 10
#set density 2400.0   ; # kg/m³

#set result [CONCRETE_DUMBBELL_SHEAR_WALL_DIFF_FLANGE_FUN_EXTRA_3D $secTag $steelType $fc $Kfc $bfT $tfT $bfB $tfB $h $tw $nFib $density]
#set height [lindex $result 0]
#set mass   [lindex $result 1]
#puts "Height = $height mm, Mass = $mass kg/mm"