# PRESTRESS_CONCRETE_I_SECTION_DIFF_FLANGE_FUN
# Create a prestressed concrete I‑section fiber model (OpenSees) with different flange widths.
# Uses Concrete02 for unconfined cover (and confined core, though not applied to any patch),
# Hysteretic or Elastic steel for rebars and prestressing cables.
# Returns the section height (mm) and mass per unit length (kg/mm).
#
# Arguments:
#   secTag            - section identifier (integer)
#   matTag            - (unused, kept for compatibility)
#   STEEL_TYPE        - "ELASTIC" or "INELASTIC"
#   fc                - unconfined concrete compressive strength (MPa, negative)
#   Kfc               - ratio of confined to unconfined concrete strength
#   bfT, tfT          - top flange width and thickness (mm)
#   bfB, tfB          - bottom flange width and thickness (mm)
#   h                 - total section height (mm)
#   tw                - web thickness (mm)
#   nFib              - number of fiber divisions in each direction for rectangular patches
#   CONCRETE_DENSITY  - density in kg/m³ (used directly: density * area -> kg/mm)
#
# NOTE: The mass calculation follows the original Python code (density [kg/m³] * area [mm²])
#       which is physically incorrect (should be multiplied by 1e-9), but kept for compatibility.
#       The prestressing cable materials are defined but no fibers are added, as in the original.
# THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)

proc PRESTRESS_CONCRETE_I_SECTION_DIFF_FLANGE_FUN {secTag matTag STEEL_TYPE fc Kfc bfT tfT bfB tfB h tw nFib CONCRETE_DENSITY} {
    # ------------------------------------------------------------------
    #  Material properties – Concrete
    # ------------------------------------------------------------------
    set fc_pos [expr {abs($fc)}]                    ; # positive magnitude
    # Cover (unconfined)
    set fpc_cover [expr {-$fc_pos}]                 ; # negative compressive strength
    set Ec_cover  [expr {4700.0 * sqrt($fc_pos)}]
    set ec0_cover [expr {2.0 * $fpc_cover / $Ec_cover}]
    set fpcu_cover [expr {0.2 * $fpc_cover}]
    set ecu_cover  [expr {5.0 * $ec0_cover}]
    set lambda_cover 0.1
    set ft_cover   [expr {0.7 * sqrt($fc_pos)}]
    set Ets_cover  [expr {$ft_cover / abs($ec0_cover)}]

    # Core (confined) – defined but not used in patches (same as original)
    set fpc_core   [expr {-$Kfc * $fc_pos}]         ; # negative, stronger
    set Ec_core    [expr {4700.0 * sqrt($Kfc * $fc_pos)}]
    set ec0_core   [expr {2.0 * $fpc_core / $Ec_core}]
    set fpcu_core  [expr {0.65 * $fpc_core}]
    set ecu_core   [expr {15.0 * $ec0_core}]
    set lambda_core 0.1
    set ft_core    [expr {0.7 * sqrt($Kfc * $fc_pos)}]
    set Ets_core   [expr {$ft_core / abs($ec0_core)}]

    # ------------------------------------------------------------------
    #  Steel materials – reinforcing bars
    # ------------------------------------------------------------------
    set fy   400.0
    set Es   200000.0
    set ey   [expr {$fy / $Es}]
    set fu   [expr {1.1818 * $fy}]
    set esu  0.09
    set Esh  [expr {($fu - $fy) / ($esu - $ey)}]
    set Bs   [expr {$Esh / $Es}]

    # Steel for prestressing cables
    set fyC  1200.0
    set EsC  200000.0
    set eyC  [expr {$fyC / $EsC}]
    set fuC  [expr {1.1818 * $fyC}]
    set esuC 0.06
    set EshC [expr {($fuC - $fyC) / ($esuC - $eyC)}]
    set BsC  [expr {$EshC / $EsC}]

    # Common hysteretic parameters
    set pinchX   0.8
    set pinchY   0.5
    set damage1  0.0
    set damage2  0.0
    set beta     0.1

    # Material tags
    set coreTag   [expr {$secTag + 100}]
    set coverTag  [expr {$secTag + 200}]
    set steelTag  [expr {$secTag + 300}]   ; # rebar steel
    set cableTag  [expr {$secTag + 400}]   ; # base cable material
    set cableITag [expr {$secTag + 500}]   ; # cable with initial strain

    # ------------------------------------------------------------------
    #  Create uniaxial materials
    # ------------------------------------------------------------------
    # Concrete
    uniaxialMaterial Concrete02 $coreTag  $fpc_core $ec0_core $fpcu_core $ecu_core $lambda_core $ft_core $Ets_core
    uniaxialMaterial Concrete02 $coverTag $fpc_cover $ec0_cover $fpcu_cover $ecu_cover $lambda_cover $ft_cover $Ets_cover

    # Steel (rebars and cables)
    if {$STEEL_TYPE eq "ELASTIC"} {
        uniaxialMaterial Elastic $steelTag $Es
        uniaxialMaterial Elastic $cableTag $EsC
    } elseif {$STEEL_TYPE eq "INELASTIC"} {
        uniaxialMaterial Hysteretic $steelTag \
                $fy $ey $fu $esu [expr {0.2*$fu}] [expr {1.1*$esu}] \
               -$fy -$ey -$fu -$esu [expr {-0.2*$fu}] [expr {-1.1*$esu}] \
                $pinchX $pinchY $damage1 $damage2 $beta
        uniaxialMaterial Hysteretic $cableTag \
                $fyC $eyC $fuC $esuC [expr {0.2*$fuC}] [expr {1.1*$esuC}] \
               -$fyC -$eyC -$fuC -$esuC [expr {-0.2*$fuC}] [expr {-1.1*$esuC}] \
                $pinchX $pinchY $damage1 $damage2 $beta
    } else {
        puts "ERROR: STEEL_TYPE must be 'ELASTIC' or 'INELASTIC'"
        return
    }

    # Prestressing initial strain material (compression negative)
    set InitialStrain -0.001   ; # mm/mm
    uniaxialMaterial InitStrainMaterial $cableITag $cableTag $InitialStrain

    # ------------------------------------------------------------------
    #  Section properties (area and neutral axis)
    # ------------------------------------------------------------------
    set hw [expr {$h - $tfT - $tfB}]
    set area_top_flange    [expr {$bfT * $tfT}]
    set area_bottom_flange [expr {$bfB * $tfB}]
    set area_web           [expr {$tw * $hw}]
    set total_area         [expr {$area_top_flange + $area_bottom_flange + $area_web}]
    puts "Total Section Area = [format "%.2f" $total_area] mm²"

    set ybar [expr {($area_top_flange * ($tfB + $hw + 0.5*$tfT) +
                     $area_web       * ($tfB + 0.5*$hw) +
                     $area_bottom_flange * $tfB) / $total_area}]
    puts "Neutral axis in the Y‑direction = [format "%.2f" $ybar] mm"

    # ------------------------------------------------------------------
    #  Fiber section
    # ------------------------------------------------------------------
    section Fiber $secTag

    # Helper to add a rectangular patch
    # (uses global nFib for both directions)
    proc addRect {matTag y_bot y_top x_left x_right {nFibVar 0}} {
        upvar nFib nFibVal
        if {$nFibVar == 0} { set nFibVal $nFibVal }
        patch rect $matTag $nFibVal $nFibVal $x_left $y_bot $x_right $y_top
    }

    # Concrete patches (all using unconfined cover material, as in original)
    # Top flange
    addRect $coverTag [expr {$h - $tfT}] $h [expr {-$bfT/2.0}] [expr {$bfT/2.0}]
    # Bottom flange
    addRect $coverTag 0.0 $tfB [expr {-$bfB/2.0}] [expr {$bfB/2.0}]
    # Web
    addRect $coverTag $tfB [expr {$h - $tfT}] [expr {-$tw/2.0}] [expr {$tw/2.0}]

    # ------------------------------------------------------------------
    #  Reinforcing bars (fibers)
    # ------------------------------------------------------------------
    # Each rebar: diameter (mm), y-coord (mm), x-coord (mm)
    set rebars {
        {16.0 [expr {$h - 1.0*$tfT/3.0}] [expr {-$bfT/2.0 + 0.05*$bfT}]}
        {16.0 [expr {$h - 1.0*$tfT/3.0}] [expr { $bfT/2.0 - 0.05*$bfT}]}
        {16.0 [expr {$h - 1.0*$tfT/3.0}] [expr {-$bfT/2.0 + 0.20*$bfT}]}
        {16.0 [expr {$h - 1.0*$tfT/3.0}] [expr { $bfT/2.0 - 0.20*$bfT}]}
        {16.0 [expr {$h - 1.0*$tfT/3.0}] [expr {-$bfT/2.0 + 0.60*$bfT}]}
        {16.0 [expr {$h - 1.0*$tfT/3.0}] [expr { $bfT/2.0 - 0.60*$bfT}]}

        {16.0 [expr {$h - 2.0*$tfT/3.0}] [expr {-$bfT/2.0 + 0.05*$bfT}]}
        {16.0 [expr {$h - 2.0*$tfT/3.0}] [expr { $bfT/2.0 - 0.05*$bfT}]}
        {16.0 [expr {$h - 2.0*$tfT/3.0}] [expr {-$bfT/2.0 + 0.20*$bfT}]}
        {16.0 [expr {$h - 2.0*$tfT/3.0}] [expr { $bfT/2.0 - 0.20*$bfT}]}
        {16.0 [expr {$h - 2.0*$tfT/3.0}] [expr {-$bfT/2.0 + 0.60*$bfT}]}
        {16.0 [expr {$h - 2.0*$tfT/3.0}] [expr { $bfT/2.0 - 0.60*$bfT}]}

        {12.0 [expr {$tfB + 0.05*$hw}] [expr {-$tw/4.0}]}
        {12.0 [expr {$tfB + 0.05*$hw}] [expr { $tw/4.0}]}
        {12.0 [expr {$tfB + 0.15*$hw}] [expr {-$tw/4.0}]}
        {12.0 [expr {$tfB + 0.15*$hw}] [expr { $tw/4.0}]}
        {12.0 [expr {$tfB + 0.25*$hw}] [expr {-$tw/4.0}]}
        {12.0 [expr {$tfB + 0.25*$hw}] [expr { $tw/4.0}]}
        {12.0 [expr {$tfB + 0.35*$hw}] [expr {-$tw/4.0}]}
        {12.0 [expr {$tfB + 0.35*$hw}] [expr { $tw/4.0}]}

        {12.0 [expr {$tfB + 0.45*$hw}] [expr {-$tw/4.0}]}
        {12.0 [expr {$tfB + 0.45*$hw}] [expr { $tw/4.0}]}
        {12.0 [expr {$tfB + 0.55*$hw}] [expr {-$tw/4.0}]}
        {12.0 [expr {$tfB + 0.55*$hw}] [expr { $tw/4.0}]}
        {12.0 [expr {$tfB + 0.65*$hw}] [expr {-$tw/4.0}]}
        {12.0 [expr {$tfB + 0.65*$hw}] [expr { $tw/4.0}]}
        {12.0 [expr {$tfB + 0.75*$hw}] [expr {-$tw/4.0}]}
        {12.0 [expr {$tfB + 0.75*$hw}] [expr { $tw/4.0}]}
        {12.0 [expr {$tfB + 0.85*$hw}] [expr {-$tw/4.0}]}
        {12.0 [expr {$tfB + 0.85*$hw}] [expr { $tw/4.0}]}
        {12.0 [expr {$tfB + 0.95*$hw}] [expr {-$tw/4.0}]}
        {12.0 [expr {$tfB + 0.95*$hw}] [expr { $tw/4.0}]}

        {16.0 [expr {1.0*$tfB/4.0}] [expr {-$bfB/2.0 + 0.10*$bfB}]}
        {16.0 [expr {1.0*$tfB/4.0}] [expr { $bfB/2.0 - 0.10*$bfB}]}
        {16.0 [expr {1.0*$tfB/4.0}] [expr {-$bfB/2.0 + 0.35*$bfB}]}
        {16.0 [expr {1.0*$tfB/4.0}] [expr { $bfB/2.0 - 0.35*$bfB}]}
        {16.0 [expr {2.0*$tfB/4.0}] [expr {-$bfB/2.0 + 0.10*$bfB}]}
        {16.0 [expr {2.0*$tfB/4.0}] [expr { $bfB/2.0 - 0.10*$bfB}]}
        {16.0 [expr {2.0*$tfB/4.0}] [expr {-$bfB/2.0 + 0.35*$bfB}]}
        {16.0 [expr {2.0*$tfB/4.0}] [expr { $bfB/2.0 - 0.35*$bfB}]}
        {16.0 [expr {3.0*$tfB/4.0}] [expr {-$bfB/2.0 + 0.10*$bfB}]}
        {16.0 [expr {3.0*$tfB/4.0}] [expr { $bfB/2.0 - 0.10*$bfB}]}
        {16.0 [expr {3.0*$tfB/4.0}] [expr {-$bfB/2.0 + 0.35*$bfB}]}
        {16.0 [expr {3.0*$tfB/4.0}] [expr { $bfB/2.0 - 0.35*$bfB}]}
    }

    foreach rebar $rebars {
        lassign $rebar dia y x
        set area [expr {3.141592653589793 * $dia * $dia / 4.0}]
        fiber $x $y $area $steelTag
    }

    # ------------------------------------------------------------------
    #  Prestressing cables – materials defined but no fibers added (as in original)
    # ------------------------------------------------------------------
    set cables {
        {25.0 [expr {1.0*$tfB/4.0}] 0.0}
        {40.0 [expr {2.0*$tfB/4.0}] 0.0}
        {25.0 [expr {3.0*$tfB/4.0}] 0.0}
    }
    # (No fibers added for cables – kept for compatibility)

    # ------------------------------------------------------------------
    #  Mass per unit length (kg/mm) – following original Python formula
    # ------------------------------------------------------------------
    set ele_mass [expr {$CONCRETE_DENSITY * $total_area}]

    return [list $h $ele_mass]
}

#set secTag 1
#set matTag 0          ; # unused
#set steelType "INELASTIC"
#set fc -30.0          ; # unconfined concrete strength (negative, MPa)
#set Kfc 1.3
#set bfT 400.0
#set tfT 80.0
#set bfB 300.0
#set tfB 100.0
#set h 800.0
#set tw 50.0
#set nFib 10
#set density 2400.0    ; # kg/m³

#set result [PRESTRESS_CONCRETE_I_SECTION_DIFF_FLANGE_FUN $secTag $matTag $steelType $fc $Kfc $bfT $tfT $bfB $tfB $h $tw $nFib $density]
#set height [lindex $result 0]
#set mass   [lindex $result 1]
#puts "Height = $height mm, Mass = $mass kg/mm"