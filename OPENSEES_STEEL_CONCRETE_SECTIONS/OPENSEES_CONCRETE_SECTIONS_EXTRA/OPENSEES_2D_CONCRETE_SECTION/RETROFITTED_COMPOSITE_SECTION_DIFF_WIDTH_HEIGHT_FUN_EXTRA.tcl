# RETROFITTED_COMPOSITE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA
# Create a retrofitted composite steel-concrete fiber section (OpenSees) with confined core,
# unconfined cover, steel plates, cross‑shaped steel profiles, and rebar.
# Returns section height (mm) and mass per unit length (kg/mm).
#
# Arguments:
#   secTag            - section identifier (integer)
#   STEEL_TYPE        - "ELASTIC" or "INELASTIC"
#   fc                - unconfined concrete compressive strength (MPa, positive value)
#   Kfc               - ratio of confined to unconfined concrete strength
#   CONCRETE_DENSITY  - density in kg/m³ (used directly: density * area -> kg/mm)
#
# NOTE: The mass calculation follows the original Python code and is physically
#       incorrect (should be multiplied by 1e-9), but kept for compatibility.
# THIS PYTHON SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)

proc RETROFITTED_COMPOSITE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA {secTag STEEL_TYPE fc Kfc CONCRETE_DENSITY} {
    # ------------------------------------------------------------------
    #  Material properties (units: mm, N)
    # ------------------------------------------------------------------
    set fc_abs [expr {abs($fc)}]                     ; # positive magnitude

    # Cover concrete (unconfined)
    set fpc_cover    [expr {-$fc_abs}]
    set Ec_cover     [expr {4700.0 * sqrt($fc_abs)}]
    set ec0_cover    [expr {2.0 * $fpc_cover / $Ec_cover}]
    set fpcu_cover   [expr {0.2 * $fpc_cover}]
    set ecu_cover    [expr {5.0 * $ec0_cover}]
    set lambda_cover 0.1
    set ft_cover     [expr {0.7 * sqrt($fc_abs)}]
    set Ets_cover    [expr {$ft_cover / abs($ec0_cover)}]

    # Core concrete (confined)
    set fpc_core     [expr {-$Kfc * $fc_abs}]
    set Ec_core      [expr {4700.0 * sqrt($Kfc * $fc_abs)}]
    set ec0_core     [expr {2.0 * $fpc_core / $Ec_core}]
    set fpcu_core    [expr {0.65 * $fpc_core}]
    set ecu_core     [expr {15.0 * $ec0_core}]
    set lambda_core  0.1
    set ft_core      [expr {0.7 * sqrt($Kfc * $fc_abs)}]
    set Ets_core     [expr {$ft_core / abs($ec0_core)}]

    # Steel for reinforcing bars
    set fy   400.0
    set Es   200000.0
    set ey   [expr {$fy / $Es}]
    set fu   [expr {1.1818 * $fy}]
    set esu  0.09
    set Esh  [expr {($fu - $fy) / ($esu - $ey)}]
    set Bs   [expr {$Esh / $Es}]

    # Steel for plates / I‑section (retrofit steel)
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
    set steelTag  [expr {$secTag + 300}]
    set steelITag [expr {$secTag + 400}]

    # ------------------------------------------------------------------
    #  Create uniaxial materials
    # ------------------------------------------------------------------
    # Concrete
    uniaxialMaterial Concrete02 $coreTag  $fpc_core $ec0_core $fpcu_core $ecu_core $lambda_core $ft_core $Ets_core
    uniaxialMaterial Concrete02 $coverTag $fpc_cover $ec0_cover $fpcu_cover $ecu_cover $lambda_cover $ft_cover $Ets_cover

    # Steel
    if {$STEEL_TYPE eq "ELASTIC"} {
        uniaxialMaterial Elastic $steelTag  $Es
        uniaxialMaterial Elastic $steelITag $EsI
    } elseif {$STEEL_TYPE eq "INELASTIC"} {
        set pinchX   0.8
        set pinchY   0.5
        set damage1  0.0
        set damage2  0.0
        set beta     0.1
        # Rebar steel
        uniaxialMaterial Hysteretic $steelTag \
                $fy $ey $fu $esu [expr {0.2*$fu}] [expr {1.1*$esu}] \
               -$fy -$ey -$fu -$esu [expr {-0.2*$fu}] [expr {-1.1*$esu}] \
                $pinchX $pinchY $damage1 $damage2 $beta
        # Plate / I‑section steel
        uniaxialMaterial Hysteretic $steelITag \
                $fyI $eyI $fuI $esuI [expr {0.2*$fuI}] [expr {1.1*$esuI}] \
               -$fyI -$eyI -$fuI -$esuI [expr {-0.2*$fuI}] [expr {-1.1*$esuI}] \
                $pinchX $pinchY $damage1 $damage2 $beta
    } else {
        puts "ERROR: STEEL_TYPE must be 'ELASTIC' or 'INELASTIC'"
        return
    }

    # ------------------------------------------------------------------
    #  Fiber section definition (no -GJ, same as original Python)
    # ------------------------------------------------------------------
    section Fiber $secTag

    # ------------------------------------------------------------------
    #  Geometry parameters
    # ------------------------------------------------------------------
    set CD   500.0     ; # overall depth / width of concrete core (mm)
    set PT   10.0      ; # plate thickness (mm)
    set N    7         ; # number of divisions in each direction (core grid)
    set step [expr {$CD / double($N)}]   ; # ~71.4286 mm
    set nfy  10        ; # number of fiber divisions in Y per patch
    set nfx  10        ; # number of fiber divisions in X per patch

    # ------------------------------------------------------------------
    #  Concrete core patches (N x N grid)
    # ------------------------------------------------------------------
    set min_bottom 1e9
    set max_top   -1e9
    set total_area 0.0

    for {set i 0} {$i < $N} {incr i} {
        for {set j 0} {$j < $N} {incr j} {
            set depth $step
            set width $step
            set center_y [expr {($i + 0.5) * $step}]
            set center_x [expr {($j + 0.5) * $step}]
            set x_left [expr {$center_x - $width/2.0}]
            set x_right [expr {$center_x + $width/2.0}]
            set y_bot  [expr {$center_y - $depth/2.0}]
            set y_top  [expr {$center_y + $depth/2.0}]

            # Update extents
            if {$y_bot < $min_bottom} { set min_bottom $y_bot }
            if {$y_top > $max_top}    { set max_top    $y_top }

            # Update total area
            set total_area [expr {$total_area + $depth * $width}]

            # Create patch
            patch rect $coreTag $nfy $nfx $x_left $y_bot $x_right $y_top
        }
    }

    # ------------------------------------------------------------------
    #  Steel plates and cross‑shaped profiles (retrofit elements)
    #  Each entry: depth, width, center_y, center_x, nfy, nfx, material tag
    #  (values taken directly from the active Python mat_layers)
    # ------------------------------------------------------------------
    set steel_patches {
        {10.0 110.0 505.0 455.0 10 10}
        {10.0 110.0  -5.0  45.0 10 10}
        {100.0 10.0 450.0 505.0 10 10}
        {100.0 10.0  50.0  -5.0 10 10}
        {10.0 110.0 505.0  45.0 10 10}
        {10.0 110.0  -5.0 455.0 10 10}
        {100.0 10.0 450.0  -5.0 10 10}
        {100.0 10.0  50.0 505.0 10 10}
        {10.0 400.0 515.0 250.0 10 10}
        {10.0 400.0 -15.0 250.0 10 10}
        {400.0 10.0 250.0 515.0 10 10}
        {400.0 10.0 250.0 -15.0 10 10}
        {10.0 150.0 125.0 250.0 10 10}
        {240.0 10.0 250.0 250.0 10 10}
        {10.0 150.0 375.0 250.0 10 10}
        {10.0 105.0 250.0 192.5 10 10}
        {10.0 105.0 250.0 307.5 10 10}
        {150.0 10.0 135.0 250.0 10 10}
        {150.0 10.0 365.0 250.0 10 10}
    }

    foreach patch $steel_patches {
        lassign $patch depth width cy cx nfy nfx
        set x_left [expr {$cx - $width/2.0}]
        set x_right [expr {$cx + $width/2.0}]
        set y_bot  [expr {$cy - $depth/2.0}]
        set y_top  [expr {$cy + $depth/2.0}]

        # Update extents
        if {$y_bot < $min_bottom} { set min_bottom $y_bot }
        if {$y_top > $max_top}    { set max_top    $y_top }

        # Update total area
        set total_area [expr {$total_area + $depth * $width}]

        # Create patch
        patch rect $steelITag $nfy $nfx $x_left $y_bot $x_right $y_top
    }

    # ------------------------------------------------------------------
    #  Section height and area
    # ------------------------------------------------------------------
    set section_height [expr {$max_top - $min_bottom}]
    puts "Section height = $section_height mm"
    puts "Section area   = $total_area mm²"

    # ------------------------------------------------------------------
    #  Reinforcing bars (fibers)
    # ------------------------------------------------------------------
    # Each rebar: diameter (mm), y-coord (mm), x-coord (mm)
    set rebars {
        {25.0  50.0 50.0}
        {26.0  50.0 250.0}
        {25.0  50.0 450.0}
        {25.0 450.0 50.0}
        {16.0 450.0 250.0}
        {25.0 450.0 450.0}
        {25.0 250.0 50.0}
        {25.0 250.0 450.0}
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

    return [list $section_height $ele_mass]
}

#set secTag 1
#set steelType "INELASTIC"
#set fc 30.0          ; # unconfined concrete strength (positive MPa)
#set Kfc 1.3          ; # confinement factor
#set density 2400.0   ; # concrete density (kg/m³)

#set result [RETROFITTED_COMPOSITE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA $secTag $steelType $fc $Kfc $density]
#set height [lindex $result 0]
#set mass   [lindex $result 1]
#puts "Height = $height mm, Mass = $mass kg/mm"