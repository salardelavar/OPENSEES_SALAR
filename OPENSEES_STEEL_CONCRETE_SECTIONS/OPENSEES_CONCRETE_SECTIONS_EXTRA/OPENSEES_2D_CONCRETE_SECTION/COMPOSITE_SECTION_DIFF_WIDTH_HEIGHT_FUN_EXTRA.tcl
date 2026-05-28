# COMPOSITE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA
# Create a composite steel-concrete fiber section (OpenSees) with confined core,
# unconfined cover, steel plates, and rebar. Returns section height (mm) and
# mass per unit length (kg/mm).
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
# THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
proc CompositeSectionDiffWidthHeightFunExtra {secTag STEEL_TYPE fc Kfc CONCRETE_DENSITY} {
    # ------------------------------------------------------------------
    #  Material properties
    # ------------------------------------------------------------------
    # Concrete (units: mm, N)
    set fc_abs [expr {abs($fc)}]              ; # positive magnitude
    # Cover (unconfined)
    set fpc_cover [expr {-$fc_abs}]           ; # negative compressive strength
    set Ec_cover  [expr {4700.0 * sqrt($fc_abs)}]
    set ec0_cover [expr {2.0 * $fpc_cover / $Ec_cover}]   ; # negative strain
    set fpcu_cover [expr {0.2 * $fpc_cover}]  ; # negative
    set ecu_cover  [expr {5.0 * $ec0_cover}]  ; # negative (more compressive)
    set lambda_cover 0.1
    set ft_cover   [expr {0.7 * sqrt($fc_abs)}]          ; # positive tensile strength
    set Ets_cover  [expr {$ft_cover / abs($ec0_cover)}]  ; # tension softening stiffness

    # Core (confined)
    set fpc_core [expr {-$Kfc * $fc_abs}]     ; # negative, stronger
    set Ec_core  [expr {4700.0 * sqrt($Kfc * $fc_abs)}]
    set ec0_core [expr {2.0 * $fpc_core / $Ec_core}]
    set fpcu_core [expr {0.65 * $fpc_core}]   ; # negative
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

    # Steel for plates / I‑section
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
    set steelTag  [expr {$secTag + 300}]   ; # rebar steel
    set steelITag [expr {$secTag + 400}]   ; # plate / I‑section steel

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
    #  Fiber section definition
    # ------------------------------------------------------------------
    section Fiber $secTag

    # ------------------------------------------------------------------
    #  Concrete core patches (N x N grid covering 500x500 mm)
    # ------------------------------------------------------------------
    set CD   500.0      ; # overall depth / width of concrete core (mm)
    set N    7          ; # number of divisions in each direction
    set step [expr {$CD / double($N)}]   ; # ~71.4286 mm
    set nfy  10         ; # number of fiber divisions in Y per patch
    set nfx  10         ; # number of fiber divisions in X per patch

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
            patch rect $coreTag $nfy $nfx $x_left $y_bot $x_right $y_top
        }
    }

    # ------------------------------------------------------------------
    #  Steel plates (perimeter) and cross‑shaped steel profiles
    # ------------------------------------------------------------------
    set PT 10.0        ; # [mm] plate thickness

    # Helper to add a rectangular patch and track extents / area
    # We'll collect patches data for area and height calculation.
    set mat_layers {}   ; # list of lists: depth, width, center_y, center_x, nfy, nfx, matTag
    # (we will fill manually, then later loop to add patches and compute totals)

    # 4 perimeter plates
    lappend mat_layers [list $PT [expr {$CD + 2.0*$PT}] [expr {$CD + 0.5*$PT}] [expr {0.5*$CD}] $nfy $nfx $steelITag]
    lappend mat_layers [list $PT [expr {$CD + 2.0*$PT}] [expr {0.0 - 0.5*$PT}] [expr {0.5*$CD}] $nfy $nfx $steelITag]
    lappend mat_layers [list $CD $PT [expr {0.5*$CD}] [expr {$CD + 0.5*$PT}] $nfy $nfx $steelITag]
    lappend mat_layers [list $CD $PT [expr {0.5*$CD}] [expr {0.0 - 0.5*$PT}] $nfy $nfx $steelITag]

    # Cross‑shaped steel profiles (7 pieces)
    lappend mat_layers [list 10.0 150.0 [expr {250.0 - 125.0}] 250.0 $nfy $nfx $steelITag]
    lappend mat_layers [list 240.0 10.0 250.0 250.0 $nfy $nfx $steelITag]
    lappend mat_layers [list 10.0 150.0 [expr {250.0 + 125.0}] 250.0 $nfy $nfx $steelITag]
    lappend mat_layers [list 10.0 105.0 250.0 [expr {250.0 - 57.5}] $nfy $nfx $steelITag]
    lappend mat_layers [list 10.0 105.0 250.0 [expr {250.0 + 57.5}] $nfy $nfx $steelITag]
    lappend mat_layers [list 150.0 10.0 [expr {250.0 - 115.0}] 250.0 $nfy $nfx $steelITag]
    lappend mat_layers [list 150.0 10.0 [expr {250.0 + 115.0}] 250.0 $nfy $nfx $steelITag]

    # ------------------------------------------------------------------
    #  Add all patches and compute total area and section height
    # ------------------------------------------------------------------
    set min_bottom 1e9
    set max_top   -1e9
    set total_area 0.0

    # First add the NxN core patches (already added above, but we need to account for their area and extents)
    # We'll recompute from the same parameters to avoid duplication of patch creation.
    # Instead of re‑adding, we just compute area and extents for the core grid.
    for {set i 0} {$i < $N} {incr i} {
        for {set j 0} {$j < $N} {incr j} {
            set depth $step
            set width $step
            set center_y [expr {($i + 0.5) * $step}]
            set center_x [expr {($j + 0.5) * $step}]
            set y_bot [expr {$center_y - $depth/2.0}]
            set y_top [expr {$center_y + $depth/2.0}]
            if {$y_bot < $min_bottom} { set min_bottom $y_bot }
            if {$y_top > $max_top}    { set max_top    $y_top }
            set total_area [expr {$total_area + $depth * $width}]
        }
    }

    # Now process the additional steel patches (mat_layers)
    foreach layer $mat_layers {
        lassign $layer depth width cy cx nfy nfx matTag
        set x_left [expr {$cx - $width/2.0}]
        set x_right [expr {$cx + $width/2.0}]
        set y_bot  [expr {$cy - $depth/2.0}]
        set y_top  [expr {$cy + $depth/2.0}]
        # Update extents
        if {$y_bot < $min_bottom} { set min_bottom $y_bot }
        if {$y_top > $max_top}    { set max_top    $y_top }
        set total_area [expr {$total_area + $depth * $width}]
        # Create patch
        patch rect $matTag $nfy $nfx $x_left $y_bot $x_right $y_top
    }

    set section_height [expr {$max_top - $min_bottom}]
    puts "Section height = $section_height mm"
    puts "Section area   = $total_area mm²"

    # ------------------------------------------------------------------
    #  Reinforcing bars (fibers)
    # ------------------------------------------------------------------
    # rebar list: diameter (mm), y-coord (mm), x-coord (mm)
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

#set result [CompositeSectionDiffWidthHeightFunExtra $secTag $steelType $fc $Kfc $density]
#set height [lindex $result 0]
#set mass   [lindex $result 1]
#puts "Height = $height mm, Mass = $mass kg/mm"