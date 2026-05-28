# COMPOSITE_REC_CROSS_SHAPED_STEEL_SECTION_FUN_EXTRA
# Create a composite rectangular section with cross‑shaped steel profile,
# confined concrete core, unconfined cover, and reinforcing bars.
#
# Arguments:
#   secTag            - section identifier (integer)
#   HSec, BSec        - height and width of concrete section (mm)
#   cover             - concrete cover to reinforcement (mm)
#   STEEL_TYPE        - "ELASTIC" or "INELASTIC"
#   fc                - unconfined concrete compressive strength (positive MPa)
#   Kc                - ratio of confined to unconfined concrete strength
#   STEEL_DENSITY     - density of steel (kg/m³)
#   CONCRETE_DENSITY  - density of concrete (kg/m³)
#   PLOT              - boolean (ignored in Tcl, kept for compatibility)
#
# Returns:
#   A Tcl list containing: [overall_section_height_mm mass_per_unit_length_kg_per_mm]
# THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)

proc COMPOSITE_REC_CROSS_SHAPED_STEEL_SECTION_FUN_EXTRA {secTag HSec BSec cover STEEL_TYPE fc Kc STEEL_DENSITY CONCRETE_DENSITY PLOT} {
    # ------------------------------------------------------------------
    #  Geometry parameters
    # ------------------------------------------------------------------
    set halfH [expr {$HSec / 2.0}]
    set halfB [expr {$BSec / 2.0}]
    set NUMFIBERS 20
    set PT 10.0                     ; # plate thickness for external steel plates
    set REBAR_DIA 25.0              ; # main rebar diameter (mm)

    # ------------------------------------------------------------------
    #  Concrete material properties (units: mm, N)
    # ------------------------------------------------------------------
    set fc_abs [expr {abs($fc)}]                ; # positive magnitude
    # Unconfined (cover)
    set fpcU [expr {-$fc_abs}]                  ; # negative
    set EcU  [expr {4700.0 * sqrt($fc_abs)}]
    set ec0U [expr {2.0 * $fpcU / $EcU}]
    set fpcuU [expr {0.2 * $fpcU}]
    set ecuU [expr {5.0 * $ec0U}]
    set lambdaU 0.1
    set ftU  [expr {0.7 * sqrt($fc_abs)}]
    set EtsU [expr {$ftU / abs($ec0U)}]

    # Confined (core)
    set fpcC [expr {-$Kc * $fc_abs}]            ; # negative, stronger
    set EcC  [expr {4700.0 * sqrt($Kc * $fc_abs)}]
    set ec0C [expr {2.0 * $fpcC / $EcC}]
    set fpcuC [expr {0.65 * $fpcC}]
    set ecuC [expr {15.0 * $ec0C}]
    set lambdaC 0.1
    set ftC  [expr {0.7 * sqrt($Kc * $fc_abs)}]
    set EtsC [expr {$ftC / abs($ec0C)}]

    # ------------------------------------------------------------------
    #  Steel material properties
    # ------------------------------------------------------------------
    # Reinforcing steel
    set fy   400.0
    set Es   200000.0
    set ey   [expr {$fy / $Es}]
    set fu   [expr {1.1818 * $fy}]
    set esu  0.09
    # Plate / I‑section steel
    set fyI  240.0
    set EsI  200000.0
    set eyI  [expr {$fyI / $EsI}]
    set fuI  [expr {1.1818 * $fyI}]
    set esuI 0.25

    # Material tags
    set coreTag   [expr {$secTag + 100}]
    set coverTag  [expr {$secTag + 200}]
    set steelTag  [expr {$secTag + 300}]       ; # rebar
    set steelITag [expr {$secTag + 400}]       ; # steel plate / I‑section

    # ------------------------------------------------------------------
    #  Create uniaxial materials
    # ------------------------------------------------------------------
    # Concrete
    uniaxialMaterial Concrete02 $coreTag  $fpcC $ec0C $fpcuC $ecuC $lambdaC $ftC $EtsC
    uniaxialMaterial Concrete02 $coverTag $fpcU $ec0U $fpcuU $ecuU $lambdaU $ftU $EtsU

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
        # Rebar
        uniaxialMaterial Hysteretic $steelTag \
                $fy $ey $fu $esu [expr {0.2*$fu}] [expr {1.1*$esu}] \
               -$fy -$ey -$fu -$esu [expr {-0.2*$fu}] [expr {-1.1*$esu}] \
                $pinchX $pinchY $damage1 $damage2 $beta
        # Plate / I‑section
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

    # --- Concrete core (confined) ---
    # Core patch: from (cover - halfH, cover - halfB) to (halfH - cover, halfB - cover)
    set x1_core [expr {$cover - $halfH}]
    set y1_core [expr {$cover - $halfB}]
    set x2_core [expr {$halfH - $cover}]
    set y2_core [expr {$halfB - $cover}]
    patch rect $coreTag $NUMFIBERS $NUMFIBERS $x1_core $y1_core $x2_core $y2_core

    # --- Concrete cover (unconfined) strips ---
    # Top strip
    patch rect $coverTag $NUMFIBERS 10 [expr {-$halfH}] [expr {$halfB - $cover}] $halfH $halfB
    # Bottom strip
    patch rect $coverTag $NUMFIBERS 10 [expr {-$halfH}] [expr {-$halfB}] $halfH [expr {$cover - $halfB}]
    # Left strip
    patch rect $coverTag $NUMFIBERS 10 [expr {-$halfH}] [expr {$cover - $halfB}] [expr {$cover - $halfH}] [expr {$halfB - $cover}]
    # Right strip
    patch rect $coverTag $NUMFIBERS 10 [expr {$halfH - $cover}] [expr {$cover - $halfB}] $halfH [expr {$halfB - $cover}]

    # ------------------------------------------------------------------
    #  Steel cross‑shaped profile and external plates (steel_layers)
    # ------------------------------------------------------------------
    # List of steel layers: depth, width, center_y, center_x, nfy, nfx
    set steel_layers {
        {10.0 150.0   -125.0   0.0 10 10}
        {240.0 10.0     0.0   0.0 10 10}
        {10.0 150.0    125.0   0.0 10 10}
        {10.0 105.0     0.0  -57.5 10 10}
        {10.0 105.0     0.0   57.5 10 10}
        {150.0 10.0     0.0 -115.0 10 10}
        {150.0 10.0     0.0  115.0 10 10}
        {$PT [expr {$BSec + 2.0*$PT}] [expr {0.5*$HSec + 0.5*$PT}] 0.0 10 10}
        {$PT [expr {$BSec + 2.0*$PT}] [expr {-0.5*$HSec - 0.5*$PT}] 0.0 10 10}
        {$HSec $PT 0.0 [expr {0.5*$BSec + 0.5*$PT}] 10 10}
        {$HSec $PT 0.0 [expr {-0.5*$BSec - 0.5*$PT}] 10 10}
    }

    set steel_area 0.0
    set min_bottom 1e9
    set max_top   -1e9

    foreach layer $steel_layers {
        # Evaluate expressions inside the list (since PT, HSec, BSec are variables)
        lassign [uplevel 1 [list subst $layer]] depth width cy cx nfy nfx
        set x_left [expr {$cx - $width/2.0}]
        set x_right [expr {$cx + $width/2.0}]
        set y_bot  [expr {$cy - $depth/2.0}]
        set y_top  [expr {$cy + $depth/2.0}]

        # Update overall extents (for information only)
        if {$y_bot < $min_bottom} { set min_bottom $y_bot }
        if {$y_top > $max_top}    { set max_top    $y_top }

        # Add steel area
        set steel_area [expr {$steel_area + $depth * $width}]

        # Create patch
        patch rect $steelITag $nfy $nfx $x_left $y_bot $x_right $y_top
    }

    set steel_section_height [expr {$max_top - $min_bottom}]
    puts "Steel Section height = $steel_section_height mm"
    puts "Steel Section area   = $steel_area mm²"

    # ------------------------------------------------------------------
    #  Reinforcing bars (discrete fibers)
    # ------------------------------------------------------------------
    # Rebar list: diameter (mm), y‑coord (mm), x‑coord (mm)
    set rebars {
        [list $REBAR_DIA [expr {$halfH - $cover}] [expr {-$halfB + $cover}]]
        [list $REBAR_DIA [expr {$halfH - $cover}] [expr { $halfB - $cover}]]
        [list $REBAR_DIA [expr {-$halfH + $cover}] [expr {-$halfB + $cover}]]
        [list $REBAR_DIA [expr {-$halfH + $cover}] [expr { $halfB - $cover}]]
        [list 16.0 [expr {$halfH/2.0 - $cover}] [expr {-$halfB + $cover}]]
        [list 16.0 [expr {$halfH/2.0 - $cover}] [expr { $halfB - $cover}]]
        [list 16.0 [expr {-$halfH/2.0 + $cover}] [expr {-$halfB + $cover}]]
        [list 16.0 [expr {-$halfH/2.0 + $cover}] [expr { $halfB - $cover}]]
        [list $REBAR_DIA 0.0 [expr {-$halfB + $cover}]]
        [list $REBAR_DIA 0.0 [expr { $halfB - $cover}]]
        [list $REBAR_DIA [expr {$halfH - $cover}] 0.0]
        [list $REBAR_DIA [expr {-$halfH + $cover}] 0.0]
    }

    foreach rebar $rebars {
        lassign $rebar dia y x
        set area [expr {3.141592653589793 * $dia * $dia / 4.0}]
        fiber $x $y $area $steelTag
    }

    # ------------------------------------------------------------------
    #  Mass per unit length (kg/mm) – following original Python formula
    # ------------------------------------------------------------------
    set concrete_mass [expr {$HSec * $BSec * $CONCRETE_DENSITY}]
    set steel_mass    [expr {$steel_area * $STEEL_DENSITY}]
    set total_mass    [expr {$concrete_mass + $steel_mass}]

    set overall_height [expr {$HSec + 2.0 * $PT}]

    return [list $overall_height $total_mass]
}

#set secTag 1
#set HSec 500.0
#set BSec 500.0
#set cover 40.0
#set steelType "INELASTIC"
#set fc 30.0          ; # positive MPa
#set Kc 1.3
#set steelDensity 7850.0
#set concreteDensity 2400.0
#set plotFlag 0

#set result [COMPOSITE_REC_CROSS_SHAPED_STEEL_SECTION_FUN_EXTRA $secTag $HSec $BSec $cover $steelType $fc $Kc $steelDensity $concreteDensity $plotFlag]
#set totalHeight [lindex $result 0]
#set totalMass   [lindex $result 1]
#puts "Total section height = $totalHeight mm"
#puts "Total mass per unit length = $totalMass kg/mm"