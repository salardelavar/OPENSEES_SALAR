# STEEL_CORNER_BOX_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D
# Create a steel corner‑box fiber section (OpenSees) with elastic or hysteretic material.
# Returns section height (mm) and mass per unit length (kg/mm).
#
# Arguments:
#   secTag        - section identifier (integer)
#   STEEL_TYPE    - "ELASTIC" or "INELASTIC"
#   STEEL_DENSITY - density in kg/m³ (used directly: density * area -> kg/mm)
#
# NOTE: The mass calculation follows the original Python code and is physically
#       incorrect (should be multiplied by 1e-9), but kept for compatibility.
# THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
proc SteelCornerBoxSectionDiffWidthHeightFunExtra {secTag STEEL_TYPE STEEL_DENSITY} {
    # ------------------------------------------------------------------
    #  Material properties (steel)
    # ------------------------------------------------------------------
    set fyI   240.0      ; # N/mm²
    set EsI   200000.0   ; # N/mm²
    set eyI   [expr {$fyI / $EsI}]
    set fuI   [expr {1.1818 * $fyI}]
    set esuI  0.25
    set EshI  [expr {($fuI - $fyI) / ($esuI - $eyI)}]
    set BsI   [expr {$EshI / $EsI}]

    set coreTag   [expr {$secTag + 100}]
    set coverTag  [expr {$secTag + 200}]
    set steelTag  [expr {$secTag + 300}]
    set steelITag [expr {$secTag + 400}]

    # ------------------------------------------------------------------
    #  Create uniaxial material
    # ------------------------------------------------------------------
    if {$STEEL_TYPE eq "ELASTIC"} {
        uniaxialMaterial Elastic $steelITag $EsI
    } elseif {$STEEL_TYPE eq "INELASTIC"} {
        set pinchX   0.8
        set pinchY   0.5
        set damage1  0.0
        set damage2  0.0
        set beta     0.1
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
    #  Geometry parameters (mm)
    # ------------------------------------------------------------------
    set CD   300.0    ; # Section depth (overall)
    set ACD  100.0    ; # Corner depth (distance from outer edge to corner plate)
    set PT   10.0     ; # Plate thickness
    set nfy  10       ; # Number of fiber divisions in Y
    set nfx  10       ; # Number of fiber divisions in X

    # Each patch: depth, width, center_y, center_x, nfy, nfx, material tag
    # The 12 patches are defined exactly as in the original Python `mat_layers`.
    set patches {
        {10 110  285  235 10 10}
        {10 110   -5   45 10 10}
        {100 10  250  305 10 10}
        {100 10   50   -5 10 10}
        {10 110  285   45 10 10}
        {10 110   -5  235 10 10}
        {100 10  250   -5 10 10}
        {100 10   50  305 10 10}
        {10 200  315  150 10 10}
        {10 200  -15  150 10 10}
        {200 10  150  315 10 10}
        {200 10  150  -15 10 10}
    }

    set min_bottom  1e9
    set max_top    -1e9
    set total_area 0.0

    foreach patch $patches {
        lassign $patch depth width cy cx nfy nfx

        set x_left [expr {$cx - $width/2.0}]
        set x_right [expr {$cx + $width/2.0}]
        set y_bot  [expr {$cy - $depth/2.0}]
        set y_top  [expr {$cy + $depth/2.0}]

        # Update overall extents
        if {$y_bot < $min_bottom} { set min_bottom $y_bot }
        if {$y_top > $max_top}    { set max_top    $y_top }

        # Add area of this patch
        set total_area [expr {$total_area + $depth * $width}]

        # Create rectangular patch in fiber section
        patch rect $steelITag $nfy $nfx $x_left $y_bot $x_right $y_top
    }

    set section_height [expr {$max_top - $min_bottom}]
    puts "Section height = $section_height mm"
    puts "Section area   = $total_area mm²"

    # ------------------------------------------------------------------
    #  Mass per unit length (kg/mm) – following original Python formula
    # ------------------------------------------------------------------
    set ele_mass [expr {$STEEL_DENSITY * $total_area}]

    return [list $section_height $ele_mass]
}

#set secTag 1
#set steelType "INELASTIC"
#set density 7850.0   ; # kg/m³
#set result [SteelCornerBoxSectionDiffWidthHeightFunExtra $secTag $steelType $density]
#set height [lindex $result 0]
#set mass   [lindex $result 1]
#puts "Height = $height mm, Mass = $mass kg/mm"