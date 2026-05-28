# STEEL_UNP_BOX_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D
# Create a steel fiber section with different widths and heights (3D box section)
#
# Arguments:
#   secTag        - section identifier (integer)
#   STEEL_TYPE    - "ELASTIC" or "INELASTIC"
#   STEEL_DENSITY - density in kg/m³ (converted to N·s²/mm³ implicitly)
#
# Returns:
#   A Tcl list containing: [section_height_mm mass_per_unit_length_kg_per_mm]
#
# NOTE: The mass calculation follows the original Python code: 
#       mass = density (kg/m³) * area (mm²) -> kg/mm. This is physically incorrect
#       (should be multiplied by 1e-9), but kept for exact compatibility.
# THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
proc SteelUnpBoxSectionDiffWidthHeightFunExtra3D {secTag STEEL_TYPE STEEL_DENSITY} {
    # ------------------------------------------------------------------
    #  Material properties (steel I-section)
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
    set CD   280.0    ; # Section Depth
    set ACD  100.0    ; # Corner Depth
    set PT   10.0     ; # Plate Thickness
    set nfy  10       ; # Number of fiber divisions in Y
    set nfx  10       ; # Number of fiber divisions in X

    # Each patch: depth, width, center_y, center_x, nfy, nfx, matTag
    set patches {
        {10 110  285  235 10 10}
        {10 110   -5   45 10 10}
        {280 10  140   -5 10 10}
        {10 110  285   45 10 10}
        {10 110   -5  235 10 10}
        {280 10  140  285 10 10}
        {10 180  295  140 10 10}
        {10 180  -15  140 10 10}
    }

    set min_bottom 1e9
    set max_top   -1e9
    set total_area 0.0

    foreach patch $patches {
        lassign $patch depth width cy cx nfy nfx
        set x_left [expr {$cx - $width/2.0}]
        set x_right [expr {$cx + $width/2.0}]
        set y_bot  [expr {$cy - $depth/2.0}]
        set y_top  [expr {$cy + $depth/2.0}]

        # Update extents
        if {$y_bot < $min_bottom} { set min_bottom $y_bot }
        if {$y_top > $max_top}    { set max_top    $y_top }

        # Add area
        set total_area [expr {$total_area + $depth * $width}]

        # Create rectangular patch
        patch rect $steelITag $nfy $nfx $x_left $y_bot $x_right $y_top
    }

    set section_height [expr {$max_top - $min_bottom}]
    puts "Section height = $section_height mm"
    puts "Section area   = $total_area mm²"

    # ------------------------------------------------------------------
    #  Mass per unit length (kg/mm) – follows original Python calculation
    #  (no conversion from kg/m³ to kg/mm³)
    # ------------------------------------------------------------------
    set ele_mass [expr {$STEEL_DENSITY * $total_area}]   ; # kg/mm (physically incorrect but compatible)

    return [list $section_height $ele_mass]
}

# Create a steel section
#set secTag 1
#set steelType "INELASTIC"
#set density 7850.0   ; # kg/m³
#set result [SteelUnpBoxSectionDiffWidthHeightFunExtra $secTag $steelType $density]
#set height [lindex $result 0]
#set mass   [lindex $result 1]
#puts "Height = $height mm, Mass = $mass kg/mm"
