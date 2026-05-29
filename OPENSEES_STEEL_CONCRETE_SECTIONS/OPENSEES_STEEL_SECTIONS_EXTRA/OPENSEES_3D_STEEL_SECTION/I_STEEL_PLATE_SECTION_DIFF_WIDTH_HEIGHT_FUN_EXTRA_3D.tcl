# ----------------------------------------------------------------------
# I_STEEL_PLATE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D
#
# Creates a steel fiber section with varying widths/heights.
# Computes section height and mass per unit length.
#
# Arguments:
#   secTag        - integer, section identifier
#   STEEL_TYPE    - string, "ELASTIC" or "INELASTIC"
#   STEEL_DENSITY - float, density in kg/m³ (converted internally to kg/mm³)
#   plot          - integer (0/1), plotting not supported in Tcl, ignored
#
# Returns:
#   A list containing two values: SECTION_HEIGHT (mm), ELE_MASS (kg/mm)
#
# THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
# ----------------------------------------------------------------------
proc I_STEEL_PLATE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D {secTag STEEL_TYPE STEEL_DENSITY {plot 0}} {
    # ------------------------------------------------------------------
    # 1. Material parameters (units: N, mm, s)
    # ------------------------------------------------------------------
    set fy   240.0      ;# Yield strength [N/mm²]
    set Es   200000.0   ;# Modulus of elasticity [N/mm²]
    set ey   [expr {$fy / $Es}]
    set fu   [expr {1.1818 * $fy}]
    set esu  0.15
    set Esh  [expr {($fu - $fy) / ($esu - $ey)}]
    set Bs   [expr {$Esh / $Es}]

    set steelTag [expr {$secTag + 100}]

    # ------------------------------------------------------------------
    # 2. Define uniaxial material (elastic or hysteretic)
    # ------------------------------------------------------------------
    if {[string toupper $STEEL_TYPE] eq "ELASTIC"} {
        uniaxialMaterial Elastic $steelTag $Es
    } elseif {[string toupper $STEEL_TYPE] eq "INELASTIC"} {
        set pinchX  0.8
        set pinchY  0.5
        set damage1 0.0
        set damage2 0.0
        set beta    0.1

        uniaxialMaterial Hysteretic $steelTag \
            $fy  $ey  $fu  $esu  [expr {0.2*$fu}] [expr {1.1*$esu}] \
            [expr {-$fy}] [expr {-$ey}] [expr {-$fu}] [expr {-$esu}] \
            [expr {-0.2*$fu}] [expr {-1.1*$esu}] \
            $pinchX $pinchY $damage1 $damage2 $beta
    } else {
        puts "ERROR: STEEL_TYPE must be 'ELASTIC' or 'INELASTIC'"
        return -code error
    }

    # ------------------------------------------------------------------
    # 3. Define fiber section and patches
    # ------------------------------------------------------------------
    section Fiber $secTag -GJ 1.0e7

    # Each sublist: depth, width, center_y, center_x, numSubdivY, numSubdivZ
    set steel_layers {
        {10.0 100.0   5.0 100.0  10 10}
        {160.0 10.0  90.0 100.0  10 10}
        {10.0 100.0 175.0 100.0  10 10}
        {10.0  75.0  -5.0 100.0  10 10}
        {10.0  75.0 185.0 100.0  10 10}
        {100.0 10.0  90.0  90.0  10 10}
        {100.0 10.0  90.0 110.0  10 10}
    }

    # Initialize min/max and total area
    set min_bottom 1e9
    set max_top   -1e9
    set total_area 0.0

    foreach layer $steel_layers {
        lassign $layer depth width center_y center_x numY numZ

        set x_left [expr {$center_x - $width/2.0}]
        set x_right [expr {$center_x + $width/2.0}]
        set y_bot [expr {$center_y - $depth/2.0}]
        set y_top [expr {$center_y + $depth/2.0}]

        # Add patch
        patch rect $steelTag $numY $numZ $x_left $y_bot $x_right $y_top

        # Update min/max and area
        if {$y_bot < $min_bottom} {set min_bottom $y_bot}
        if {$y_top > $max_top}    {set max_top    $y_top}
        set total_area [expr {$total_area + $depth * $width}]
    }

    set section_height [expr {$max_top - $min_bottom}]

    puts "Section height = $section_height mm"
    puts "Section area   = $total_area mm²"

    # ------------------------------------------------------------------
    # 4. Compute mass per unit length (kg/mm)
    #    Convert density from kg/m³ to kg/mm³ (1 kg/m³ = 1e-9 kg/mm³)
    # ------------------------------------------------------------------
    set density_kg_per_mm3 [expr {$STEEL_DENSITY * 1.0e-9}]
    set ele_mass [expr {$density_kg_per_mm3 * $total_area}]

    # Plotting is not supported in native Tcl/OpenSees – skip or output info
    if {$plot} {
        puts "Plotting not available in Tcl. Use external post‑processing if needed."
    }

    return [list $section_height $ele_mass]
}


#set result [I_STEEL_PLATE_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D 1 "INELASTIC" 7850.0 0]
#lassign $result height massPerUnitLength
#puts "Height = $height mm, Mass per length = $massPerUnitLength kg/mm"