# STEEL_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA_3D
# Create a steel fiber section (OpenSees) with given geometry using elastic or hysteretic material.
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
proc SteelSectionDiffWidthHeightFunExtra3D {secTag STEEL_TYPE STEEL_DENSITY} {
    # ------------------------------------------------------------------
    #  Steel material properties
    # ------------------------------------------------------------------
    set fy   240.0      ; # N/mm²
    set Es   200000.0   ; # N/mm²
    set ey   [expr {$fy / $Es}]
    set fu   [expr {1.1818 * $fy}]
    set esu  0.15
    set Esh  [expr {($fu - $fy) / ($esu - $ey)}]
    set Bs   [expr {$Esh / $Es}]

    set steelTag [expr {$secTag + 100}]

    # ------------------------------------------------------------------
    #  Create uniaxial material
    # ------------------------------------------------------------------
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
    #  Steel layers: depth, width, center_y, center_x, nfy, nfx
    #  (exactly as in the original Python `steel_layers`)
    # ------------------------------------------------------------------
    set steel_layers {
        {10.0 100.0   5.0 100.0 10 10}
        {160.0 10.0  90.0 100.0 10 10}
        {10.0 200.0 175.0 100.0 10 10}
        {6.0  70.0  -3.0 100.0 10 10}
        {6.0 150.0 183.0 100.0 10 10}
        {120.0 8.0   90.0  91.0 10 10}
        {120.0 8.0   90.0 109.0 10 10}
        {80.0  8.0   10.0  46.0 10 10}
        {80.0  8.0   10.0 154.0 10 10}
        {80.0  8.0  166.0   4.0 10 10}
        {80.0  8.0  166.0 196.0 10 10}
    }

    set min_bottom  1e9
    set max_top    -1e9
    set total_area 0.0

    foreach layer $steel_layers {
        lassign $layer depth width cy cx nfy nfx

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
        patch rect $steelTag $nfy $nfx $x_left $y_bot $x_right $y_top
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
#set result [SteelSectionDiffWidthHeightFunExtra3D $secTag $steelType $density]
#set height [lindex $result 0]
#set mass   [lindex $result 1]
#puts "Height = $height mm, Mass = $mass kg/mm"