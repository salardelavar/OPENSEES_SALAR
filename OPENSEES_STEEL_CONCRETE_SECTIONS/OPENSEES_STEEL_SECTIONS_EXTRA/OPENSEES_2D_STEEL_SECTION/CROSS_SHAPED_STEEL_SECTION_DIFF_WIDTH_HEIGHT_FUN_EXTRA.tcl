# CROSS_SHAPED_STEEL_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA
# Create a steel fiber section (OpenSees) with 7 rectangular patches.
# Returns section height (mm) and mass per unit length (kg/mm).
#
# Arguments:
#   secTag        - section identifier (integer)
#   STEEL_TYPE    - "ELASTIC" or "INELASTIC"
#   STEEL_DENSITY - density in kg/m³ (used directly: density * area -> kg/mm)
#
# NOTE: For STEEL_TYPE = "ELASTIC", the original Python uses Steel01
#       (bilinear with hardening). This is replicated exactly.
# THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)

proc CROSS_SHAPED_STEEL_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA {secTag STEEL_TYPE STEEL_DENSITY} {
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
        # Steel01: bilinear steel material with hardening ratio Bs
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
    #  Fiber section definition
    # ------------------------------------------------------------------
    section Fiber $secTag

    # ------------------------------------------------------------------
    #  Steel layers: depth, width, center_y, center_x, nfy, nfx
    #  (exactly as in the original Python `steel_layers`)
    # ------------------------------------------------------------------
    set steel_layers {
        {10.0 100.0   5.0 100.0 10 10}
        {160.0 10.0  90.0 100.0 10 10}
        {10.0 100.0 175.0 100.0 10 10}
        {10.0  75.0  90.0  57.5 10 10}
        {10.0  75.0  90.0 142.5 10 10}
        {100.0 10.0  90.0  15.0 10 10}
        {100.0 10.0  90.0 185.0 10 10}
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
#set result [CROSS_SHAPED_STEEL_SECTION_DIFF_WIDTH_HEIGHT_FUN_EXTRA $secTag $steelType $density]
#set height [lindex $result 0]
#set mass   [lindex $result 1]
#puts "Height = $height mm, Mass = $mass kg/mm"