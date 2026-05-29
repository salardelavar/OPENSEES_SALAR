###############################################################################
# COMPOSITE_BEAM_REC_I_SECTION_FUN_EXTRA_3D
# Creates a composite beam‑slab fiber section in OpenSees (Tcl version).
#
# Arguments:
#   secTag            - section identifier (integer)
#   Hsec, Bsec        - slab height and width [mm]
#   STEEL_TYPE        - "ELASTIC" or "INELASTIC"
#   nFib              - number of fiber divisions per direction (for patches)
#   fc                - unconfined concrete compressive strength [MPa] (negative)
#   Kfc               - ratio confined/unconfined concrete strength
#   bf, tf            - I‑section flange width and thickness [mm]
#   h, tw             - I‑section total height and web thickness [mm]
#   STEEL_DENSITY     - steel density [kg/m³]  (used directly, see note)
#   CONCRETE_DENSITY  - concrete density [kg/m³] (used directly, see note)
#   PLOT              - ignored in Tcl (no plotting)
#
# Returns:
#   A Tcl list: {total_height mass_per_unit_length}
#   total_height = Hsec + h  [mm]
#   mass_per_unit_length = Hsec*Bsec*CONCRETE_DENSITY + AREA*STEEL_DENSITY
#   (Note: Units are inconsistent in original Python; this replicates that.)
#
# THIS TCL SCRIPT WRITTEN BY SALAR DELAVAR GHASHGHAEI (QASHQAI)
###############################################################################
proc COMPOSITE_BEAM_REC_I_SECTION_FUN_EXTRA_3D {secTag Hsec Bsec STEEL_TYPE nFib fc Kfc bf tf h tw STEEL_DENSITY CONCRETE_DENSITY PLOT} {

    # -------------------- Material tags --------------------
    set coreTag   [expr {$secTag + 100}]
    set coverTag  [expr {$secTag + 200}]
    set steelTag  [expr {$secTag + 300}]
    set steelITag [expr {$secTag + 400}]

    # ----- Concrete material parameters (unconfined / cover) -----
    set fcU   [expr {-$fc}]                      # positive compressive strength [N/mm²]
    set EcU   [expr {4700 * sqrt($fcU)}]         # [N/mm²]
    set ec0U  [expr {2.0 * $fcU / $EcU}]         # strain at peak stress
    set fcUU  [expr {0.2 * $fcU}]                # ultimate compressive strength
    set ecuU  [expr {5.0 * $ec0U}]               # ultimate strain
    set lambdaU 0.1
    set ftU   [expr {0.7 * sqrt($fcU)}]          # tensile strength [N/mm²]
    set EtsU  [expr {$ftU / abs($ec0U)}]         # tension softening stiffness

    # ----- Concrete material parameters (confined / core) -----
    set fcC   [expr {$Kfc * $fcU}]               # confined compressive strength
    set EcC   [expr {4700 * sqrt($fcC)}]
    set ec0C  [expr {2.0 * $fcC / $EcC}]
    set fcUC  [expr {0.65 * $fcC}]               # ultimate compressive strength (confined)
    set ecuC  [expr {15.0 * $ec0C}]              # ultimate strain (confined)
    set lambdaC 0.1
    set ftC   [expr {0.7 * sqrt($fcC)}]
    set EtsC  [expr {$ftC / abs($ec0C)}]

    # ----- Steel material parameters (rebar) -----
    set fy    400.0      # [N/mm²]
    set Es    200000.0   # [N/mm²]
    set ey    [expr {$fy / $Es}]
    set fu    [expr {1.1818 * $fy}]
    set esu   0.09
    set Esh   [expr {($fu - $fy) / ($esu - $ey)}]
    set Bs    [expr {$Esh / $Es}]

    # ----- Steel I‑section parameters -----
    set fyI   240.0
    set EsI   200000.0
    set eyI   [expr {$fyI / $EsI}]
    set fuI   [expr {1.1818 * $fyI}]
    set esuI  0.25
    set EshI  [expr {($fuI - $fyI) / ($esuI - $eyI)}]
    set BsI   [expr {$EshI / $EsI}]

    # ----- Create materials -----
    if {$STEEL_TYPE == "ELASTIC"} {
        uniaxialMaterial Elastic $steelTag   $Es
        uniaxialMaterial Elastic $steelITag  $EsI
    } elseif {$STEEL_TYPE == "INELASTIC"} {
        set pinchX  0.8
        set pinchY  0.5
        set damage1 0.0
        set damage2 0.0
        set beta    0.1
        # Rebar
        uniaxialMaterial Hysteretic $steelTag \
            $fy $ey $fu $esu [expr {0.2*$fu}] [expr {1.1*$esu}] \
            [expr {-$fy}] [expr {-$ey}] [expr {-$fu}] [expr {-$esu}] \
            [expr {-0.2*$fu}] [expr {-1.1*$esu}] \
            $pinchX $pinchY $damage1 $damage2 $beta
        # I‑section
        uniaxialMaterial Hysteretic $steelITag \
            $fyI $eyI $fuI $esuI [expr {0.2*$fuI}] [expr {1.1*$esuI}] \
            [expr {-$fyI}] [expr {-$eyI}] [expr {-$fuI}] [expr {-$esuI}] \
            [expr {-0.2*$fuI}] [expr {-1.1*$esuI}] \
            $pinchX $pinchY $damage1 $damage2 $beta
    } else {
        puts "WARNING: STEEL_TYPE must be ELASTIC or INELASTIC"
    }

    # Concrete02 (includes tension)
    uniaxialMaterial Concrete02 $coreTag  $fcC $ec0C $fcUC $ecuC $lambdaC $ftC $EtsC
    uniaxialMaterial Concrete02 $coverTag $fcU $ec0U $fcUU $ecuU $lambdaU $ftU $EtsU

    # -------------------- Fiber section --------------------
    section Fiber $secTag -GJ 1.0e7

    # ----- Slab (concrete core) -----
    # patch rect matTag nY nX yI zI yJ zJ
    # Note: y is vertical, x is horizontal. Here slab runs from y=0 to y=Hsec.
    patch rect $coreTag $nFib $nFib [expr {-$Bsec/2.0}] 0.0 [expr {$Bsec/2.0}] $Hsec

    # ----- I‑section (steel) -----
    # Top flange (y = -tf to 0)
    patch rect $steelITag $nFib $nFib [expr {-$bf/2.0}] [expr {-$tf}] [expr {$bf/2.0}] 0.0
    # Bottom flange (y = -h to -h+tf)
    patch rect $steelITag $nFib $nFib [expr {-$bf/2.0}] [expr {-$h}] [expr {$bf/2.0}] [expr {-$h+$tf}]
    # Web (y = -h+tf to -tf)
    patch rect $steelITag $nFib $nFib [expr {-$tw/2.0}] [expr {-$h+$tf}] [expr {$tw/2.0}] [expr {-$tf}]

    # ----- Compute I‑section area (for mass) -----
    set area_top_flange    [expr {$bf * $tf}]
    set area_bottom_flange [expr {$bf * $tf}]
    set area_web           [expr {$tw * ($h - 2.0*$tf)}]
    set AREA_I             [expr {$area_top_flange + $area_bottom_flange + $area_web}]
    puts "Area of Top Flange: $area_top_flange mm^2"
    puts "Area of Bottom Flange: $area_bottom_flange mm^2"
    puts "Area of Web: $area_web mm^2"
    puts "Total I Section Area: $AREA_I mm^2"

    # ----- Reinforcing bars (rebars) -----
    set RD   16.0    # [mm] bar diameter
    set DIST 150.0   # [mm] spacing offset

    # Coordinates (y, x) for 18 bars as in the original code
    set rebar_coords {
        { [expr {$Hsec/4.0}] [expr {-$Bsec/2.0 + 1*$DIST}] }
        { [expr {$Hsec/4.0}] [expr { $Bsec/2.0 - 1*$DIST}] }
        { [expr {$Hsec/4.0}] [expr {-$Bsec/2.0 + 2*$DIST}] }
        { [expr {$Hsec/4.0}] [expr { $Bsec/2.0 - 2*$DIST}] }
        { [expr {$Hsec/4.0}] [expr {-$Bsec/2.0 + 3*$DIST}] }
        { [expr {$Hsec/4.0}] [expr { $Bsec/2.0 - 3*$DIST}] }
        { [expr {$Hsec/4.0}] [expr {-$Bsec/2.0 + 4*$DIST}] }
        { [expr {$Hsec/4.0}] [expr { $Bsec/2.0 - 4*$DIST}] }
        { [expr {$Hsec/4.0}] 0.0 }
        { [expr {3.0*$Hsec/4.0}] [expr {-$Bsec/2.0 + 1*$DIST}] }
        { [expr {3.0*$Hsec/4.0}] [expr { $Bsec/2.0 - 1*$DIST}] }
        { [expr {3.0*$Hsec/4.0}] [expr {-$Bsec/2.0 + 2*$DIST}] }
        { [expr {3.0*$Hsec/4.0}] [expr { $Bsec/2.0 - 2*$DIST}] }
        { [expr {3.0*$Hsec/4.0}] [expr {-$Bsec/2.0 + 3*$DIST}] }
        { [expr {3.0*$Hsec/4.0}] [expr { $Bsec/2.0 - 3*$DIST}] }
        { [expr {3.0*$Hsec/4.0}] [expr {-$Bsec/2.0 + 4*$DIST}] }
        { [expr {3.0*$Hsec/4.0}] [expr { $Bsec/2.0 - 4*$DIST}] }
        { [expr {3.0*$Hsec/4.0}] 0.0 }
    }

    set bar_area [expr {3.141592653589793 * ($RD/2.0)*($RD/2.0)}]   # mm²
    foreach coord $rebar_coords {
        lassign $coord yCoord xCoord
        fiber $xCoord $yCoord $bar_area $steelTag
    }

    # ----- Mass per unit length (consistent with original Python) -----
    set MASS [expr {$Hsec * $Bsec * $CONCRETE_DENSITY + $AREA_I * $STEEL_DENSITY}]
    set total_height [expr {$Hsec + $h}]

    # ----- Optional plotting (not implemented in Tcl) -----
    if {$PLOT} {
        puts "Plotting is not supported in the Tcl version. Use Python for visualization."
    }

    return [list $total_height $MASS]
}

#set result [COMPOSITE_BEAM_REC_I_SECTION_FUN_EXTRA_3D 1 200 1500 "INELASTIC" 10 -30 1.3 200 20 400 12 7850 2400 0]
#set totalH [lindex $result 0]
#set mass   [lindex $result 1]
