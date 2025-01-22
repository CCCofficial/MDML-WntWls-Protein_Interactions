proc change_transparency {new_alpha} {
        # This will always get the correct colors even if VMD
        # gains new definitions in the future
        set color_start [colorinfo num]
        set color_end [expr $color_start * 2]
        # Go through the list of colors (by index) and 
        # change their transp. value
        for {set color $color_start} {$color < $color_end} {incr color} {
                color change alpha $color $new_alpha
        }
}

proc draw_contacts {bonds colsel} {
set myfile $bonds
set a [open $myfile]
set lines [split [read -nonewline $a] "\n"]
close $a;                          # Saves a few bytes :-)

set countersWNT {}
set counterWNTLESS {}
foreach line $lines {
    set respair [ split $line "_" ]
    set l1 [ lindex $respair 0 ]
    set l2 [ lindex $respair 1 ]
    dict incr countersWNT $l1
    dict incr countersWNTLESS $l2
}

dict for {item count} $countersWNT {
    set sel3 [atomselect top "resid $item and name CA and segid PROA"]
    set atom1 [lindex [$sel3 get { x y z }] 0]
    draw color $colsel
    #graphics top sphere $atom1 radius [ expr $count/12 ]
    graphics top sphere $atom1 radius 0.5
}
dict for {item count} $countersWNTLESS {
    set sel3 [atomselect top "resid $item and name CA and segid PROB"]
    set atom1 [lindex [$sel3 get { x y z }] 0]
    draw color $colsel
    #graphics top sphere $atom1 radius [ expr $count/12 ]
    graphics top sphere $atom1 radius 0.5
}
foreach line $lines {
    set respair [ split $line "_" ]
    set l1 [ lindex $respair 0 ]
    set l2 [lindex $respair 1 ]	
    set sel3 [atomselect top "resid $l1 and name CA and segid PROA"]
    set atom1 [lindex [$sel3 get { x y z }] 0]
    set sel4 [atomselect top "resid $l2 and name CA and segid PROB"]
    set atom2 [lindex [$sel4 get { x y z }] 0]
    set idx1 [$sel3 get index]
    set idx2 [$sel4 get index]
    draw color $colsel
    graphics top cylinder $atom1 $atom2 radius 0.25
}
}

proc run {} {
mol new ../00_map/input/Wnt8a_align.pdb
draw_contacts ../07_preprocess/output/region_red.txt red
draw_contacts ../07_preprocess/output/region_green.txt green
draw_contacts ../07_preprocess/output/region_blue.txt blue
draw_contacts ../07_preprocess/output/region_purple.txt purple
draw_contacts ../07_preprocess/output/region_gray.txt gray
}
