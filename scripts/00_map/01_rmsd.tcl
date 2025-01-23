# This is a script to calculate the RMSD for all frames of a molecule
# Usage: rmsd <molecule id> <seltext>
# 19 to 352 3a
# 22 to 337 8a
proc rmsd { molid1 fitstruct fitwntstart fitwntend seltext1 range outfilename} {
    # jie: no need to load ref structure multiple times
    #mol new input/Wnt8a_align.pdb
    mol new $fitstruct

    set molid2 [ molinfo index 1 ]
    set refwntstart 22
    set refwntend 337


    set separator " "  
    set rms {}
    
#    set currbegin [ expr $seltext1 - $range *2 ]
#    set currend [ expr $seltext1 + $range *2]
#
#    if { $currend > $fitwntend } {
#	set rangestart [expr $refwntstart + 2*$range]
#	set rangeend [expr $refwntend]
#    } else {			 
#	set rangestart [expr $refwntstart] 
#	set rangeend [expr $refwntend - 2*$range]
#    }


    set currbegin [ expr $seltext1 - $range  ]
    set currend [ expr $seltext1 + $range ]
    
    if {$currbegin < $fitwntstart  } { # N-ter, look at next 2*range AA
	set rangestart [expr $refwntstart] 
	set rangeend [expr $refwntend - 2*$range]
    } elseif {$currend > $fitwntend } { # C-ter, look at the previous 2*range AA
	set rangestart [expr $refwntstart + 2*$range]
	set rangeend [expr $refwntend]
    } else { # look at (i-range, i+range) AA
	set rangestart [expr $refwntstart + $range]
	set rangeend [expr $refwntend - $range ]
    }
    

    
    puts "Searching for:${seltext1}\nCurrent Begin:${rangestart}\nCurrent End:${rangeend}"
    for { set i $rangestart } { $i <= $rangeend } { incr i } {

#	if { $currend > $fitwntend } {
#	    set end1 [expr $i]
#	    set start1 [expr $i - 2*$range]
#	    set end2 [expr $seltext1]
#	    set start2 [expr $seltext1 - 2*$range]
#	} else {
#	    set start1 [expr $i]
#	    set end1 [expr $i + 2*$range]
#	    set start2 [expr $seltext1]
#	    set end2 [expr $seltext1 + 2*$range]
#	}

	# jie: re-write the range condition
	if {$currbegin < $fitwntstart  } { # N-ter, look at next 2*range AA
	    set start1 [expr $i]
	    set end1 [expr $i + 2*$range]
	    set start2 [expr $seltext1]
	    set end2 [expr $seltext1 + 2*$range]
	} elseif {$currend > $fitwntend } { # C-ter, look at the previous 2*range AA
	    set end1 [expr $i]
	    set start1 [expr $i - 2*$range]
	    set end2 [expr $seltext1]
	    set start2 [expr $seltext1 - 2*$range]	    
	} else {		# look at (i-range, i+range) AA
	    set end1 [expr $i + $range]
	    set start1 [expr $i - $range]
	    set end2 [expr $seltext1 + $range]
	    set start2 [expr $seltext1 - $range]   
	}
    

	
	puts "Currently ending at: ${currend}\n"
	puts "Searching for:${seltext1}\nCurrent Begin:${rangestart}\nCurrent End:${rangeend}\nSearch Range:${start1} - ${end1}\nCurr Range:${start2} - ${end2}\n"

	# jie: keeping PROB (WLS) here to avoid local fitting between two different secondary structure
	set ref2 [atomselect $molid1 "(resid $start1 to $end1 and name CA and segname PROA) or (segname PROB and name CA and resid 4 to 496)"]
	set sel2 [atomselect $molid2 "(resid $start2 to $end2 and name CA and segname PROA) or (segname PROB and name CA and resid 4 to 496)"]

	set transformation_matrix [measure fit $sel2 $ref2] 
	$sel2 move $transformation_matrix
    #set rmssingle [measure rmsd $sel2 $ref2]
	#puts "Start to end: ${start1} - ${end1}\nStart to end: ${start2} - ${end2}\nRMSD: ${rmssingle}\n"

	# jie: actual rms is calculated excluding PROB
	set ref3 [atomselect $molid1 "(resid $start1 to $end1 and name CA and segname PROA)"]
	set sel3 [atomselect $molid2 "(resid $start2 to $end2 and name CA and segname PROA)"]
	
	lappend rms [measure rmsd $sel3 $ref3 ]
    }
  
    set outfile1 [open $outfilename a]
#    puts -nonewline $outfile1 $rangestart$separator     
    # jie: rangestart is not informative enough
    puts -nonewline $outfile1 $seltext1$separator$rangestart$separator     
    puts $outfile1 $rms
    close $outfile1
    #mol delete $molid1
    mol delete $molid2
}

# 32 to 369 1
# 19 to 352 3a
# 22 to 337 8a
proc rmsd_1 { } {
    # jie only load ref struc once
    mol new input/Wnt8a_align.pdb
    set molid1 [ molinfo index 0 ]

    set wnt1astart 32
    set wnt1aend 369
    for { set i $wnt1astart } { $i <= $wnt1aend } {incr i} {
	puts $i   
	rmsd $molid1  "input/Wnt1_align.pdb" $wnt1astart $wnt1aend $i 20 "output/1_to_8a_jie.out"
    }
    mol delete $molid1
}
proc rmsd_3 { } {
    mol new input/Wnt8a_align.pdb
    set molid1 [ molinfo index 0 ]

    set wnt3astart 19
    set wnt3aend 352
    for { set i $wnt3astart } { $i <= $wnt3aend } {incr i} {
	puts $i
	rmsd $molid1  "input/Wnt3a_align.pdb" $wnt3astart $wnt3aend $i 20 "output/3a_to_8a_jie.out"
    }
    mol delete $molid1
}

proc rmsd_4 { } {
    mol new input/Wnt8a_align.pdb
    set molid1 [ molinfo index 0 ]

    set wnt4astart 24
    set wnt4aend 352
    for { set i $wnt4astart } { $i <= $wnt4aend } {incr i} {
	puts $i
	rmsd $molid1  "input/Wnt4a_align.pdb" $wnt4astart $wnt4aend $i 20 "output/4_to_8a_jie.out"
    }
    mol delete $molid1
}

proc rmsd_5 { } {
    mol new input/Wnt8a_align.pdb
    set molid1 [ molinfo index 0 ]

    set wnt5astart 44
    set wnt5aend 380
    for { set i $wnt5astart } { $i <= $wnt5aend } {incr i} {
	puts $i
	rmsd $molid1  "input/Wnt5a_align.pdb" $wnt5astart $wnt5aend $i 20 "output/5a_to_8a_jie.out"
    }
    mol delete $molid1
}

proc rmsd_8 { } {
    mol new input/Wnt8a_align.pdb
    set molid1 [ molinfo index 0 ]

    set wnt8astart 22
    set wnt8aend 337
    for { set i $wnt8astart } { $i <= $wnt8aend } {incr i} {
	puts $i
	rmsd $molid1  "input/Wnt8a_align.pdb" $wnt8astart $wnt8aend $i 20 "output/8a_to_8a_jie.out"
    }
    mol delete $molid1
}
