#;- Perl -*-
#
########   Alberto Roldan
#
#
# 	Adsorption reaction coordinate = Z, then to obtain the sticky coeficient on Adsorption process
#   		without a determined transition state, free molecule frequencies on X and Y will be used.
#
#  	Quantum effect according to the harmonic Wigner correction
#      		[Henkelman, G.; Arnaldsson, A.; Jónsson, H., Theoretical Calculations of CH4 and H2 Associative Desorption from Ni(111):
#       	Could Subsurface Hydrogen Play an Important Role? J. Chem. Phys. 2006, 124, 044706.]
#
#
#   version 0.0 - from Mathematica_v4.2
#   version 2.0 - SSURFACE is obsolete, now ISITES
#               - ACAT only for surfaces and COVERAGE for intermediates
#   version 3.x - V active [EEXCH] = number of exchanged electrons   ---   E=E0+n_e * KbT*ln[Q_p/Q_r]
#		- Degree of Rate and Selectivity control (DRC & DSC)
#   version 4.x - It reads all the needed information from file.mk.in so it is DFT-Software independent
#
#
#
#
#
eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;


use Math::Trig;
 use Cwd;
 use Term::ANSIColor qw( colored );
 use List::Util qw(min max);
 use Scalar::Util qw(looks_like_number);

#------------------------------------------------------------------------------------------------------------------------------
       $numargs=$#ARGV+1;
   if ($numargs lt 1 ) { print "\n ---  Please, the input file has to be after the executable --- \n \n" ; &example_sub();  exit 0 ; };
   if (!-s $ARGV[0]) { print "\n ---  Please, provide a proper file --- \n \n" ; exit 0 ; };
#------------------------------------------------------------------------------------------------------------------------------  PATHs
# Alberto 05/2020

#   system("pwd |cut -b23- >>pathroot2");
#   open IN, "./pathroot2"; while (<IN>) { $tmpath.=$_;}; close IN; system("rm pathroot2"); @tmpath=split(/\n/,$tmpath);
#       $rootpath2="C:/Users/user/Desktop/ALBERTO/$tmpath[0]";
       $rootpath=getcwd;
#--------------------------------------------------------------------------------------------------------------------------------- open input file
  $filename=$ARGV[0];
print "\nReading file $ARGV[0]\n";
open IN, $ARGV[0]; while (<IN>) {$input.= $_ ;}; close IN ;
     @input=split(/\n/,$input);
     @in=grep(!/^#/,@input);    # weed out comments
#=================================================================================================================================
#
#                                                          READING INPUT
#
#=================================================================================================================================
print "Reading thermodynamics";
    foreach $line (@in) {
        if ($line) {
            @tag=split(/\s+/,$line); foreach $t (@tag) { if ($t) { push(@Nline,$t); }; }; @tag=@Nline; @Nline=();
          if (@tag[0] eq 'RTEMP') { $itemp=@tag[2];
		 	if ($tag[3]) { $ftemp=@tag[3]; };
			if ($tag[4]) { $ttemp=@tag[4]; }; };
          if (@tag[0] eq 'RTIME') { if ($tag[3]) { $itime=0.0; $ftime=@tag[2]; $ttime=@tag[3]; 
				    }else{ $itime=@tag[2]; }; };	
          if (@tag[0] eq 'RVEXT') {    $nVext=@tag[2];
		 	if ($tag[3]) { $pVext=@tag[3]; };
		 	if ($tag[4]) { $sVext=@tag[4]; }; };
          if (@tag[0] eq 'RpH') {      $ipH=@tag[2]; 
		  	if ($tag[3]) { $fpH=@tag[3]; };
			if ($tag[4]) { $spH=@tag[4]; }; };
        }; #--> line
    }; #-->foreach
#==================================================================================================================================
#
#                                                       READING SYSTEMS
#
      $nsystems=0; $nend=0; $nmolecules=0; $countersystems=0; $nprocess=0; 
#==================================================================================================================================
  foreach $line (@in) {  
    if ( $line ) { @tag=split(/\s+/,$line);
      foreach $t (@tag) { if (($t) or (looks_like_number($t))) { push(@Nline,$t); }; }; @tag=@Nline; @Nline=();
            @eqtmp=split(/= /,$line); foreach $e (@eqtmp) { if ($t) { push(@Nline,$t); }; }; @eqtmp=@Nline; @Nline=();
            if (@tag[0] eq 'SYSTEM') { @nsys[$nsystems]=@tag[2]; $countersystems++; };
            if (@tag[0] eq 'XYZPATH') { @ixyz{@nsys[$nsystems]}=@tag[2]; };
            if (@tag[0] eq 'E0') { $e0{@nsys[$nsystems]}=@tag[2]; };
            if (@tag[0] eq 'DEGENERATION') { @degeneration{@nsys[$nsystems]}=@tag[2]; };
            if (@tag[0] eq 'FREQ') { $t=2; @ifreq=''; @imafreq=''; while ((@tag[$t] ne '#') and (@tag[$t])) {
                if (@tag[$t] > 0) { push(@ifreq,@tag[$t]) }else { push(@imafreq,@tag[$t]); }; $t++; };
					@freq = sort {$b <=> $a} @ifreq; @freq{$nsys[$nsystems]}="@freq";
                    @imafr{$nsys[$nsystems]}="@imafreq"; };
            if (@tag[0] eq 'FREQ2D') { $t=2; @ifreq=''; while (($tag[$t] ne '#') and ($tag[$t])) { push(@ifreq,$tag[$t]); $t++; };
                                        @freq2D = sort {$b <=> $a} @ifreq;
                                        @freq2D{$nsys[$nsystems]}="@freq";
					push(@molecules,@nsys[$nsystems]); };


#---------------------------------------------------------------------------------------------------------------------------------molecule
            if (@molecules[$#molecules] eq @nsys[$nsystems]) {
               if (@tag[0] eq 'IMASS') { $t=2; $Av=6.022139922973909E+023; @iimass='';
                  while ((@tag[$t] ne '#') and (@tag[$t])) { push(@iimass,@tag[$t]/($Av*1000)); $t++; };
                                                @imass{@nsys[$nsystems]}="@iimass"; };
               if (@tag[0] eq 'INATOMS') { $t=2; @iiatoms=''; while ((@tag[$t] ne '#') and (@tag[$t])) {
                                         push(@iiatoms,@tag[$t]); $t++; }; @inatoms{@nsys[$nsystems]}="@iiatoms";
                                         $natom=0; ($natom+=$_) for @iiatoms; @Tnatoms{@nsys[$nsystems]}=$natom;
                                         $iTmass=0;
                                        if ($#iimass ne $#iiatoms) {
                                        print "\n-- The number of atomic masses and number of atoms for @nsys[$nsystems] disagree --\n\n";
                                                exit 0 ;
                                        }else{  for ($i=0; $i<=$#iimass; $i++) {
                                                $iTmass+=$iiatoms[$i]*$iimass[$i]; };};
                                        @Tmass{@nsys[$nsystems]}=$iTmass; };
               if (@tag[0] eq 'SYMFACTOR') { @sym{@nsys[$nsystems]}=@tag[2]; };
               if (@tag[0] eq 'INERTIA') { $t=2; @iIM=(); while (@tag[$t] ne '#') { push(@iIM,@tag[$t]); $t++; };
                                         if (($iIM[0] == $iIM[1]) or ($iIM[1] == $iIM[2])) { @linear{@nsys[$nsystems]}="yes";  
					 }else{ @linear{@nsys[$nsystems]}="no"; };
					 if (@linear{@nsys[$nsystems]} eq "yes") {  if ($iIM[1] >= $iIM[0]) { @inertia{@nsys[$nsystems]}=$iIM[1]; 
						 			   }else{ @inertia{@nsys[$nsystems]}=$iIM[0]; };
					 }else{ @inertia{@nsys[$nsystems]}=$iIM[0]*$iIM[1]*$iIM[2]; };


# Alberto
#print "system=@nsys[$nsystems]\t|||linear=@linear{@nsys[$nsystems]}\t||inertia=@inertia{@nsys[$nsystems]}\t||||\n";



 };

               if (@tag[0] eq 'ISITES') { @tmpsites=split(/=/,$line); @tmp=split(/\s+/,@tmpsites[1]);
                  foreach $t (@tmp) { if (($t) or (looks_like_number($t))) { push(@Nline,$t); }; }; @tmp=@Nline; @Nline=();
                  if ( @tmp[0] =~ /^[0-9]+$/) { @nsitetype{@nsys[$nsystems]}=@tmp[0];
                                                @sitetype{@nsys[$nsystems]}=@tmp[1];
                                         }else{ @sitetype{@nsys[$nsystems]}=@tmp[0];
                                                @nsitetype{@nsys[$nsystems]}=@tmp[1]; };
                  if (!$nsitetype{@nsys[$nsystems]}) { @nsitetype{@nsys[$nsystems]}=1; }; };
               if (@tag[0] eq 'IPRESSURE') { @pressure{@nsys[$nsystems]}=@tag[2]; };
               if (@tag[0] eq 'RPRESSURE') { $ipressure{@nsys[$nsystems]}=@tag[2];
                                             $fpressure{@nsys[$nsystems]}=@tag[3];
                                             $spressure{@nsys[$nsystems]}=@tag[4]; };
            }; #-->if molecules
#---------------------------------------------------------------------------------------------------------------------------------catalyst
            if (@molecules[$#molecules] ne @nsys[$nsystems]) { 
               if (@tag[0] eq 'ISITES') { @tmpsites=split(/=/,$line); @tmp=split(/\s+/,@tmpsites[1]);
		  foreach $t (@tmp) { if (($t) or (looks_like_number($t))) { push(@Nline,$t); }; }; @tmp=@Nline; @Nline=();     
                  if ( @tmp[0] =~ /^[0-9]+$/) { @nsitetype{@nsys[$nsystems]}=@tmp[0];
		                                @sitetype{@nsys[$nsystems]}=@tmp[1];
            		                 }else{ @sitetype{@nsys[$nsystems]}=@tmp[0];
					        @nsitetype{@nsys[$nsystems]}=@tmp[1]; };
	          if (!$nsitetype{@nsys[$nsystems]}) { @nsitetype{@nsys[$nsystems]}=1; }; };
	       if (@tag[0] eq 'IACAT') { @Acat{@sitetype{@nsys[$nsystems]}}=@tag[2]; push(@surfaces,@nsys[$nsystems]); };
               if (@tag[0] eq 'RCOVERAGE') { $icov{@nsys[$nsystems]}=@tag[2];
                                             $fcov{@nsys[$nsystems]}=@tag[3];
                                             $scov{@nsys[$nsystems]}=@tag[4]; };
               if (@tag[0] eq 'ICOVERAGE') { $coverage{@nsys[$nsystems]}=@tag[2]; };
            }; #-->if no molecule
#---------------------------------------------------------------------------------------------------------------------------------end system
           if ((@tag[0] eq 'end') or (@tag[0] eq 'END')) { $nend++; $nsystems++; };
     }; #--> line
  }; #-->foreach line
       $nsystems=$countersystems;

#--------------------------------------------------------------------------------------------------------------------------------- assigns Mol, Sur and Cat
    foreach $s (@nsys) { $new="yes";
 	   foreach $mol (@molecules) { if ($s eq $mol) { $new="no"; }; };
 	   if ($new eq "yes") { foreach $sur (@surfaces) { if ($s eq $sur) { $new="no"; }; }; };
 	   if ($new eq "yes") { push(@catalysts,$s); };
    };

#================================================================================================================================= 
#
#                                               SYSTEMS PROPERTIES
#
#=================================================================================================================================




#--------------------------------------------------------------------------------------------------------------------------------- fails
   if ($nend != $nsystems) {  print "\n -  Please, check the number of SYSTEMs and ENDs - \n" ; exit 0; };
print"\n";
   foreach $mol (@molecules) {
      if (!$e0{$mol}) { print "  $mol energy\t... fail\n"; exit 0;
#		 }else{ print "  $mol energy ($e0{$mol})\t... OK\n";
	};   
      if (!@freq{$mol}) { print "  $mol frequencies (@nf{$mol})\t... fail\n"; exit 0; 
#	   	   }else{ print "  $mol frequencies (@nf{$mol})----@freq{$mol}\t... OK\n";
	};
      if (!@degeneration{$mol}) { print "  $mol degeneration (@deg_comment{$mol})\t... fail\n"; exit 0; 
#			   }else{ print "  $mol degeneration (@degeneration{$mol})\t... OK\n";
	};
      if (!@sitetype{$mol}) { print "  $mol adsorption sites\t... fail\n"; exit 0;
#		      }else{ print "  $mol adsorption sites (@sitetype{$mol})\t... OK\n";
	};
      if (!@Tmass{$mol}) { print"  $mol atomic mass\t... fail\n"; exit 0;
#		    }else{ print"  $mol atomic mass (@imass{$mol}; @Tmass{$mol})\t... OK\n";
	};
      if (!@inatoms{$mol}) { print"  $mol number of atoms\t... fail\n"; exit 0;
#		      }else{ print"  $mol number of atoms (@inatoms{$mol})\t... OK\n";
	};
      if (!@sym{$mol}) { print"  $mol symmetry factor\t... fail\n"; exit 0;
#		  }else{ print"  $mol symmetry factor (@sym{$mol})\t... OK\n";
	};
      if (!@inertia{$mol}) { print"  $mol inertia detection\t... fail\n"; exit 0;
#		      }else{ print"  $mol inertia detection\t... OK\n";
	};
   }; # foreach molecules
   foreach $sur (@surfaces) {
      if (!$e0{$sur}) { print "  $sur energy\t... fail\n"; exit 0; };
      if (!@freq{$sur}) { print "\n  $sur frequencies (@nf{$sur})\t... fail (continue)\n";};
      if (!@degeneration{$sur}) { print "  $sur degeneration (@deg_comment{$sur})\t... fail\n"; exit 0; };
      if (!@Acat{@sitetype{$sur}}) { print "  $sur surface area for @sitetype{$sur} site\t... fail\n"; exit 0; };
      if (!@sitetype{$sur}) { print "  $sur active sites selection\t... fail\n"; exit 0; };
   }; # foreach surfaces   
   foreach $cat (@catalysts) {
      if (!$e0{$cat}) { print "  $cat energy\t... fail\n"; exit 0; };
      if (!@freq{$cat}) { print "  $cat frequencies (@nf{$cat})\t... fail\n"; exit 0; };
      if (!@degeneration{$cat}) { print "  $cat degeneration (@deg_comment{$cat})\t... fail\n"; exit 0; };
      if (!@Acat{@sitetype{$cat}}) { print "  $cat active sites selection\t... fail\n"; exit 0; };
      if (!@sitetype{$cat}) { print "  $cat active sites selection\t... fail\n"; exit 0; };
   }; # foreach      
print "\t\t... done\n";   
#============================================================================================================================ 
#
#                                                   THEMODYNAMICS PROPERTIES
#
#============================================================================================================================
print "  Writing thermodynamics.m";
   foreach $mol (@molecules) { ()=&thermo_sub($mol,"mol"); };
   foreach $sur (@surfaces) { ()=&thermo_sub($sur,"sur"); };
   foreach $cat (@catalysts) { ()=&thermo_sub($cat,"cat"); };
#------------------------------------------------------------------------------------------------------------------------------ creating foldes
       system("mkdir -p $rootpath/THERMODYNAMICS/DATA $rootpath/KINETICS/DATA $rootpath/KINETICS/PROCESS $rootpath/KINETICS/PLOTS");
       mkdir("$rootpath/XYZ/");
print "\t... done\n";
#==================================================================================================================================
#
#
#
#
#                                                       READING KINETICS
#
#
#
#
#==================================================================================================================================
    $nprocess=1;
#==================================================================================================================================
print "Reading kinetics";
  foreach $line (@in) {
    if ( $line ) { @tag=split(/\s+/,$line);
       foreach $t (@tag) { if (($t) or (looks_like_number($t))) { push(@Nline,$t); }; }; @tag=@Nline; @Nline=();
       	if ((@tag[0] eq 'PROCESS') or (@tag[0] eq 'PROCESS=')) {
	   	@tmp=split(/=/,$line); @itmp=split(/\s+/,@tmp[1]); $t=0; $c=0; @proc=();
           foreach $t (@itmp) { if ($t) { push(@Nline,$t); }; }; @itmp=@Nline; @Nline=();
           while (($c ne '#') and (@itmp[$t])) { push(@proc,@itmp[$t]); $t++; @cc=split(//,@itmp[$t]); $c=@cc[0]; };
              @process[$nprocess]="@proc"; $nprocess++; };# if PROCESS
       	if ((@tag[0] eq 'STOICHIOMETRY') or (@tag[0] eq 'STOICHIOMETRY=')) { 
	       @tmp=split(/=/,$line); @itmp=split(/\s+/,@tmp[1]); $t=0;$c=0; @stoi=();
           foreach $t (@itmp) { if ($t) { push(@Nline,$t); }; }; @itmp=@Nline; @Nline=();
	   while (($c ne '#') and (@itmp[$t])) { push(@stoi,@itmp[$t]); $t++; @cc=split(//,@itmp[$t]); $c=@cc[0];};
	      @stoichiometry[$#process]="@stoi"; };# if STOICHIOMETRY
	if (@tag[0] eq 'EEXCH') { @ineexch[$#process]=@tag[2]; };
    }; #--> if line
  }; #--> line
    print "\t\t... done\n";
#---------------------------------------------------------------------------------------------------------------------------------coments
   foreach $pro (@process) { @tmp=split(/\s+/,@pro); $typeproc=$tmp[0];
      if (($typeproc != 'A') or ($typeproc != 'a') or ($typeproc != 'D') or ($typeproc != 'd') or ($typeproc != 'R') or ($typeproc != 'r')) {
         print "\n -  Please, provide the process type for $pro - \n\n" ; exit 0 ;}; }; #--> foreach process
#=================================================================================================================== STARTING MICROKINETIC ANALYSIS
#                LOOKING FOR CONSTANT RATES DEPENDING ON:
#                                                       TEMPERATURE (T)
#                                                       EXTERNA POTENTIAL (V)
#                                                       pH BATH CONDITIONS (pH)
#
#===================================================================================================================
#-------------------------------------------------------------------------------------------------------------------------------- react_species 
          @react_species=(); @transitions=();
   for ($pr=1; $pr<=$#process; $pr++ ) { 
       ($PR,$PTS,$PP,$type)=&process_sub(@process[$pr],@stoichiometry[$pr]);
          @ProcessReactants[$pr]="@$PR"; @ProcessTS[$pr]="@$PTS"; @ProcessProducts[$pr]="@$PP"; @typeproc[$pr]=$type;
          @PR=split(/\s+/,@ProcessReactants[$pr]); @PTS=split(/\s+/,@ProcessTS[$pr]); @PP=split(/\s+/,@ProcessProducts[$pr]);
      foreach $p (@PR) { $new="yes"; foreach $rs (@react_species) { if ($rs eq $p) { $new="no"; };};
         if ($new eq "yes") { push(@react_species,$p); };};
      foreach $p (@PTS) { $new="yes"; foreach $rs (@react_species) { if ($rs eq $p) { $new="no"; };};
         if ($new eq "yes") { push(@react_species,$p); push(@transitions,$p); };};
      foreach $p (@PP) { $new="yes"; foreach $rs (@react_species) { if ($rs eq $p) { $new="no"; };};
         if ($new eq "yes") { push(@react_species,$p); };};
   }; # for process

#-------------------------------------------------------------------------------------------------------------------------------- Import into processes
    ()=&Introduction_sub("processes");
         @done=();
      foreach $mol (@molecules) { $tmp=1; foreach $did (@done) { if ($did eq $mol) { $tmp=0; };};
         if ($tmp eq 1) { ()=&Import_sub($mol,"processes"); push(@done,$mol); };};
      foreach $sur (@surfaces) { $tmp=1; foreach $did (@done) { if ($did eq $sur) { $tmp=0; };};
         if ($tmp eq 1) { ()=&Import_sub($sur,"processes"); push(@done,$sur); };};
      foreach $cat (@catalysts) { $tmp=1; foreach $did (@done) { if ($did eq $cat) { $tmp=0; };};
         if ($tmp eq 1) { ()=&Import_sub($cat,"processes"); push(@done,$cat); };};
#------------------------------------------------------------------------------------------------------------------------------- Interpolation
            @interpolated=();
      foreach $species (@react_species) { $tmp=0; 
          foreach $mol (@molecules) { if ($mol eq $species) { $tmp=1; };};
          foreach $sur (@surfaces) { if ($sur eq $species) { $tmp=1; };};
	  foreach $cat (@catalysts) { if ($cat eq $species) { $tmp=1; };};
	  foreach $in (@interpolated) { if ($in eq $species) { $tmp=1; };};
          if ($tmp eq 0) { ($E,$Q3D)=&Interpolate_sub($species); 
		  push(@interpolated,$species); @en{$species}=$E; @q{$species}=$Q3D; };};
#-------------------------------------------------------------------------------------------------------------------------------- processes
print "  Writing processes.m";
           $row=1; %y=();
         foreach $mol (@molecules) { $go="no"; 	foreach $rs (@react_species) { if ($rs eq $mol) { $go="yes"; };};
		 				foreach $t (@transitions) { if ($mol eq $t) { $go="no"; };};
		 if ($go eq "yes") { if (!@y{$mol}) { @y{$mol}=$row; $row++; };};}; 
         foreach $rs (@react_species) { $go="yes"; $tr="no"; foreach $mol (@molecules) { if ($rs eq $mol) { $go="no"; };};
					                     foreach $t (@transitions) { if ($rs eq $t) { $go="no"; };};
							     foreach $sur (@surfaces) { if ($rs eq $sur) { $go="no"; };};
                if ($go eq "yes") { if (!@y{$rs}) { @y{$rs}=$row; $row++; };};}; 
         foreach $sur (@surfaces) { if (!@y{$sur}) { @y{$sur}=$row; $row++; };};
     for ($pr=1; $pr<=$#process; $pr++ ) {
           @PR=split(/\s+/,@ProcessReactants[$pr]); @PTS=split(/\s+/,@ProcessTS[$pr]); @PP=split(/\s+/,@ProcessProducts[$pr]);
                 	         ()=&stoichiometries_sub("processes");
                                 ()=&ProcessE_sub($typeproc[$pr],$pr);
	                         ()=&ProcessQ_sub($typeproc[$pr],$pr);
                                 ()=&ProcessK_sub($typeproc[$pr],$pr);
			         ()=&ProcessParameters_sub($typeproc[$pr],$pr); }; 
print "\t\t... done\n";			 
#=================================================================================================================== DIFFERENTIAL EQUATION SYSTEM
##
##
##
##
##===================================================================================================================			 
       @experiments=("const_TEMP", "variable_TEMP", "TPR", "RateControl");
#-------------------------------------------------------------------------------------------------------------------------------- Import into ODE
     foreach $exp (@experiments) {
print "  Writing $exp.m";	   
	if (($exp eq "const_TEMP") or ($exp eq "variable_TEMP")) { 
          	()=&Introduction_sub($exp);
	        ()=&ODE_import_sub($exp);
         	($IC)=&variables_sub($exp); $IC="@$IC"; @IniCon=split(/\s+/,$IC);
       		foreach $i (@IniCon) { if ($i) { push(@Nline,$i); }; }; @IniCon=@Nline; @Nline=();
          	foreach $v (@variables) { 
             		foreach $rs (@react_species) { $go="no";
                		if ($v eq $rs) {foreach $mol (@molecules) { if ($rs eq $mol) { $go="no"; };};
						foreach $t (@transitions) { if ($rs eq $t) { $go="no"; };};
                                 		foreach $sur (@surfaces) { if ($rs eq $sur) { $go="no"; };};
                                 	if ($go eq "yes") { $new="yes"; foreach $ic (@IniCon) { if ($ic eq $v) {$new="no"; };};
                                    	if ($new eq "yes") { push(@IniCon,$v); };};};};};
        	()=&ODE_call_sub($exp);
         	()=&ODE_solution_sub($exp);
         	()=&ODE_closingloops_sub($exp);
	 	()=&ODE_printing_sub($exp);
		if ($exp eq "const_TEMP") {
print "\t\t... done\n";	
		}else{
print "\t... done\n";
                };			
	}elsif ($exp eq "TPR") {
                ()=&Introduction_sub($exp);
                ()=&ODE_import_sub($exp);
		@TPRmol=();
                for ($pr=1; $pr<=$#process; $pr++ ) {
                        if ((@typeproc[$pr] eq "A") or (@typeproc[$pr] eq "a")) {
                                @PR=split(/\s+/,@ProcessReactants[$pr]);
                                foreach $R (@PR) { foreach $mol (@molecules) { if ($R eq $mol) { push(@TPRmol,$mol); };};};};};
		foreach $TPRmolecule (@TPRmol) {
                	($IC)=&variables_sub($exp); $IC="@$IC"; @IniCon=split(/\s+/,$IC);
                	foreach $i (@IniCon) { if ($i) { push(@Nline,$i); }; }; @IniCon=@Nline; @Nline=();
        	        foreach $v (@variables) {
	                        foreach $rs (@react_species) { $go="no";
                        	        if ($v eq $rs) {foreach $mol (@molecules) { if ($rs eq $mol) { $go="no"; };};
                      	                  		foreach $t (@transitions) { if ($rs eq $t) { $go="no"; };};
                                        		foreach $sur (@surfaces) { if ($rs eq $sur) { $go="no"; };};
                                	        if ($go eq "yes") { $new="yes"; foreach $ic (@IniCon) { if ($ic eq $v) {$new="no"; };};
                        	                if ($new eq "yes") { push(@IniCon,$v); };};};};};
                	()=&ODE_call_sub($exp);
                	()=&ODE_solution_sub($exp);
                	()=&ODE_closingloops_sub($exp);
        	        ()=&ODE_printing_sub($exp);		
		}; #foreach TPR
print "\t\t\t... done\n";

# Alberto 02/2019
        }elsif ($exp eq "RateControl") {
                ()=&Introduction_sub($exp);
		()=&ODE_import_sub($exp);
		for ($pr=1; $pr<=$#process; $pr++ ) {
                        @PR=split(/\s+/,@ProcessReactants[$pr]); @PTS=split(/\s+/,@ProcessTS[$pr]); @PP=split(/\s+/,@ProcessProducts[$pr]);
                        ()=&stoichiometries_sub($exp);
                };
                ()=&DRCrates_sub($exp,$pr);
                ($IC)=&variables_sub($exp); $IC="@$IC"; @IniCon=split(/\s+/,$IC);
                foreach $i (@IniCon) { if ($i) { push(@Nline,$i); }; }; @IniCon=@Nline; @Nline=();
                foreach $v (@variables) {
                        foreach $rs (@react_species) { $go="no";
                                if ($v eq $rs) {foreach $mol (@molecules) { if ($rs eq $mol) { $go="no"; };};
                                        	foreach $t (@transitions) { if ($rs eq $t) { $go="no"; };};
                                        	foreach $sur (@surfaces) { if ($rs eq $sur) { $go="no"; };};
                                        if ($go eq "yes") { $new="yes"; foreach $ic (@IniCon) { if ($ic eq $v) {$new="no"; };};
                                        if ($new eq "yes") { push(@IniCon,$v); };};};};};


                ()=&ODE_call_sub($exp);
                ()=&ODE_solution_sub($exp);
		($EQ)=&equations_sub($exp); $EQ="@$EQ"; @Equa=split(/\s+/,$EQ);
		()=&DRC_sub($exp,@Equa); 

print "\t\t... done\n";
	}; # if exp 
    }; # foreach exp
       ()=&ODE_sub("noTPR");
       ()=&ODE_sub("TPR");  



#=================================================================================================================================





#==============================================================================================================================
#
#                               EEEEEEE	N     N DDD
#				E	NN    N D  D
#				E	N N   N D   D
#				EEE	N  N  N D    D
#				E	N   N N D   D
#				E	N    NN D  D
#				EEEEEEE	N     N DDD
#
#==============================================================================================================================

	     
	
	
	
#=================================================================================================================================	
#==============================================================================================================================
# This sub looks for the symmetry group assigned in the input (e.g. OUTCAR)

   sub symmetry_sub {
      ($ipath)=@_;
     system ("grep 'The static configuration has the point symmetry' $ipath |tail -1 >>temp.prop"); 
     if (-z temp.prop) { system ("grep ' symmetry ' $ipath |tail -1 >>temp.prop"); };
        @input=(); $prop='';
     open IN, "temp.prop"; while (<IN>) { $prop.= $_ ;} close IN;
         @stmp=split(/\s+/,$prop);
       foreach $l (@stmp) { if ($l) { push(@Nline,$l); }; }; @stmp=@Nline; @Nline=();
       if (@stmp[0] eq 'The') { $Gsym=@stmp[7]; }elsif (@stmp[0] eq 'NOSYMM:'){ $Gsym=@stmp[9]; };
       system ("rm temp.prop");
     return ($GlobalSym);
   }; #--> sub symmetry
#==============================================================================================================================
# This sub looks for the symmetry factor according to symmetry_sub and the linearity of the structure.

   sub symfactor_sub {
       ($ipath,$inatoms,$linear)=@_;
         @natoms=split(/\s+/,$inatoms);
         ($Globalsym)=&symmetry_sub($ipath);
         $tatoms=$#natoms+1;
         if ($linear eq "yes") { if ($tatoms == 1) { $symfac=2; };
                                 if ($tatoms >= 2) { $symfac=1; };
         }else{ @tmp=split(//,$symme); if (@tmp[2] == 1) { 
			         if ((@tmp[3] eq "h") or (@tmp[3] eq "H")) { $symfac=@tmp[2]+1; 
		                 }else{ $symfac=@tmp[2]; };};};#--> if
	  return ($symfac);
   }; #--> sub symfactor
#==============================================================================================================================
#==============================================================================================================================
# This sub prints the properties and look for the themodynamic elements of the system.

   sub thermo_sub {
        ($sys,$imol)=@_; @ifreq=split(/\s+/,@freq{$sys}); @ifreq2D=split(/\s+/,@freq2D{$sys}); $fsum=0;
	      foreach $f (@ifreq) { if (($f) and (looks_like_number($f)) and ($f != 0)) { push(@Nline,$f); }; }; @ifreq=@Nline; @Nline=();
              foreach $f (@ifreq2D) { if (($f) or (looks_like_number($f))) { push(@Nline,$f); }; }; @ifreq2D=@Nline; @Nline=();
      open OUT, ">>thermodynamics.m";
          print OUT "%  ------------------- $sys ----------------\n";
          print OUT "   mkdir THERMODYNAMICS/DATA $sys; mkdir THERMODYNAMICS/PLOTS $sys;\n";
	  print OUT "clearvars;\n\t\t\t%_________________General Constants__________________\n";
	  print OUT "h=6.62607588705515e-34; kb=1.380658045430573e-23; c=29979299840.00; R=8.31451034545898;\n";
          print OUT "Av=6.022139922973909e23; Fa=96485.296875000; toeV=1.60217648740e-19; electron=-1.60217648740e-19; JtoeV=6.24150974e18;\n\n";
	  print OUT "%________________________________________________T is considered as variable (symbolic)\n";
	  print OUT "syms T;\n"; @variable=();
	  print OUT "     fileID=fopen(\'./THERMODYNAMICS/DATA/$sys/summary.dat\','a+');\n";
	  print OUT "     fprintf(fileID, \'SYSTEM = $sys\\n\');\n";
        if ($imol eq "mol") {
	     print OUT "\t\t\t%___________________System properties extracted from the $sys 's input file___________________\n";	
      	     print OUT "pressure$sys=101325;\t\t\tfprintf(fileID, \' pressure = %d Pa\\n\',pressure$sys);\n";
             print OUT "Mass$sys=@Tmass{$sys};\tfprintf(fileID, \' molecular mass = %.15E Kg\\n\',Mass$sys);\n";
             print OUT "symmetry$sys=@sym{$sys};\t\t\tfprintf(fileID, \' molecular symmetry = %d\\n\',symmetry$sys);\n"; 
             print OUT "Vol$sys=kb*T/pressure$sys;\t\t%fprintf(fileID, \' molecular volume = %.15E dm^3 ~eV/Pa\\n\',Vol$sys);\n";
	         push(@variable,"Vol$sys"); # eV/Pa ~ dm^3
	     print OUT "\t\t\t\t\tfprintf(fileID, \' linear molecule = @linear{$sys}\\n\');\n";
	     print OUT "inertia$sys=@inertia{$sys};\tfprintf(fileID, \' molecular moment of inertia = %.15E kg m^2\\n\',inertia$sys);\n";
	     print OUT "  area@sitetype{$sys}=@Acat{@sitetype{$sys}};\t\tfprintf(fileID, \' area on the surface = %.15E m^2\\n\',area@sitetype{$sys});\n"; 
             print OUT "  Nsites$sys=@nsitetype{$sys};\t\t\tfprintf(fileID, \' sites on the surface = %f\\n\',Nsites$sys);\n";
	     print OUT "\t\t\t%_________________Translational Partition Function (qtrans) as defined in DOI:10.1002/3527602658_________________\n";
             print OUT "    qtrans3D$sys=Vol$sys*((2*pi*Mass$sys*kb*T)^(3/2))/(h^3);	% free molecule moving in 3 dimensions\n";
             print OUT "    qtrans2D$sys=((area@sitetype{$sys}*Nsites$sys)*2*pi*Mass$sys*kb*T)/(h^2);	% free molecule moving in 2 dimensions (hovering on surface)\n";
                 push(@variable,"qtrans3D$sys");  push(@variable,"qtrans2D$sys");
	     print OUT "\t\t\t%_____________Rotational Partition Function (qrot) as defined in DOI:10.1002/3527602658_____________\n";
            if (@linear{$sys} eq "yes"){
             print OUT "    qrot3D$sys=(8*pi^2*inertia$sys*kb*T)/(symmetry$sys*h^2);\n";
#              print OUT "    qrot3D$sys=(sqrt(pi)/symmetry$sys)*(8*pi^2*kb*T*inertia$sys/(h^2))^(3/2);\n";
           }elsif (@linear{$sys} eq "no"){
             print OUT "    qrot3D$sys=(1/symmetry$sys)*(8*pi^2*kb*T/(h^2))^(3/2)*(sqrt(pi*inertia$sys));\n"; };
                 push(@variable,"qrot3D$sys");
        }else{
           foreach $sur (@surfaces) { if (@sitetype{$sys} eq @sitetype{$sur}) {
             print OUT "area@sitetype{$sys}=@Acat{@sitetype{$sur}};\tfprintf(fileID, \'surface area = %.15E\\n\',area@sitetype{$sys});\n";
	     $done=1; };
             print OUT "\t\t\t%_____________qtrans and qrots for a solid that does not move/rotate is one_____________\n";
             print OUT "    qtrans3D$sys=1;\n"; push(@variable,"qtrans3D$sys");  
     	     print OUT "    qrot3D$sys=1;\n";   push(@variable,"qrot3D$sys"); };
        }; # if molecule
	     print OUT "\t\t\t%_____________Electronic Partition Function (qelec) is the spin multiplicity\n";
             print OUT "\t\t\t%_____________qelec consideres that the working temperatures are small enough\n";
             print OUT "\t\t\t%_____________to neglect excited states. Qelec is the same in the gas and in the TS.\n";
             print OUT "    qelec$sys=@degeneration{$sys};\t\t\tfprintf(fileID, \'electronic ground state multiplicity = %f\\n\',@degeneration{$sys});\n";
		     push(@variable,"qelec$sys");
        if ($#ifreq >= 50) { $nqvib=int($#ifreq/50); }else{ $nqvib=0; };
	if ($imol eq "mol") {
	     print OUT "\t\t\t%_____________Vibrational Partition Function (qvib) as defined in DOI:10.1002/3527602658_____________\n";
	     print OUT "\t\t\t%_____________qvib considers the temperature independent contribution as ZPE,\n";
             print OUT "\t\t\t%_____________which is directly added to the Energy and removed from qvib\n";
		if ($nqvib > 0) { $tmp=0;
			while ( $tmp <= $nqvib ) {
				print OUT "     qvib3D$sys$tmp="; $nf=0;
				for ($k=50*$tmp; $k<=49+50*$tmp; $k++) {
				       if (($k < 49+50*$tmp) and ( $nf< $#ifreq)) {
                                       		print OUT "1/(1-exp(-h*c*$ifreq[$k]/(kb*T)))*"; $nf++; }else{ if (@ifreq[$k]) {
                                                print OUT "1/(1-exp(-h*c*$ifreq[$k]/(kb*T)));\n"; $nf++; };};};
				$tmp++; };
			print OUT "  qvib3D$sys="; for ($k=0; $k<=$nqvib-1; $k++) {
						print OUT "qvib3D$sys$k*"; };
						print OUT "qvib3D$sys$nqvib;\n";
		}else{  print OUT "  qvib3D$sys="; $nf=0; foreach $f (@ifreq) { if ( $nf < $#ifreq) {
                                                     print OUT "1/(1-exp(-h*c*$f/(kb*T)))*"; $nf++; }else{
                                                     print OUT "1/(1-exp(-h*c*$f/(kb*T)));\n"; $nf++;};};};
			                push(@variable,"qvib3D$sys");
	        if ($#ifreq2D >= 50) { $nqvib2D=int($#ifreq2D/50); }else{ $nqvib2D=0; };
                if ($nqvib2D > 0) { $tmp=0;
	    print OUT "\t\t\t%_____________2D Rotational Partition Function consider only frequencies that does not\n";
            print OUT "\t\t\t%_____________shift the centre of mass along the reaction coordinate (Z)\n";
                        while ( $tmp <= $nqvib2D ) {
                                print OUT "     qvib2D$sys$tmp="; $nf=0;
                                for ($k=50*$tmp; $k<=49+50*$tmp; $k++) {
                                       if (($k < 49+50*$tmp) and ( $nf < $#ifreq2D)) {
                                                print OUT "1/(1-exp(-h*c*$ifreq2D[$k]/(kb*T)))*"; $nf++; }else{ if (@ifreq2D[$k]) {
                                                print OUT "1/(1-exp(-h*c*$ifreq2D[$k]/(kb*T)));\n"; $nf++; };};};
                                $tmp++; };
                        print OUT "  qvib2D$sys="; for ($k=0; $k<=$nqvib2D-1; $k++) {
                                                print OUT "qvib2D$sys$k*"; };
                                                print OUT "qvib2D$sys$nqvib;\n";
		}else{  print OUT "  qvib2D$sys="; $nf=0; if (!$ifreq2D[0]) { print OUT "1;\n"; }else{
	                                foreach $f (@ifreq2D) { if ($nf < $#ifreq2D) {
	                                             print OUT "1/(1-exp(-h*c*$f/(kb*T)))*"; $nf++; }else{
	                                             print OUT "1/(1-exp(-h*c*$f/(kb*T)));\n"; $nf++; };};};};
					                push(@variable,"qvib2D$sys");  
	     print OUT "\t\t\t%_____________Q3Dnotrans is for calculating molecules transition (TS) towards an\n";
             print OUT "\t\t\t%_____________anchored site, i.e. a trasition state in indirect adsorption mode.\n";
             print OUT "\t\t\t%_____________In case of an existent TS defined, it will use it to model direct adsorptions.\n";
             print OUT "  Q3Dnotrans$sys=qvib3D$sys*qrot3D$sys;\n";   #------------ no *qelec@nsys[$i] because qelec(#)=qelec(gas)
	         push(@variable,"Q3Dnotrans$sys");
         }else{
                if ($nqvib > 0) { $tmp=0;
                        while ( $tmp <= $nqvib ) {
                                print OUT "     qvib3D$sys$tmp="; $nf=0;
                                for ($k=50*$tmp; $k<=49+50*$tmp; $k++) { if (@ifreq[$k]) {
                                       if (($k < 49+50*$tmp) and ( $nf < $#ifreq)) {
                                                print OUT "1/(1-exp(-h*c*$ifreq[$k]/(kb*T)))*"; $nf++; }else{
                                                print OUT "1/(1-exp(-h*c*$ifreq[$k]/(kb*T)));\n"; $nf++; };};};
                                $tmp++; };
                        print OUT "  qvib3D$sys="; for ($k=0; $k<=$nqvib-1; $k++) {
                                                print OUT "qvib3D$sys$k*"; };
                                                print OUT "qvib3D$sys$nqvib;\n";
                }else{  print OUT "  qvib3D$sys="; if ((!$ifreq[0]) or ($ifreq[0] == 0)) { print OUT "1;\n";
			}else{ $nf=0; foreach $f (@ifreq) { if ($nf < $#ifreq) {
                                                print OUT "1/(1-exp(-h*c*$f/(kb*T)))*"; $nf++; }else{
                                                print OUT "1/(1-exp(-h*c*$f/(kb*T)));\n"; $nf++; };};};};
                                        push(@variable,"qvib3D$sys");
       	};
             print OUT "\t\t\t%_____________The total partition function in the three dimensions (Q3D)\n";
             print OUT "\t\t\t%_____________is the productory of previous q\n";
             print OUT "  Q3D$sys=qvib3D$sys*qrot3D$sys*qtrans3D$sys*qelec$sys;\n";
	        push(@variable,"Q3D$sys");
#-------------------------------------------------------------------------------------------------------------------------------- ZPE          
	if ($imol eq "mol") {
             $fsum2=0;
                 if ($nqvib2D > 0) { $tmp=0;
                        while ( $tmp <= $nqvib2D ) {
             print OUT "\t\t\t%_____________The ZPE is the vibrational contribution to the energy at T=0K_____________\n";
             print OUT "\t\t\t%_____________ZPE2D derives from these vibrations parallel to the adsorbent (X and Y)\n";
             print OUT "\t\t\t%_____________as the third dimension (Z) correspond to the reaction coordinate of an\n";
             print OUT "\t\t\t%_____________indirect adsorption mode.\n";
                                print OUT "     ZPE2D$sys$tmp="; $nf=0;
                                for ($k=50*$tmp; $k<=49+50*$tmp; $k++) {
                                       if (($k < 49+50*$tmp) and ( $nf< $#ifreq2D)) { 
 					        print OUT "(sinh($ifreq2D[$k]*(h*c)/(2*kb*T))/($ifreq2D[$k]*(h*c)/(2*kb*T)))*";
					       	$fsum2=$fsum2+$ifreq2D[$k]; $nf++; }else{ if (@ifreq2D[$k]) {
					        print OUT "(sinh($ifreq2D[$k]*(h*c)/(2*kb*T))/($ifreq2D[$k]*(h*c)/(2*kb*T)));\n";
					       	$fsum2=$fsum2+$ifreq2D[$k]; $nf++; };};};
				$tmp++; };
			print OUT "  ZPE2D$sys=(kb*T/toeV)*log("; for ($k=0; $k<=$nqvib2D-1; $k++) {
                                                print OUT "ZPE2D$sys$k*"; };
                                                print OUT "ZPE2D$sys$nqvib2D);\n";
		}else{  print OUT "  ZPE2D$sys=(kb*T/toeV)*log(";
 	     if (!$ifreq2D[0]) { print OUT "1);\n"; }else{ $nf=0; foreach $f (@ifreq2D) { if ($nf < $#ifreq2D) {
                                                print OUT "(sinh($f*(h*c)/(2*kb*T))/($f*(h*c)/(2*kb*T)))*"; $fsum2=$fsum2+$f; $nf++; }else{
                                                print OUT "(sinh($f*(h*c)/(2*kb*T))/($f*(h*c)/(2*kb*T))));\n"; $fsum2=$fsum2+$f; $nf++; };};};};
             print OUT "\t\t\t%_____________ZPEClassic does not considers the Wigner Zero Point Correction as in DOI: 10.1063/1.2161193\n";
   	     print OUT "  ZPE2Dclassic$sys=(1/2)*(h*c/toeV)*$fsum2;\n\t\t\t\t\tfprintf(fileID, \'  classic ZPE 2D = %.15E\\n\',ZPE2Dclassic$sys);\n"; 
	     push(@variable,"ZPE2D$sys"); 
     	};
         $fsum=0;
         if ($nqvib > 0) { $tmp=0;
        	while ( $tmp <= $nqvib ) {
                	print OUT "     ZPE$sys$tmp="; $nf=0;
                                for ($k=50*$tmp; $k<=49+50*$tmp; $k++) {
                                       if (($k < 49+50*$tmp) and ( $nf < $#ifreq)) { 
                                                print OUT "(sinh($ifreq[$k]*(h*c)/(2*kb*T))/($ifreq[$k]*(h*c)/(2*kb*T)))*";
                                                $fsum=$fsum+$ifreq[$k]; $nf++; }else{ if (@ifreq[$k]) {
                                                print OUT "(sinh($ifreq[$k]*(h*c)/(2*kb*T))/($ifreq[$k]*(h*c)/(2*kb*T)));\n";
                                                $fsum=$fsum+$ifreq[$k]; $nf++; };};};
                                $tmp++; };
                        print OUT "  ZPE$sys=(kb*T/toeV)*log(";for ($k=0; $k<=$nqvib-1; $k++) {
                                                print OUT "ZPE$sys$k*"; };
                                                print OUT "ZPE$sys$nqvib);\n";
                }else{  print OUT "  ZPE$sys=(kb*T/toeV)*log("; $fsum=0;
		     if (!$ifreq[0]) { print OUT "1);\n"; }else{ $nf=0; foreach $f (@ifreq) { if ($nf < $#ifreq) {
					        print OUT "(sinh($f*(h*c)/(2*kb*T))/($f*(h*c)/(2*kb*T)))*";
					       	$fsum=$fsum+$f; $nf++; }else{
           					print OUT "(sinh($f*(h*c)/(2*kb*T))/($f*(h*c)/(2*kb*T))));\n";
					       	$fsum=$fsum+$f; $nf++; };};};};
           push(@variable,"ZPE$sys");
	   print OUT "  ZPEclassic$sys=(1/2)*(h*c/toeV)*$fsum;\n\t\t\t\t\tfprintf(fileID, \'  classic ZPE 3D = %.15E\\n\',ZPEclassic$sys);\n";   

#-------------------------------------------------------------------------------------------------------------------------------- @S Entropy
# pag 88; Chorkendorff, I. a. N., J. W. Concept... f Modern Catalysis and Kinetics; Wiley-VCH: Weinheim, 2005. ---> S={Kb*ln(Q)}
           print OUT "\t\t\t%_____________Entroy (S) is S={Kb*ln(Q)} as as defined in DOI:10.1002/3527602658 (page 88)_____________\n";
#           print OUT "S$sys=(kb*log(Q3D$sys)+(kb*T*diff(log(Q3D$sys),T)))/toeV;\n";  # eV/K
# re-Formulation from explicit derivatives:: https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html
    if ($imol eq "mol") {
           print OUT "   Strans$sys=kb*(log((2*pi*Mass$sys*kb*T)^(3/2)/(h^3)*Vol$sys)+5/2)*JtoeV;\n";
           if (@linear{$sys} eq "yes") {
               print OUT "   Srot$sys=kb*(log((1/symmetry$sys)*(8*pi^2*kb*T*inertia$sys/(h^2)))+1)*JtoeV;\n";
           }
           elsif (@linear{$sys} eq "no") {
               print OUT "   Srot$sys=kb*(log((sqrt(pi*inertia$sys)/symmetry$sys)*(8*pi^2*kb*T/(h^2))^(3/2))+3/2)*JtoeV;\n";};
    }else{ print OUT "   Strans$sys=0;\n   Srot$sys=0;\n"; };
 		   print OUT "   Svib$sys=";
                if ((! @ifreq) or ($#ifreq == 0)) { print OUT "0;\n"; }else{ $nf=0;
                             print OUT "kb*("; foreach $f (@ifreq) { if ($nf < $#ifreq) {
                             print OUT "(h*c*$f)/(kb*T*(exp(h*c*$f/(kb*T)))) - log(1-exp(-h*c*$f/(kb*T)))+"; $nf++;
                                   }else{
                             print OUT "(h*c*$f)/(kb*T*(exp(h*c*$f/(kb*T)))) - log(1-exp(-h*c*$f/(kb*T))) )*JtoeV;\n"; $nf++; };
                             };};
           print OUT "   Selec$sys=kb*(log(2*@degeneration{$sys}))*JtoeV;\n";

           print OUT "S$sys=Strans$sys+Srot$sys+Svib$sys+Selec$sys;     % -kb*log(P/P0)\n";
	      push(@variable,"S$sys");
#-------------------------------------------------------------------------------------------------------------------------------- @Cp
# pag 551; Principles of Physical Chemistry, By Hans Kuhn, Horst-Dieter Försterling, David Hennessey Waldeck   ---> Cp=T*[dS/dT](N,P)
           print OUT "\t\t\t%_____________Entroy (S) is Cp=T*[dS/dT](N,P) as as defined in ISBN: 9780470089644 (page 551)_____________\n";
           print OUT "Cp$sys=T*diff(S$sys,T); funcH$sys=matlabFunction(Cp$sys);\n";
	      push(@variable,"Cp$sys"); 
#------------------------------------------------------------------------------------------------------------------------------- @H Enthalpy
# pag 91; PHYSICAL CHEMISTRY for Chemical and Biological Science, By Raymand Chang                             ---> H=[dCp/dT](T1,T2)
#           print OUT "    H0$sys=@e0{$sys}+ZPE$sys+int(Cp$sys,T);\n";  push(@variable,"H0$sys");
           print OUT "\t\t\t%_____________Enthalpy (H) is H=[dCp/dT](T1,T2) as as defined in ISBN: 1-891389-06-8 (page 91)_____________\n";
           print OUT "\t\t\t%_____________The integral is solved with a relative and absolute tolerance of 1E-3 and 1E-4 respectively\n";
           print OUT "\t\t\t%_____________In the case of treating with molecular O2(g), the DFT bond\n";
           print OUT "\t\t\t%_____________overbinding correction is applied according to DOI: 10.1039/c4cp00529e\n";
          if ($sys eq "O2") { print OUT "% H$sys=$e0{$sys}+0.9/2+MZPE$sys(i,2)+integral(funcH$sys,0,T,'RelTol', 1e-3, 'AbsTol', 1e-4); % O2 overbinding correction (O2=-5.17eV)\n";
  	    push(@variable,"H$sys"); }else{ print OUT "% H$sys=$e0{$sys}+MZPE$sys(i,2)+integral(funcH$sys,0,T,'RelTol', 1e-3, 'AbsTol', 1e-4);\n";  push(@variable,"H$sys");};
           print OUT "% G$sys=H$sys-T*S$sys;\n\n";     push(@variable,"G$sys");
#-------------------------------------------------------------------------------------------------------------------------------       
        if ($ttemp) { $size=($ftemp-$itemp)/$ttemp; }else{ $size=1; };
        for ($j=0; $j<=$#variable; $j++) { print OUT "  M$variable[$j]=zeros($size,2); "; }; print OUT "\n";
		if ($ttemp) { 	print OUT "i=1;\nfor T = $itemp:$ttemp:$ftemp\n"; 
		}else{ 		print OUT "i=1;\n T=$itemp;\n"; };
  	for ($j=0; $j<=$#variable-2; $j++) { print OUT "      M$variable[$j](i,:)=[T, subs($variable[$j])];\n"; };
           print OUT "   if T == 300\n\tfprintf(fileID, \' ZPE 3D (T=300K) = %.15E\\n\',MZPE$sys(i,2));\n";
                           print OUT "\tfprintf(fileID, \'  Qrotational 3D (T=300K) = %.15E\\n\',Mqrot3D$sys(i,2));\n";
                           print OUT "\tfprintf(fileID, \'  Qtranslational 3D (T=300K) = %.15E\\n\',Mqtrans3D$sys(i,2));\n";
                           print OUT "\tfprintf(fileID, \'  Qvibrational 3D (T=300K) = %.15E\\n\',Mqvib3D$sys(i,2));\n";
                           print OUT "\tfprintf(fileID, \'  Qelectronical (T=300K) = %.15E\\n\',Mqelec$sys(i,2));\n";
                           print OUT "\tfprintf(fileID, \' Q 3D (T=300K) = %.15E\\n\',MQ3D$sys(i,2));\n";
                           print OUT "\tfprintf(fileID, \' Specific hear (Cp) (T=300K) = %.15E eV/K\\n\',MCp$sys(i,2));\n";
                           print OUT "\tfprintf(fileID, \' Entropy (S) (T=300K) = %.15E eV/K\\n\',MS$sys(i,2));\n   end\n";

       	if ($sys eq "O2") { print OUT "  MH$sys(i,:)=[T, subs($e0{$sys}+0.9/2+MZPE$sys(i,2)+integral(funcH$sys,0,T,'RelTol', 1e-3, 'AbsTol', 1e-4))];\n";
        	}else{ 	if (!$ifreq[0]) {
		     		print OUT "  MH$sys(i,:)=[T, subs($e0{$sys}+MZPE$sys(i,2))]; % Cp = 0 -> no vibrations\n"; 
			}else{  print OUT "  MH$sys(i,:)=[T, subs($e0{$sys}+MZPE$sys(i,2)+integral(funcH$sys,0,T,'RelTol', 1e-3, 'AbsTol', 1e-4))];\n"; };};
           print OUT "  MG$sys(i,:)=[T, subs(MH$sys(i,2)-T*MS$sys(i,2))];\n";

           print OUT "   if T == 300\n\tfprintf(fileID, \'Enthalpy (H) (T=300K) = %.15E eV\\n\',MH$sys(i,2));\n";
                           print OUT "\tfprintf(fileID, \'Free energy (G) (T=300K) = %.15E eV\\n\',MG$sys(i,2));\n   end\n";
	if ($ttemp) {  print OUT "i=i+1;\nend\nfclose(fileID);\n";
		}else{ print OUT "fclose(fileID);\n"; };
        for ($j=0; $j<=$#variable; $j++) {
           print OUT "dlmwrite(\'./THERMODYNAMICS/DATA/$sys/$variable[$j].dat\',M$variable[$j],\'delimiter\',\'\\t\',\'precision\',\'%1.15E\');\n";
        }; # for variable 
           print OUT "\n";
      close OUT;
    return ();
   }; #--> sub thermo    
#==============================================================================================================================
   sub Introduction_sub {
	   ($exp)=@_;
      open OUT, ">>$exp.m";
         if ($exp eq "const_TEMP") {
            print OUT "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n%\n%\n";
	    print OUT "%      This script runs the kinetics at CONSTANT TEMPERATURE    \n%\n%\n";
	    print OUT "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
	    print OUT "mkdir KINETICS/DATA/$exp; mkdir KINETICS/PLOTS/$exp;\n";
	    print OUT "\nclearvars;\n    sol=[];\n";
         }elsif ($exp eq "variable_TEMP" ) {
            print OUT "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n%\n%\n";
            print OUT "%      This script runs the kinetics at VARIABLE TEMPERATURE    \n%\n%\n";
            print OUT "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
            print OUT "mkdir KINETICS/DATA/$exp; mkdir KINETICS/PLOTS/$exp;\n";
            print OUT "\nclearvars;\n    sol=[];\n";
         }elsif ($exp eq "TPR") {
            print OUT "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n%\n%\n";
            print OUT "%                       This script runs the                    \n";
            print OUT "%                TEMPERATURE PROGRAMMED DESORPTION               \n";
            print OUT "%                      at variable temperature                  \n%\n%\n";
            print OUT "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
            print OUT "mkdir KINETICS/DATA/$exp; mkdir KINETICS/PLOTS/$exp;\n";
            print OUT "\nclearvars;\n    sol=[];\n";
         }elsif ($exp eq "RateControl") {
            print OUT "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n%\n%\n";
            print OUT "%                       This script runs the                    \n";
            print OUT "%                DEGREE OF RATE and SELECTIVITY CONTROL         \n";
            print OUT "%                      at constant temperature                  \n%\n%\n";
            print OUT "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
            print OUT "mkdir KINETICS/DATA/$exp; mkdir KINETICS/PLOTS/$exp;\n";
            print OUT "\nclearvars;\n    sol=[];\n";
    	}elsif ($exp eq "processes") {
            print OUT "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n%\n%\n";
            print OUT "%                    This script calculates the                 \n";
            print OUT "%                ARRHENIUS PARAMETERS and CONSTANTS             \n%\n%\n";
            print OUT "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
            print OUT "clearvars;\n";
         };
            print OUT "h=6.62607588705515e-34; kb=1.380658045430573e-23; c=29979299840.00; R=8.31451034545898;\n";
            print OUT "Av=6.022139922973909e23; Fa=96485.296875000; toeV=1.60217648740e-19; electron=-1.60217648740e-19;\n\n";
	    print OUT "syms T V pH ;\n";
       close OUT;
     return ();
   }; #--> sub Introduction
#==============================================================================================================================
   sub process_sub {
           ($proc,$stoich)=@_; 
	   @tmp=split(/\s+/,$proc); $type=$tmp[0]; delete (@tmp[0]);
      foreach $t (@tmp) { if ($t) { push(@Nline,$t); }; }; @tmp=@Nline; @Nline=();
#------------------------------------------------------------------------------------------------------------------ components
           $r=0; $p=0; $p2=0; $step=''; $step2=''; $cha=''; @PR=(); @PTS=(); @PP=();
      foreach $cha (@tmp) { if ($cha) {
         if (($step ne '>') and ($step ne '>>')) { $step=$cha;
            if (($cha ne '+') and ($cha ne '>>') and ($cha ne '>')) { push(@PR,$cha); $r++; };
         }elsif(($step eq '>') and ($step2 ne '>')) { $step2=$cha;      ############# in case of TS for reactions
            if (($cha ne '+') and ($cha ne '>')) { push(@PTS,$cha); $ts++; };
         }elsif(($step eq '>>') or ($step2 eq '>')) {
	    if ($cha ne '+') { push(@PP,$cha); $p++; }; }; }; }; # foreach cha
#----------------------------------------------------------------------------------------------------------------- stoichiometry
           $r=0; $p=0;  $step=''; @tmp=split(/\s+/,$stoich);
      foreach $t (@tmp) { if ($t) { push(@Nline,$t); }; }; @tmp=@Nline; @Nline=();
      foreach $cha (@tmp) { if ($cha) {
         if (($step ne '>') and ($step ne '>>')) {  $step=$cha;
            if (($cha ne '+') and ($cha ne '>') and ($step ne '>>')){ $stoichio{@PR[$r]}->[$pr]=$cha; $r++; }; }else{
            if (($cha ne '+') and ($cha ne '>') and ($step ne '>>')){ $stoichio{@PP[$p]}->[$pr]=$cha; $p++; }; }; }; };
      open OUT2, ">>Processes";
         print OUT2 "Process$pr: $type\t"; $processWidth=0;
        foreach $r (@PR) { printf OUT2 "%3s*%s ",$stoichio{$r}->[$pr],$r; $processWidth=$processWidth+max(map length, $r)+3; };
     	 print OUT2 colored(">","bold red");
        foreach $ts (@PTS) { printf OUT2 colored(" $ts ","bright_yellow"); $processWidth=$processWidth+max(map length, $ts)+3; };
         print OUT2 colored(">","bold green");
        foreach $p (@PP) { printf OUT2 " $stoichio{$p}->[$pr]*%-s",$p; $processWidth=$processWidth+max(map length, $p)+3; };
	 printf OUT2 "%-*s\t",40-($processWidth);
	  $rE=0; $tE=(); $pE=0; 
        foreach $r (@PR) { printf OUT2 "%.2f ",$e0{$r}*$stoichio{$r}->[$pr]; $rE=$rE+$e0{$r}*$stoichio{$r}->[$pr]; }; print OUT2 colored("> ","bold red");
        foreach $ts (@PTS) { printf OUT2 "%.2f ",$e0{$ts}; $tE=$tE+$e0{$ts};}; print OUT2 colored("> ","bold green");
        foreach $p (@PP) { printf OUT2 "%.2f ",$e0{$p}*$stoichio{$p}->[$pr]; $pE=$pE+$e0{$p}*$stoichio{$p}->[$pr]; }; 
	$nS=$#PR+$#PTS+$#PP+3;
	if ($nS eq 6) { printf OUT2 " "; }elsif ($nS eq 5) { printf OUT2 "\t%5s";}elsif ($nS eq 4) { printf OUT2 "\t\t%5s";
	}elsif ($nS eq 3) { printf OUT2 "\t\t\t%5s"; }elsif ($nS eq 2) { printf OUT2 "\t\t\t\t%5s"; };
        if ($tE) { print OUT2 "Ea=";
		if (($tE-$rE > 3) or ($tE-$rE < 0)) { print OUT2 colored( sprintf("%.2f eV ",$tE-$rE),"bold red on_white");
	       	}else{ print OUT2 colored( sprintf("%.2f eV ",$tE-$rE),"bright_blue"); };
       	  }else{ printf OUT2 "%11s";};
	 print OUT2 "Er=";
	if (($pE-$rE < -3) or ($pE-$rE > 3)) { print OUT2 colored( sprintf("%.2f eV\n",$pE-$rE),"bold red on_white"); 
	  }else{ print OUT2 colored( sprintf("%.2f eV\n",$pE-$rE),"bright_blue"); };
      close OUT2;
      open OUT3, ">>Processes.txt";
         print OUT3 "Process$pr: $type\t"; $processWidth=0;
        foreach $r (@PR) { 	printf OUT3 "%3s*%s ",$stoichio{$r}->[$pr],$r;	$processWidth=$processWidth+max(map length, $r)+3; };	print OUT3 ">";
        foreach $ts (@PTS) { 	printf OUT3 " $ts "; 				$processWidth=$processWidth+max(map length, $ts)+3; };	print OUT3 ">";
        foreach $p (@PP) { 	printf OUT3 " $stoichio{$p}->[$pr]*%-s",$p; 	$processWidth=$processWidth+max(map length, $p)+3; };
         printf OUT3 "%-*s\t",40-($processWidth);
          $rE=0; $tE=(); $pE=0;
        foreach $r (@PR) { 	printf OUT3 "%.2f ",$e0{$r}*$stoichio{$r}->[$pr]; 	$rE=$rE+$e0{$r}*$stoichio{$r}->[$pr]; };	print OUT3 "> ";
        foreach $ts (@PTS) { 	printf OUT3 "%.2f ",$e0{$ts}; 				$tE=$tE+$e0{$ts}; }; 				print OUT3 "> ";
        foreach $p (@PP) { 	printf OUT3 "%.2f ",$e0{$p}*$stoichio{$p}->[$pr]; 	$pE=$pE+$e0{$p}*$stoichio{$p}->[$pr]; };
        $nS=$#PR+$#PTS+$#PP+3;
        if ($nS eq 6) { print OUT3 " "; }elsif ($nS eq 5) { printf OUT3 "\t%5s"; }elsif ($nS eq 4) { printf OUT3 "\t\t%5s";
         }elsif ($nS eq 3) { printf OUT3 "\t\t\t%5s"; }elsif ($nS eq 2) { printf OUT3 "\t\t\t\t%5s"; };
        if ($tE) { printf OUT3 "Ea=%.2f eV ",$tE-$rE; }else{ printf OUT3 "%11s";};
         printf OUT3 "Er=%.2f eV\n",$pE-$rE;
      close OUT3;
   return (\@PR, \@PTS, \@PP, $type);
   }; #--> sub process 
#==============================================================================================================================
  sub Import_sub {
         ($sys,$file)=@_;
         @done_import=(); $tmp=1;
     foreach $did (@done_import) { $tmp=1; if ($sys eq $did) { $tmp=0; };};
     if ($tmp eq 1) { push(@done_import,$sys);
        open OUT, ">>$file.m";
            print OUT "ENERGY$sys=readtable(\'./THERMODYNAMICS/DATA/$sys/G$sys.dat\');\n";
	    print OUT "PARTITION3D$sys=readtable(\'./THERMODYNAMICS/DATA/$sys/Q3D$sys.dat\');\n";
#	    print OUT "    IR$sys=readtable(\'./IRs/originals/$sys/intensities/IRSPECTRA\');\n";
	   foreach $mol (@molecules) { if ($sys eq $mol) {
	       print OUT "q3Dnotrans$sys=readtable(\'./THERMODYNAMICS/DATA/$sys/Q3Dnotrans$sys.dat\');\n";
	       print OUT "qt$sys=readtable(\'./THERMODYNAMICS/DATA/$sys/qtrans2D$sys.dat\');\n";
	       print OUT "qv$sys=readtable(\'./THERMODYNAMICS/DATA/$sys/qvib2D$sys.dat\');\n";
	       print OUT "ZPE2D$sys=readtable(\'./THERMODYNAMICS/DATA/$sys/ZPE2D$sys.dat\');\n";
               print OUT "ZPE$sys=readtable(\'./THERMODYNAMICS/DATA/$sys/ZPE$sys.dat\');\n";
               print OUT "Mass$sys=@Tmass{$sys}; area@sitetype{$sys}=@Acat{@sitetype{$sys}}; Nsites$sys=@nsitetype{$sys};\n"; };}; # if mol
               print OUT "\n"; 
	close OUT; }; # if tmp
     return();
   }; #--> sub Import
#==============================================================================================================================
   sub Interpolate_sub {
         ($sys)=@_;
        open OUT, ">>processes.m";
               print OUT "x$sys=[0, 0.25, 0.50, 0.75, 1];\n   Ey$sys=[ENERGY$sys, ENERGY$sys, ENERGY$sys, ENERGY$sys, ENERGY$sys];";
	       print OUT " Qy$sys=[PARTITION3D$sys, PARTITION3D$sys, PARTITION3D$sys, PARTITION3D$sys, PARTITION3D$sys];\n";
	       $E="Interp1(x$sys, Ey$sys, cov$sys, 'pchip','extrap')";
               $Q3D="Interp1(x$sys, Qy$sys, cov$sys, 'pchip', 'extrap')";
        close OUT; 
     return($E,$Q3D);
   }; #--> sub Interpolate
#==============================================================================================================================
   sub ProcessE_sub {
          ($typeP,$pr)=@_;
        open OUT, ">>processes.m"; print OUT "\n";
              @Etmp=(); @Esyms=();
           foreach $R (@PR) { if (@en{$R}) { push(@Etmp,"+stoichio$pr$R*@en{$R}");
		                      }else{ push(@Etmp,"+stoichio$pr$R*E$R"); push(@Esyms,"E$R"); }; };
           if (@PTS) { @Etmp2=(); 
              foreach $TS (@PTS) { @imafrq=split(/\s+/,@imafreq{$TS});           # Quantum tunnelling [J. Chem. Phys. 2006, 124, 044706.]
                   print OUT "   Qtunnel$TS=(kb*T/toeV)*log("; foreach $if (@imafrq) { print OUT "sinh($if*(h*c)/(2*kb*T))/($if*(h*c)/(2*kb*T))*"; }; print OUT "1);\n"; 
                 if (@en{$TS}) { push(@Etmp2,"+@en{$TS}+Qtunnel$TS");
		          }else{ push(@Etmp2,"+E$TS+Qtunnel$TS"); push(@Esyms,"E$TS"); }; };
               print OUT "syms"; foreach $s (@Esyms) { print OUT " $s"; }; print OUT "\n";
	       print OUT "   Ereactants$pr=0@Etmp;\n";
               print OUT "   Ets$pr=0@Etmp2;\n"; 
           }elsif (!@PTS) {
               if (($typeP eq 'A') or ($typeP eq 'a')) {
                    @Etmp2=();
                  foreach $R (@PR) { if (@en{$R}) { push(@Etmp2,"+stoichio$pr$R*@en{$R}"); }else{ push(@Etmp2,"+stoichio$pr$R*E$R"); push(@Esyms,"E$R");};
                     foreach $mol (@molecules) { if ($R eq $mol) { push(@Etmp2,"+stoichio$pr$R*(-Z$R+Z2D$R)"); push(@Esyms,"Z$R Z2D$R"); }; }; };
	             print OUT "syms"; foreach $s (@Esyms) { print OUT " $s"; }; print OUT "\n";
                     print OUT "   Ereactants$pr=0@Etmp;\n";
                     print OUT "   Ets$pr=0@Etmp2;\n"; 
	       }else{
	            @Etmp2=();
	          foreach $P (@PP) { if (@en{$P}) { push(@Etmp2,"+stoichio$pr$P*@en{$P}"); }else{ push(@Etmp2,"+stoichio$pr$P*E$P"); push(@Esyms,"E$P"); };
	             foreach $mol (@molecules) { if ($P eq $mol) { push(@Etmp2,"+stoichio$pr$P*(-Z$P+Z2D$P)"); push(@Esyms,"Z$P Z2D$P"); }; }; };
	             print OUT "syms"; foreach $s (@Esyms) { print OUT " $s"; }; print OUT "\n";
	             print OUT "   Ereactants$pr=0@Etmp;\n";
                     print OUT "   Ets$pr=0@Etmp2;\n";
	       }; # if A D or R
	   }; # if TS
	print OUT "AE$pr=Ets$pr-Ereactants$pr;";
        close OUT;
    return();
  }; #--> sub ProcessE
#==============================================================================================================================  
   sub ProcessQ_sub {
	  ($typeP,$pr)=@_;
           @Qtmp=(); @Qsyms=();
       if (($typeP eq 'A') or ($typeP eq 'a')) {
          foreach $R (@PR) { $go='no'; foreach $mol (@molecules) { if ($R eq $mol) { $go='yes'; }; };
             if ($go eq 'yes') { push(@Qtmp,"*(qtrans2D$R*Q3Dnotrans$R)^stoichio$pr$R"); push(@Qsyms,"qtrans2D$R Q3Dnotrans$R");
             }else{ if (@en{$R}) { push(@Qtmp,"*@q{$R}^stoichio$pr$R"); }else{ push(@Qtmp,"*Q3D$R^stoichio$pr$R"); push(@Qsyms,"Q3D$R"); }; }; };
       }else{ foreach $R (@PR) { if ($q{$R}) { push(@Qtmp,"*@q{$R}^stoichio$pr$R"); }else{ push(@Qtmp,"*Q3D$R^stoichio$pr$R"); push(@Qsyms,"Q3D$R");};};};
       if (@PTS) { @Qtmp2=();
          foreach $TS (@PTS) { if ($q{$TS}) { push(@Qtmp2,"*@q{$TS}"); }else{ push(@Qtmp2,"*Q3D$TS"); push(@Qsyms,"Q3D$TS"); };};
       }elsif (!@PTS) { @Qtmp2=();
          if (($typeP eq 'A') or ($typeP eq 'a')) {
             foreach $R (@PR) { $go='no'; foreach $mol (@molecules) { if ($R eq $mol) { $go='yes'; }; };
                if ($go eq 'yes') { push(@Qtmp2,"*qvib2D$R^stoichio$pr$R"); push(@Qsyms,"qvib2D$R");
                }else{ if (@q{$R}) { push(@Qtmp2,"@q{$R}^stoichio$pr$R"); }else{ push(@Qtmp2,"*Q3D$R^stoichio$pr$R"); push(@Qsyms,"Q3D$R"); }; }; };
          }elsif (($typeP eq 'D') or ($typeP eq 'd')) {
	      foreach $P (@PP) { $go='no'; foreach $mol (@molecules) { if ($P eq $mol) {  $go='yes'; }; };
	         if ($go eq 'yes') { push(@Qtmp2,"*qvib2D$P^stoichio$pr$P"); push(@Qsyms,"qvib2D$P");
	         }else{ if (@q{$P}) { push(@Qtmp2,"@q{$P}^stoichio$pr$P"); }else{ push(@Qtmp2,"*Q3D$P^stoichio$pr$P"); push(@Qsyms,"Q3D$P"); }; }; };
	  }elsif (($typeP eq 'R') or ($typeP eq 'r')) {
	      foreach $R (@PR) { $go='no'; foreach $mol (@molecules) { if ($R eq $mol) { $go='yes'; }; };
	         if ($go eq 'yes') { push(@Qtmp2,"*qvib2D$R^stoichio$pr$R"); push(@Qsyms,"qvib2D$R"); }; };                 
	      foreach $P (@PP) { $go='no'; foreach $mol (@molecules) { if ($P eq $mol) {  $go='yes'; }; };
	         if ($go eq 'yes') { push(@Qtmp2,"*qvib2D$P^stoichio$pr$P"); push(@Qsyms,"qvib2D$P");
	         }else{ if (@q{$P}) { push(@Qtmp2,"@q{$P}^stoichio$pr$P"); }else{ push(@Qtmp2,"*Q3D$P^stoichio$pr$P"); push(@Qsyms,"Q3D$P"); };};};};};
       open OUT, ">>processes.m";
           print OUT "syms"; foreach $s (@Qsyms) { print OUT " $s"; }; print OUT "\n";
           print OUT "   Qreactants$pr=1@Qtmp;\n";
           print OUT "   Qts$pr=1@Qtmp2;\n";
        close OUT;
     return();
   }; #--> sub ProcessQ
#==============================================================================================================================
   sub ProcessK_sub {
	  ($typeP,$pr)=@_;
          @r=(); $tmp=(); $nrate="Krate$pr";

#	  $h2o="n"; $h3op="n"; $ohn="n";
#	foreach $R (@PR) {	if ($R eq "H2O") { $h2o="R"; 	push(@r,"y(@y{$R})^2"); };
#				if (($R eq "H3O") or ($R eq "H3Op")) { $h3op="R"; 	push(@r,"y(@y{$R})^1");	};
#				if (($R eq "OH") or ($R eq "OHn")) { $ohn="R"; 	push(@r,"y(@y{$R})^1"); }; };
#       	foreach $P (@PP) { 	if ($P eq "H2O") { $h2o="P"; };
#                                if (($P eq "H3O") or ($P eq "H3Op")) { $h3op="P"; };
#                                if (($P eq "OH") or ($P eq "OHn")) { $ohn="P"; }; };

# Temperature dependence of the water ionization constant.
# Based on data from "Release on the Ionization Constant of H2O". International Association for Properties of Water and Steam. August 2007
# 			pK_w(T,P)= ( pK_w[0.1Pa]+pK_w[25Pa] )/2  --> 4th degree polynomial interpolation
        open OUT, ">>processes.m";
#	$mol="H2O";
#	if (($h2o eq "R") and ($h3op eq "P") and ($ohn eq "P")) {	print OUT "Arrhenius$pr=1;	P$mol=@pressure{$mol};\n";
#									print OUT "Krate$pr=(1/P$mol)*10^(6.8522e-10*T^4-1.305e-6*T^3+9.7231e-4*T^2-0.33802*T+57.515);\n";
#	}elsif (($h2o eq "P") and ($h3op eq "R") and ($ohn eq "R")) {	print OUT "Arrhenius$pr=1;      P$mol=@pressure{$mol};\n";
#                                                                        print OUT "Krate$pr=(1/P$mol)*10^-(6.8522e-10*T^4-1.305e-6*T^3+9.7231e-4*T^2-0.33802*T+57.515);\n";
#	}else{ @r=(); $tmp=(); 
		foreach $R (@PR) { foreach $mol (@molecules) { if ($R eq $mol) {push(@r,"y(@y{$mol})^stoichio$pr$mol"); $tmp=$mol;};};
                         foreach $sur (@surfaces) { if ($R eq $sur) { 		push(@r,"y(@y{$sur})^stoichio$pr$sur"); };};
                         foreach $cat (@catalysts) { if ($R eq $cat) { 		push(@r,"y(@y{$cat})^stoichio$pr$cat"); };};};
      		if (($typeP eq 'A') or ($typeP eq 'a')) { 
             		print OUT "sticky$pr=(Qts$pr/Qreactants$pr)*exp(-(AE$pr*toeV/(kb*T)));\n";
             		print OUT "Arrhenius$pr=area$sitetype{$tmp}*1/((2*pi*Mass$tmp*kb*T)^(1/2));\n";
             		print OUT "Krate$pr=Arrhenius$pr*sticky$pr;\n";
      		}else{ 	print OUT "Arrhenius$pr=(kb*T/h)*(Qts$pr/Qreactants$pr);\n";   
             		print OUT "Krate$pr=Arrhenius$pr*exp(-(AE$pr*toeV/(kb*T)));\n"; };
#	};

      foreach $rpr (@r) { $nrate="$nrate*$rpr"; }; 
       push(@rates,"rate$pr=$nrate");
        close OUT;
    return();
  }; #--> sub ProcessK
#==============================================================================================================================
   sub DRCrates_sub {
          ($exp)=@_;
        open OUT, ">>$exp.m"; %done=();
	print OUT "\nsyms "; for ($pr=1; $pr<=$#process; $pr++ ) { print OUT " Krate$pr"; }; print OUT "\n";
	print OUT "syms "; for ($pr=1; $pr<=$#process; $pr++ ) { @PR=split(/\s+/,@ProcessReactants[$pr]); 
                foreach $R (@PR) { if (!$done{$R}) {
			foreach $mol (@molecules) { if ($R eq $mol) {  print OUT " P$R"; };};
                   	foreach $sur (@surfaces) { if ($R eq $sur) {   print OUT " coverage$R"; };};
                   	foreach $cat (@catalysts) { if ($R eq $cat) {  print OUT " coverage$R"; };};
			$done{$R}=1; };};}; print OUT "\n";
	for ($pr=1; $pr<=$#process; $pr++ ) { @DRCr=(); $tmp=(); $DRCnrate="Krate$pr";
		@PR=split(/\s+/,@ProcessReactants[$pr]); 
      		foreach $R (@PR) {
			foreach $mol (@molecules) { if ($R eq $mol) {  push(@DRCr,"P$mol^stoichio$pr$mol"); $tmp=$mol;};};
                        foreach $sur (@surfaces) { if ($R eq $sur) {   push(@DRCr,"coverage$sur^stoichio$pr$sur"); };};
                        foreach $cat (@catalysts) { if ($R eq $cat) {  push(@DRCr,"coverage$cat^stoichio$pr$cat"); };};};
		foreach $rpr (@DRCr) { $DRCnrate="$DRCnrate*$rpr"; };
             	print OUT "   rate$pr=$DRCnrate;\n";
	}; # for pr
        close OUT;
    return();
  }; #--> sub DRCrates
#==============================================================================================================================
   sub ProcessParameters_sub {
	  ($typeP,$pr)=@_;
	  if ($ttemp) { $row=1+($ftemp-$itemp)/$ttemp; }else{ $row=1; };
  	  if ($sVext) { $row=$row+($nVext+$pVext)/$sVext; };
	  if ($spH) { $row=$row+($fpH-$ipH)/$spH; };
       open OUT, ">>processes.m";
                 print OUT " fileID=fopen(\"./KINETICS/PROCESS/ReactionParameters$pr.dat\",'a+');";
          if (($typeP eq 'A') or ($typeP eq 'a')) {
		 print OUT " fprintf(fileID, \'#  T \t sticky$pr \t\t Arrhenius$pr \t\t Krate$pr\\n\');\n";
            }else{ 
                 print OUT " fprintf(fileID, \'#  T \t Arrhenius$pr \t\t Krate$pr\\n\');\n"; };  print OUT "\n";
	                if ($ttemp) {   print OUT "j=1;\nfor T = $itemp:$ttemp:$ftemp\n";
        	        }else{          print OUT "j=1;\nT=$itemp;\n"; };
		$done="no";
	     foreach $R (@PR) { @do{$R}="yes"; }; foreach $P (@PP) { @do{$P}="yes"; };  foreach $TS (@PTS) { @do{$TS}="yes"; };
         foreach $R (@PR) { if ($done eq "no") { if ($ttemp) {
				     	    print OUT "   while ENERGY$R\{j,1} ~= T ; j=j+1; end\n"; }; $done="yes"; };};

         foreach $R (@PR) { if (@do{$R} eq "yes") { $qmol="no";
		        foreach $mol (@molecules) { if ($R eq $mol) { $qmol="yes" };};
                if ($qmol eq "yes") {   print OUT "      E$R=ENERGY$R\{j,2}; Z$R=ZPE$R\{j,2}; Z2D$R=ZPE2D$R\{j,2};\n";
			                            print OUT "      Q3D$R=PARTITION3D$R\{j,2}; Q3Dnotrans$R=q3Dnotrans$R\{j,2};";
				                        print OUT " qtrans2D$R=qt$R\{j,2}; qvib2D$R=qv$R\{j,2};\n"; @do{$R}="no";
		                        }else{  print OUT "      E$R=ENERGY$R\{j,2}; Q3D$R=PARTITION3D$R\{j,2};\n"; @do{$R}="no"; };};}; #foreach PR
         foreach $TS (@PTS) { if (@do{$TS} eq "yes") {
                                        print OUT "      E$TS=ENERGY$TS\{j,2}; Q3D$TS=PARTITION3D$TS\{j,2};\n"; @do{$TS}="no";};}; #foreach PTS
#         foreach $P (@PP) { if (@do{$P} eq "yes") { $qmol="no";
#		     foreach $mol (@molecules) { if ($P eq $mol) { $qmol="yes" };};
#             if ($qmol eq "yes") {      print OUT "      E$P=ENERGY$P\{j,2}; Z$P=ZPE$P\{j,2}; Z2D$P=ZPE2D$P\{j,2};\n";
#    		                            print OUT "      Q3D$P=PARTITION3D$P\{j,2}; Q3Dnotrans$P=q3Dnotrans$P\{j,2};";
#            				            print OUT " qtrans2D$P=qt$P\{j,2}; qvib2D$P=qv$P\{j,2};\n"; @do{$P}="no";
#    		                    }else{  print OUT "      E$P=ENERGY$P\{j,2}; Q3D$P=PARTITION3D$P\{j,2};\n"; @do{$P}="no"; };};}; #foreach PR
		if (($typeP eq 'A') or ($typeP eq 'a')) {
                 print OUT " fprintf(fileID, '%.4f %1.15E %1.15E %1.15E\\n', T, subs(sticky$pr), subs(Arrhenius$pr), subs(Krate$pr));\n";
             }else{ 
                 print OUT " fprintf(fileID, '%.4f %1.15E %1.15E\\n', T, subs(Arrhenius$pr), subs(Krate$pr));\n"; };
        if ($ttemp) {  print OUT "end\nfclose(fileID);\n";
        	}else{ print OUT "fclose(fileID);\n"; };
	    print OUT "\n\n";
        close OUT;
    return();
  }; #--> sub ProcessParameters
#==============================================================================================================================
   sub variables_sub {
          ($exp)=@_;
          @variables=();  %tmpPressure=(); %tmpCov=(); 
	open OUT,">>$exp.m";  
	if (($exp eq "const_TEMP") or ($exp eq "variable_TEMP")) {
	 	print OUT "fileID=fopen(\"./KINETICS/DATA/$exp/solution$exp.dat\",'a+');\n";
         	print OUT " fprintf(fileID, \'#%2s %6s %6s %6s";
		foreach $mol (@molecules) { if (@y{$mol}) { $do="yes"; foreach $t (@transitions) { if ($t eq $mol) { $do='no'; };};
					if ($do eq "yes" ) { print OUT " %18s"; };};};
		foreach $rs (@react_species) {
        		if (@y{$rs}) { $do='yes'; foreach $mol (@molecules) { if ($mol eq $rs) { $do='no'; };};
						  foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                				  foreach $sur (@surfaces) { if ($sur eq $rs) { $do='no'; };};
               	 		if ($do eq 'yes') { print OUT " %20s";  };};};
		foreach $sur (@surfaces) { if (@y{$sur}) { print OUT " %18s";  };};
		 print OUT "\\n\',\'V\',\'pH\',\'T\',\'t\'";
        	foreach $mol (@molecules) { if (@y{$mol}) { $do="yes"; foreach $t (@transitions) { if ($t eq $mol) { $do='no'; };};
                                        if ($do eq "yes" ) { print OUT ",\'$mol\'"; };};};
        	foreach $rs (@react_species) {
                	if (@y{$rs}) { $do='yes'; foreach $mol (@molecules) { if ($mol eq $rs) { $do='no'; };};
						  foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                		  foreach $sur (@surfaces) { if ($sur eq $rs) {$do='no'; };};
                	        if ($do eq 'yes') { print OUT ",\'$rs\'";  };};};
        	foreach $sur (@surfaces) { if (@y{$sur}) { print OUT ",\'$sur\'";  };}; print OUT ");\n\n";
	}elsif ($exp eq "TPR") {
		print OUT "fileID=fopen(\"./KINETICS/DATA/$exp/TPR_$TPRmolecule.dat\",'a+');\n";
                print OUT " fprintf(fileID, \'#%2s %6s %6s %6s";
                foreach $mol (@molecules) { if (@y{$mol}) { $do="yes"; foreach $t (@transitions) { if ($t eq $mol) { $do='no'; };};
                                        if ($do eq "yes" ) { print OUT " %18s"; };};};
                foreach $rs (@react_species) {
                        if (@y{$rs}) { $do='yes'; foreach $mol (@molecules) { if ($mol eq $rs) { $do='no'; };};
						  foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                		  foreach $sur (@surfaces) { if ($sur eq $rs) {$do='no'; };};
                                if ($do eq 'yes') { print OUT " %20s";  };};};
                foreach $sur (@surfaces) { if (@y{$sur}) { print OUT " %18s";  };};
                 print OUT "\\n\',\'V\',\'pH\',\'T\',\'t\'";
                foreach $mol (@molecules) { if (@y{$mol}) { $do="yes"; foreach $t (@transitions) { if ($t eq $mol) { $do='no'; };};
                                        if ($do eq "yes" ) { print OUT ",\'$mol\'"; };};};
                foreach $rs (@react_species) {
                        if (@y{$rs}) { $do='yes'; foreach $mol (@molecules) { if ($mol eq $rs) { $do='no'; };};
						  foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                		  foreach $sur (@surfaces) { if ($sur eq $rs) {$do='no'; };};
                                if ($do eq 'yes') { print OUT ",\'$rs\'";  };};};
                foreach $sur (@surfaces) { if (@y{$sur}) { print OUT ",\'$sur\'";  };}; print OUT ");\n\n";
        }elsif ($exp eq "RateControl") {

	}; # if exp
	if ($exp eq "TPR") { 
	     	print OUT " V=0.00;\n pH=0.00;\n";
	        $TPRdone="n";	
		for ($pr=1; $pr<=$#process; $pr++ ) { 
			if ((@typeproc[$pr] eq "A") or (@typeproc[$pr] eq "a")) {
               			@PR=split(/\s+/,@ProcessReactants[$pr]);
		                foreach $R (@PR) { if ($R eq $TPRmolecule) { @PP=split(/\s+/,@ProcessProducts[$pr]);
					foreach $P (@PP) { $go="y"; 
						foreach $P (@PP) { $go="y";
                        				foreach $sur (@surfaces) { if ($P eq $sur) { $go="n"; };};
                        				foreach $mol (@molecules) {if ($P eq $mol) { $go="n"; };};
							if ($TPRdone eq "y") { $go="n"; };
                        			if ($go eq "y") { print OUT "for coverage$P = 0.1:0.3:1.0\n";
						       push(@variables,$P); $TPRdone="y"; };};};};};};};

                ($IniCon)=&InitialConcentrations_sub($exp); @IC="@$IniCon"; open OUT,">>$exp.m";
	        if ($ttemp) {   print OUT "\nfor T = $itemp:$ttemp:$ftemp\n";    #################################### there was a i=1; before T loop
                }else{          print OUT "\nT=$itemp;\n"; };
	}elsif ($exp eq "const_TEMP") {
                if (($nVext) and ($pVext) and ($sVext) and ($nVext != $pVext)) {
			print OUT "fileID2=fopen(\"./KINETICS/DATA/$exp/voltammogram.dat\",'a+');\n";
         		print OUT " fprintf(fileID2, \'# %3s %7s %12s %22s\\n\',\'T\',\'(V/s)\',\'Voltage (V)\',\'Intensity (A/m^2)\');\n\n";
                        print OUT "for V = $nVext:$sVext:$pVext\n";
                }elsif(($nVext) and ((!$pVext) or (!$sVext) or ($nVext != $pVext))) {
			print OUT "fileID2=fopen(\"./KINETICS/DATA/$exp/voltammogram.dat\",'a+');\n";
                        print OUT " fprintf(fileID2, \'# %3s %7s %12s %22s\\n\',\'T\',\'(V/s)\',\'Voltage (V)\',\'Intensity (A/m^2)\');\n\n";
                        print OUT " V=$nVext;\n";
                }else{  print OUT " V=0.00;\n"; };
                if (($ipH) and ($fpH) and ($spH) and ($ipH != $fpH)) {
                        print OUT "for pH = $ipH:$spH:$fpH\n";  @tmpPressure{"H3Op"}="10^-pH"; @tmpPressure{"OHn"}="10^-(14-pH)";
                }elsif(($ipH) and ((!$fpH) or (!$fpH) or ($ipH != $fpH))) {
                        print OUT " pH=$ipH;\n";  @tmpPressure{"H3Op"}="10^-pH"; @tmpPressure{"OHn"}="10^-(14-pH)";
                }else{ 	print OUT " pH=0.00;\n"; };
                if ($ttemp) {   print OUT "\nfor T = $itemp:$ttemp:$ftemp\n";
                }else{          print OUT "\nT=$itemp;\n"; };
            	foreach $mol (@molecules) { $go="y"; foreach $t (@transitions) { if ($t eq $mol) { $go="n"; };};
			if ($go eq "y") {
				if (($ipressure{$mol}) and ($fpressure{$mol})) { push(@variables,$mol);
               				print OUT "for P$mol = $ipressure{$mol}:$spressure{$mol}:$fpressure{$mol}\n"; };};};
            	foreach $cat (@catalysts) { if (($icov{$cat}) and ($fcov{$cat})) { push(@variables,$cat);
               		print OUT "for coverage$cat = $icov{$cat}:$scov{$cat}:$fcov{$cat}\n"; }; };
                ($IniCon)=&InitialConcentrations_sub($exp); @IC="@$IniCon"; open OUT,">>$exp.m"; 
         }elsif ($exp eq "variable_TEMP") {
                if (($nVext) and ($pVext) and ($sVext) and ($nVext != $pVext)) {
                        print OUT "for V = $nVext:$sVext:$pVext\n";
                }elsif(($nVext) and ((!$pVext) or (!$sVext) or ($nVext != $pVext))) {
                        print OUT " V=$nVext;\n";
                }else{  print OUT " V=0.00;\n"; };
                if (($ipH) and ($fpH) and ($spH) and ($ipH != $fpH)) {
                        print OUT "for pH = $ipH:$spH:$fpH\n";  @tmpPressure{"H3Op"}="10^-pH";  @tmpPressure{"OHn"}="10^-(14-pH)";
                }elsif(($ipH) and ((!$fpH) or (!$fpH) or ($ipH != $fpH))) {
                        print OUT " pH=$ipH;\n";   @tmpPressure{"H3Op"}="10^-pH";  @tmpPressure{"OHn"}="10^-(14-pH)";
                }else{ 	print OUT " pH=0.00;\n"; };
            	foreach $mol (@molecules) { $go="y"; foreach $t (@transitions) { if ($t eq $mol) { $go="n"; };};
                        if ($go eq "y") {
				if (($ipressure{$mol}) and ($fpressure{$mol})) { push(@variables,$mol);
               			print OUT "for P$mol = $ipressure{$mol}:$spressure{$mol}:$fpressure{$mol}\n"; };};};
            	foreach $cat (@catalysts) { if (($icov{$cat}) and ($fcov{$cat})) { push(@variables,$cat);
               		print OUT "for coverage$cat = $icov{$cat}:$scov{$cat}:$fcov{$cat}\n"; }; };
                ($IniCon)=&InitialConcentrations_sub($exp); @IC="@$IniCon"; open OUT,">>$exp.m";
                if ($ttemp) {   print OUT "\nfor T = $itemp:$ttemp:$ftemp\n";
                }else{          print OUT "\nT=$itemp;\n"; };
         }elsif ($exp eq "RateControl"){print OUT "\n V=0.00;\n"; 
                if ($ipH) {		print OUT " pH = $ipH;\n"; @tmpPressure{"H3Op"}="10^-pH"; @tmpPressure{"OHn"}="10^-(14-pH)";
                }else{  		print OUT " pH=0.00;\n"; };		 
                if ($itemp) {   	print OUT " T = $itemp;\n"; };
                foreach $mol (@molecules) { $go="y"; foreach $t (@transitions) { if ($t eq $mol) { $go="n"; };};
                        if ($go eq "y") {
				 if (($ipressure{$mol}) and ($fpressure{$mol})) { push(@variables,$mol);
                        	print OUT " P$mol=0.5;\n"; };};};
                foreach $cat (@catalysts) { if (($icov{$cat}) and ($fcov{$cat})) { push(@variables,$cat);
                        print OUT " coverage$cat=0.5\n"; }; };
                print OUT "\nfile1=sprintf(\"./KINETICS/DATA/$exp/DRC_%G.dat\",T);\n";
 		print OUT "file2=sprintf(\"./KINETICS/DATA/$exp/DSC_%G.dat\",T);\n";
                print OUT "DRCfile=fopen(file1,'a+');\n";
	 	print OUT "DSCfile=fopen(file2,'a+');\n\n";
                ($IniCon)=&InitialConcentrations_sub($exp); @IC="@$IniCon"; open OUT,">>$exp.m";
	}; # if exp
        close OUT;
    return(\@IC);
   }; #--> sub variables
#==============================================================================================================================
   sub InitialConcentrations_sub {
           ($exp)=@_;
           @IniCon=(); %tmp=(); 
       open OUT, ">>$exp.m";
         foreach $mol (@molecules) { $i="yes"; 	foreach $v (@variables) { if ($mol eq $v) { $i="no"; };};
 		                     $go="no"; 	foreach $rs (@react_species) { if ($rs eq $mol) { $go="yes"; };};
				     $tr="no";	foreach $t (@transitions) { if ($mol eq $t) { $tr="yes"; };};
              if (($i eq "yes") and ($go eq "yes") and ($tr eq "no")) { 
		      	push(@IniCon,"P$mol ");
			if (@tmpPressure{$mol}) { print OUT "P$mol=@tmpPressure{$mol}; ";
			}else{	if ($exp eq "TPR") { print OUT "P$mol=0.0; ";
				}else{	print OUT "P$mol=@pressure{$mol}; ";};};};};

         foreach $rs (@react_species) { $go="yes"; $tr="no"; $m="no";
	     foreach $mol (@molecules) { if ($rs eq $mol) { $go="no"; $m="yes"; };}; 
	     foreach $t (@transitions) { if ($rs eq $t) { $go="no"; $tr="yes"; };};
	     foreach $sur (@surfaces) { if ($rs eq $sur) { $go="no"; 
					}elsif ((@sitetype{$sur} eq @sitetype{$rs}) and ($tr eq "no") and ($m eq "no")) {
		       				@tmp{$sur}="@tmp{$sur} coverage$rs"; };};
             if ($go eq "yes") { $i="yes";  
# done with python	if (@freq{$rs}) { &localIR_sub($rs,@ipath{$rs}); };
			foreach $v (@variables) { if ($rs eq $v) { $i="no"; };}; 
             		if ($i eq "yes") {	  if ($tmpCov{$rs}) { print OUT "coverage$rs=@tmpCov{$rs}; ";
		  				  }else{ print OUT "coverage$rs=@coverage{$rs}; ";};
			push(@IniCon,"coverage$rs "); };};};
         foreach $sur (@surfaces) { $tmp2=(); @tmp3=split(/\s+/,@tmp{$sur});
	    foreach $c (@tmp3) { if (!$tmp2) { $tmp2="$c"; }else{ $tmp2="$tmp2+$c"; };};
            print OUT "coverage$sur=1-($tmp2);\n"; push(@IniCon,"coverage$sur "); };
       close OUT;
    return(\@IniCon);
  }; #--> sub InitialConcentrations
#==============================================================================================================================
   sub ODE_import_sub {
       ($exp)=@_;
       open OUT, ">>$exp.m";
         for ($pr=1; $pr<=$#process; $pr++ ) { print OUT "process$pr=readtable(\'./KINETICS/PROCESS/ReactionParameters$pr.dat\');\n"; }; # for process			
         print OUT "\n";
       close OUT;
    return();
   }; #--> sub ODE_import
#==============================================================================================================================
   sub ODE_call_sub {
        ($exp)=@_;
      open OUT, ">>$exp.m";
           @transfer=();
# Alberto 18/01/2019
 	   print OUT "\nj=1;\n";
	if ($ttemp) { print OUT "   while (process1\{j,1} ~= T)\n      j=j+1;\n   end\n"; };
        for ($pr=1; $pr<=$#process; $pr++ ) {
# 	Phys. Rev. Lett. 2007, 99, 126101                             DOI:https://doi.org/10.1103/PhysRevLett.99.126101
# 	J. Phys. Chem. C, 2010, 114 (42), pp 18182–18197              DOI: 10.1021/jp1048887
#	The hydrogen coverage will be dependent on the potential via the reaction:
#       		        H+ + e- + *  --> H*
# 	At standard conditions (298 K, pH 0, 1 bar H2) and U = 0 V vs NHE,
# 		the left-hand side is in equilibrium with hydrogen gas.
# 	At finite bias, U, the chemical potential of the electron will be linearly dependent on the bias.
# 	The reaction free energy can be written as:
#       		        ΔGH* = AG + AG(U)= AG + −eU 
# 	defines the chemical potential of H*.
	if ($exp ne "TPR") {
		if ($iVext) {
			if ($ineexch[$pr]) {
				print OUT "        AEv$pr=@ineexch[$pr]*V;  ";
			}else{ 	@PR=split(/\s+/,@ProcessReactants[$pr]); $H="no";
                             	foreach $l (@PR) { if ($l) { push(@Nline,$l); }; }; @PR=@Nline; @Nline=();
                            	foreach $R (@PR) { if ($R eq "H") { $H="R"; $a=$R; };};
                            	@PP=split(/\s+/,@ProcessProducts[$pr]);
                            	foreach $l (@PP) { if ($l) { push(@Nline,$l); }; }; @PP=@Nline; @Nline=();
                            	foreach $P (@PP) { if ($P eq "H") { $H="P"; $a=$P; };};
				if ($H eq "R") {    	print OUT "        AEv$pr=-$stoichio{$a}->[$pr]*V;  \n";
	                        }elsif ($H eq "P") {	print OUT "        AEv$pr=$stoichio{$a}->[$pr]*V;  \n"; }; };
                }else{  print OUT "        AEv$pr=0;  "; };

# 	J. Phys. Chem. B, 2004, 108 (46), pp 17886–17892    DOI: 10.1021/jp047349j
# 	At a pH different from 0, we can correct the free energy of H+ ions by the concentration dependence of the entropy:
#       		        G = H -TS + kT ln(Products/Reactants)   ;    pH = -log[H3O+]
#		               G(pH) = −kT ln[H+]= kT ln (10) × pH.

                if ($ipH) { @PR=split(/\s+/,@ProcessReactants[$pr]); $H3O="no";
       	                    foreach $l (@PR) { if ($l) { push(@Nline,$l); }; }; @PR=@Nline; @Nline=();
               	            foreach $R (@PR) { if (($R eq "H3O") or ($R eq "H3Op")) { $H3O="R"; };};
                       	    @PP=split(/\s+/,@ProcessProducts[$pr]);
                            foreach $l (@PP) { if ($l) { push(@Nline,$l); }; }; @PP=@Nline; @Nline=();
       	                    foreach $P (@PP) { if (($P eq "H3O") or ($P eq "H3Op")) { $H3O="P"; };};
               	        if ($H3O eq "R") {         print OUT " AEph$pr=(kb*T*log(10)*pH)/toeV;\n";
                       	    }elsif ($H3O eq "P") { print OUT " AEph$pr=-(kb*T*log(10)*pH)/toeV;\n";
                            }else{                 print OUT " AEph$pr=0;\n"; };
       	        }else{ print OUT " AEph$pr=0;\n"; };
           if (($typeproc[$pr] eq "A") or ($typeproc[$pr] eq "a")) {
		  print OUT "     Krate$pr=process$pr\{j,2}*process$pr\{j,4}*exp(-((AEv$pr+AEph$pr)*toeV)/(kb*T));\n"; push(@transfer,"Krate$pr"); 
           }else{ print OUT "     Krate$pr=process$pr\{j,3}*exp(-((AEv$pr+AEph$pr)*toeV)/(kb*T));\n"; push(@transfer,"Krate$pr"); };
	}else{	# TPR?
           if (($typeproc[$pr] eq "A") or ($typeproc[$pr] eq "a")) {
                  print OUT "     Krate$pr=process$pr\{j,2}*process$pr\{j,4};\n"; push(@transfer,"Krate$pr");
           }else{ print OUT "     Krate$pr=process$pr\{j,3};\n"; push(@transfer,"Krate$pr"); };	
	}; # TPR?
	}; # for process
#  	Angewandte Chemie (2016) Volume 128, Issue 26, Page 7627, DOI: 10.1002/ange.201511804
# 	 current intensity=[[electron charge * number of electrons transferred * area of 
# 	 active site (m^2) * reaction rates]] x all the processes exchanging electrons
# 	  units = Ampere/m^2
	if ($exp eq "const_TEMP") { print OUT "     RedoxIntensity=";   
	        for ($pr=1; $pr<=$#process; $pr++ ) {
			if (@ineexch[$pr] != 0) {
				@PR=split(/\s+/,@ProcessReactants[$pr]);
	                        foreach $l (@PR) { if ($l) { push(@Nline,$l); }; }; @PR=@Nline; @Nline=();
				print OUT "(@ineexch[$pr]*Fa*@Acat{@sitetype{$PR[0]}}*Krate$pr/Av)+"; };
		}; print OUT "0;\n\n";
	}; # if const_TEMP

        if (!@Rspecies) {
       		foreach $mol (@molecules) { if (@y{$mol}) { $go="yes";
			foreach $t (@transitions) { if ($mol eq $t) { $go="no"; };};
		        if ($go eq "yes") {push(@Rspecies,"P$mol"); };}};
       		foreach $rs (@react_species) {
	   		if (@y{$rs}) { $do='yes'; foreach $mol (@molecules) { if ($mol eq $rs) { $do='no'; };};
		        			  foreach $sur (@surfaces) { if ($sur eq $rs) { $do='no'; };};
                                                  foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
	      		if ($do eq 'yes') { push(@Rspecies,"coverage$rs");  };};};
       		foreach $sur (@surfaces) { if (@y{$sur}) { push(@Rspecies,"coverage$sur"); };};
	}; # if Rspecies 

          print OUT "\n   tspan=[$itime $ftime];\n   IC=[@Rspecies];\n";
# alberto 09/2019
	if ($exp ne "TPR") { $ode="myode";  
  	  	printf OUT "   options=odeset('Refine',5,'NonNegative',(1:%d),'RelTol',1e-7,'AbsTol',1e-5);\n",$#Rspecies+1;
	}else{ $ode="mytpr";
       		printf OUT "   options=odeset('Refine',5,'NonNegative',(1:%d),'RelTol',1e-7,'AbsTol',1e-5);\n",$#Rspecies+1; };

          print OUT "\nsolution=ode15s(\@(t,y) $ode(t,y,";
       foreach $trn (@transfer) { if ($trn ne @transfer[$#transfer]) { print OUT "$trn,"; }else{ print OUT "$trn),tspan,IC,options);\n\n"; };};
     close OUT;  
    return();
      }; #--> sub ODE_call
#==============================================================================================================================
   sub ODE_solution_sub {
 	($exp)=@_;     
	 $timeSteps=int(($ftime-$itime)/$ttime);  
      open OUT, ">>$exp.m";
#             print OUT "   variables=zeros($timeSteps,3); variables(:,1)=V; variables(:,2)=pH; variables(:,3)=T;\n";
	      print OUT "   concentrations=transpose(deval(solution,transpose(linspace($itime,$ftime,$timeSteps))));\n";
#             print OUT "   sol=vertcat(sol,[variables,transpose(linspace($itime,$ftime,$timeSteps)),concentrations]);\n\n";
             print OUT "   results=[transpose(linspace($itime,$ftime,$timeSteps)),concentrations];\n\n";
	if ($exp ne "const_TEMP") {    
	$n=1;
	foreach $mol (@molecules) { foreach $rs (@react_species) {
		if ($mol eq $rs) {  $go="y"; foreach $t (@transitions) { if ($t eq $mol) { $go="n"; };};
                        if ($go eq "y") { print OUT "     P$mol=concentrations(end,$n); "; $n++; };};};}; print OUT "\n"; 
        foreach $rs (@react_species) { $do='yes'; foreach $mol (@molecules) { if ($mol eq $rs) { $do='no'; };};
                                                  foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
    		if ($do eq 'yes') { print OUT "     coverage$rs=concentrations(end,$n); "; $n++; };}; print OUT "\n";
	foreach $sur (@surfaces) { if ($y{$sur}) { print OUT "     coverage$sur=concentrations(end,$n); "; $n++; };}; print OUT "\n";
         }; # if const_TEMP
	close OUT;
     return();
       }; #--> sub ODE_solution


#==============================================================================================================================
   sub ODE_closingloops_sub {
          ($exp)=@_;
      open OUT, ">>$exp.m";
      if ($exp eq "TPR") {
                print OUT "\n   fprintf(fileID,\' %.2f  %.2f  %.2f  %.6f";
                foreach $Rs (@Rspecies) { print OUT "  %.15E"; };
                                          print OUT "\\n\',V,pH,T,results(end,:));\n";
                if ($ttemp) { print OUT " end % Temperature\n"; };
		print OUT " end % coverage from adsorbed $TPRmolecule\n";
      }elsif ($exp eq "const_TEMP") { 
                print OUT "\n   for i = 1:$timeSteps\n      fprintf(fileID,\' %.2f  %.2f  %.2f  %.6f";
                foreach $Rs (@Rspecies) { print OUT "  %.15E"; };
                                          print OUT "\\n\',V,pH,T,results(i,:));\n   end\n";
                if ($nVext) { print OUT "   if T == 300\n     fprintf(fileID2, \' %.3f %.3f %.6f %.15E\\n\',T, abs(subs(V/$ftime)),V,RedoxIntensity);\n   end\n\n"; };
		foreach $mol (@molecules) { if (($ipressure{$mol}) and ($fpressure{$mol})) { print OUT " end % P$mol\n"; }; };
		foreach $cat (@catalysts) { if (($icov{$cat}) and ($fcov{$cat})) { print OUT " end % coverage$cat\n"; }; };
  		if ($ttemp) { print OUT " end % Temperature\n"; };
        	if (($ipH) and ($fpH) and ($ipH != $fpH)) { print OUT " end % pH\n"; };
        	if (($nVext) and ($pVext) and ($nVext != $pVext)) { print OUT " end % Vext\n"; };
      }elsif ($exp eq "variable_TEMP") {
                print OUT "\n   for i = 1:$timeSteps\n      fprintf(fileID,\' %.2f  %.2f  %.2f  %.6f";
                foreach $Rs (@Rspecies) { print OUT "  %.15E"; };
                                          print OUT "\\n\',V,pH,T,results(i,:));\n   end\n";	      
		if ($ttemp) { print OUT " end % Temperature\n"; };
		foreach $mol (@molecules) { if (($ipressure{$mol}) and ($fpressure{$mol})) { print OUT " end % P$mol\n"; }; };
	        foreach $cat (@catalysts) { if (($icov{$cat}) and ($fcov{$cat})) { print OUT " end % coverage$cat\n"; }; };	
        	if (($ipH) and ($fpH) and ($ipH != $fpH)) { print OUT " end % pH\n"; };
        	if (($nVext) and ($pVext) and ($nVext != $pVext)) { print OUT " end % Vext\n"; };
      }elsif ($exp eq "RateControl") {
#                if ($ttemp) { print OUT " end % Temperature\n"; };
#                foreach $mol (@molecules) { if (($ipressure{$mol}) and ($fpressure{$mol})) { print OUT " end % P$mol\n"; }; };
#                foreach $cat (@catalysts) { if (($icov{$cat}) and ($fcov{$cat})) { print OUT " end % coverage$cat\n"; }; };
#                if (($ipH) and ($fpH) and ($ipH != $fpH)) { print OUT " end % pH\n"; };
#                if (($nVext) and ($pVext) and ($nVext != $pVext)) { print OUT " end % Vext\n"; };
#		print OUT "fclose(DRCfile,DSCfile);\n";
	}; # if exp
     close OUT;
    return();
  }; #--> sub ODE_closingloops
#==============================================================================================================================
   sub ODE_printing_sub {
        ($exp)=@_;
      open OUT, ">>$exp.m";
#    	 print OUT "\ndlmwrite(\'./KINETICS/DATA/$exp/solution$exp.dat\',sol,\'precision\',\'%1.15E\',\'delimiter\',\'\ \');\n";
        	print OUT "fclose(fileID);\n";			
          @gas=(); @su=(); @Rmolecules=();
	foreach $mol (@molecules) { foreach $rs (@react_species) {  $go="y"; foreach $t (@transitions) { if ($t eq $mol) { $go="n"; };};
                        if ($go eq "y") { if ($mol eq $rs) { push(@Rmolecules,$rs); };};};};
	foreach $rs (@react_species) { $do='yes'; foreach $mol (@molecules) { if ($mol eq $rs) { $do='no'; };};
                                             	  foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
# alberto 05/2019
                                             	  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
          if ($do eq 'yes') { push(@su,$rs); };};
        if ($exp eq "const_TEMP") { if ($nVext) { print OUT "fclose(fileID2);\n";}; };
	print OUT "\n\n";
      close OUT;
      open OUT, ">>plot.sh";
  	if ($exp eq "const_TEMP") { 
#  	  	print OUT "\n\n\nperl ~/software/KINETICS/IRPlotting.pl $exp";
		foreach $mol (@Rmolecules) { $n=4+@y{$mol}; print OUT " $n,$mol,mol"; };
	        foreach $s (@su) { $n=4+@y{$s};
	       	if ($s ne @su[-1]) { print OUT " $n,$s,surf"; }else{ print OUT " $n,$s,surf\n"; };};};
	if ($exp ne "TPR") {
#		print OUT "\n\nperl ~/software/KINETICS/PlotConcentrations.pl $exp";
                foreach $mol (@Rmolecules) { $n=4+@y{$mol}; print OUT " $n,$mol,mol"; };
                foreach $s (@su) { $n=4+@y{$s};
                if ($s ne @su[-1]) { print OUT " $n,$s,surf"; }else{ print OUT " $n,$s,surf\n"; };};};
      close OUT;
    return();
  }; #--> sub ODE_pringing
#==============================================================================================================================
   sub ODE_sub {
	($exp)=@_; if ($exp ne "TPR") { $fileout="myode"; }else{ $fileout="mytpr"; };   
       open OUT, ">>$fileout.m";
          print OUT "\nfunction dydt=$fileout(t,y,";
         foreach $trn (@transfer) { if ($trn ne @transfer[$#transfer]) { print OUT "$trn,"; }else{ print OUT "$trn)\n\n"; };};
         for ($pr=1; $pr<=$#process; $pr++ ) {
	    @PR=split(/\s+/,@ProcessReactants[$pr]); @PTS=split(/\s+/,@ProcessTS[$pr]); @PP=split(/\s+/,@ProcessProducts[$pr]);
	    ()=&stoichiometries_sub($fileout); };
       open OUT, ">>$fileout.m";
          print OUT "\n\n";
	  foreach $nr (@rates) { print OUT "$nr;\n"; };
	  print OUT "\n\n";
	    ()=&equations_sub($fileout);
          print OUT "end\n";
       close OUT;
    return();
   }; #--> sub ODE
#==============================================================================================================================
   sub stoichiometries_sub {
          ($fileout)=@_;
       open OUT, ">>$fileout.m";
         foreach $R (@PR) { @done{$R}="no"; }; foreach $P (@PP) { @done{$P}="no"; };
         foreach $R (@PR) { if (@done{$R} eq "no") { print OUT "stoichio$pr$R=$stoichio{$R}->[$pr];\t";
                                                     push(@trans,"stoichio$pr$R"); @done{$R}="yes";};};
         foreach $P (@PP) { if (@done{$P} eq "no") { printf OUT "stoichio$pr$P=$stoichio{$P}->[$pr];\t";
                                                     push(@trans,"stoichio$pr$P"); @done{$P}="yes";};};
            print OUT "\t%  process $pr\n";
       close OUT;
    return();
   }; #--> sub stoichiometries
#==============================================================================================================================
   sub equations_sub {
          ($fileout)=@_;
	  %Appequation=(); %Desequation=();
#-------------------------------------------------------------------------------------------------------------------------------- appear / desappear
   for ($pr=1; $pr<=$#process; $pr++ ) {
          @PR=split(/\s+/,@ProcessReactants[$pr]); @PTS=split(/\s+/,@ProcessTS[$pr]); @PP=split(/\s+/,@ProcessProducts[$pr]);
      if (($typeproc[$pr] eq 'A') or ($typeproc[$pr] eq 'a')) { $tmpmol=(); $tmpcat=();
        if ($exp ne "TPR") { 
            foreach $Rname (@PR) {
               foreach $mol (@molecules) { if ($Rname eq $mol) { $tmpmol=$mol;
		   if (!$Desequation{$mol}) { $Desequation{$mol}="rate$pr"; }else{ $Desequation{$mol}="$Desequation{$mol}+rate$pr"; };};};	      
               foreach $cat (@catalysts) { if ($Rname eq $cat) { $tmpcat=$cat;
		   if (!$Desequation{$cat}) { $Desequation{$cat}="rate$pr"; }else{ $Desequation{$cat}="$Desequation{$cat}+rate$pr"; };};}; };
            foreach $fpro (@PP) {
	       foreach $mol (@molecules) { if ($fpro eq $mol) { 
                   if (!$Appequation{$mol}) { $Appequation{$mol}="(stoichio$pr$fpro/stoichio$pr$tmpmol)*rate$pr";
		   }else{ $Appequation{$mol}="$Appequation{$mol}+(stoichio$pr$fpro/stoichio$pr$tmpmol)*rate$pr"; };};};
               foreach $cat (@catalysts) { if ($fpro eq $cat) { if ($tmpcat) { $tmp=$tmpcat; }else{ $tmp=$tmpmol; };
	           if (!$Appequation{$cat}) { $Appequation{$cat}="(stoichio$pr$fpro/stoichio$pr$tmpmol)*rate$pr";
                   }else{ $Appequation{$cat}="$Appequation{$cat}+(stoichio$pr$fpro/stoichio$pr$tmpmol)*rate$pr"; };};};};};
      }elsif (($typeproc[$pr] eq 'D') or ($typeproc[$pr] eq 'd')) { $tmpcat=();
            foreach $Rname (@PR) {
	       foreach $mol (@molecules) { if ($Rname eq $mol) {
	   	   if (!$Desequation{$mol}) { $Desequation{$mol}="rate$pr"; }else{ $Desequation{$mol}="$Desequation{$mol}+rate$pr"; };};}		       
               foreach $cat (@catalysts) { if ($Rname eq $cat) { $tmpcat=$cat;
		   if (!$Desequation{$cat}) { $Desequation{$cat}="rate$pr"; }else{ $Desequation{$cat}="$Desequation{$cat}+rate$pr"; };};}; };
            foreach $fpro (@PP) {
               foreach $mol (@molecules) { if ($fpro eq $mol) { $Pstoirate="$Pstoirate*stoichio$fpro";
                   if (!$Appequation{$mol}) { $Appequation{$mol}="(stoichio$pr$fpro/stoichio$pr$tmpcat)*rate$pr";
                   }else{ $Appequation{$mol}="$Appequation{$mol}+(stoichio$pr$fpro/stoichio$pr$tmpcat)*rate$pr"; };};};
               foreach $cat (@catalysts) { if ($fpro eq $cat) {
                   if (!$Appequation{$cat}) { $Appequation{$cat}="(stoichio$pr$fpro/stoichio$pr$tmpcat)*rate$pr";
                   }else{ $Appequation{$cat}="$Appequation{$cat}+(stoichio$pr$fpro/stoichio$pr$tmpcat)*rate$pr"; };};};};			   
      }else{ $tmpcat=();
            foreach $Rname (@PR) { $tmpcat=$Rname;
                   if (!$Desequation{$Rname}) { $Desequation{$Rname}="rate$pr"; }else{ $Desequation{$Rname}="$Desequation{$Rname}+rate$pr"; };};  
            foreach $fpro (@PP){
                   if (!$Appequation{$fpro}) { $Appequation{$fpro}="(stoichio$pr$fpro/stoichio$pr$tmpcat)*rate$pr";
                   }else{ $Appequation{$fpro}="$Appequation{$fpro}+(stoichio$pr$fpro/stoichio$pr$tmpcat)*rate$pr"; };};};  
   }; # for process
#--------------------------------------------------------------------------------------------------------------------------------
	  %tmp=(); $rows=0;
      foreach $rs (@react_species) { $go="no"; if (($Appequation{$rs}) or ($Desequation{$rs})) { $go="yes"; };
         if ($go eq "yes") { if (!$Appequation{$rs}) { $Appequation{$rs}=0; };
		             if (!$Desequation{$rs}) { $Desequation{$rs}=0; };
				     $equation{$rs}="$Appequation{$rs}-($Desequation{$rs})";
				     @EQ[@y{$rs}]=$equation{$rs}; $rows++; };};
	foreach $sur (@surfaces) { if (@y{$sur}) { $rows++; };};
    $rows=$rows-1;
    if ($fileout ne "RateControl") {
    	open OUT, ">>$fileout.m";
    };
              if ($fileout ne "const_TEMP") { print OUT " dydt=zeros($rows,1);\n"; };  @DRCeq=();
         foreach $mol (@molecules) { if ($equation{$mol}) { print OUT "   dydt(@y{$mol})=$equation{$mol}; % P$mol\n"; }; };
         foreach $rs (@react_species) { $go="yes"; $tr="no"; foreach $mol (@molecules) { if ($rs eq $mol) { $go="no"; $tr="yes"; };};
                                                             foreach $t (@transitions) { if ($rs eq $t) { $go="no"; $tr="yes"; };};
		                                             foreach $sur (@surfaces) { if ($rs eq $sur) { $go="no"; };};
			     if ($go eq "yes") { print OUT "   dydt(@y{$rs})=$equation{$rs}; % coverage$rs\n"; 
						 @tmp{@sitetype{$rs}}="@tmp{@sitetype{$rs}} dydt(@y{$rs})"; };};
         foreach $sur (@surfaces) { @balan=split(/\s+/,@tmp{@sitetype{$sur}});
             foreach $b (@balan) { if ($b) { push(@Nline,$b); }; }; @balan=@Nline; @Nline=();
	              print OUT "   dydt(@y{$sur})=-("; foreach $b (@balan) { if ($b ne @balan[$#balan]) { 
        	      print OUT "$b+"; }else{ print OUT "$b); % coverage$sur\n"; };};};
	print OUT "end % ends function\n"; 
     if ($fileout ne "RateControl") {
     close OUT;
     };
    return(\@EQ);
   }; #--> sub equations
#=============================================================================================================================
   sub Reset_Conditions_0_sub {
	($fileout)=@_;
	open OUT, ">>$fileout.m";
        $n=1;
        foreach $mol (@molecules) { foreach $rs (@react_species) {
                if ($mol eq $rs) {  $go="y"; foreach $t (@transitions) { if ($t eq $mol) { $go="n"; };};
                        if ($go eq "y") { print OUT "     rP$mol=concentrations(end,$n); "; $n++; };};};}; print OUT "\n";
        foreach $rs (@react_species) { $do='yes'; foreach $mol (@molecules) { if ($mol eq $rs) { $do='no'; };};
                                                  foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') { print OUT "     rcoverage$rs=concentrations(end,$n); "; $n++; };}; print OUT "\n";
        foreach $sur (@surfaces) { if ($y{$sur}) { print OUT "     rcoverage$sur=concentrations(end,$n); "; $n++; };}; print OUT "\n";
#	close OUT;
    return();
   }; #--> sub Reset_Conditions
#=============================================================================================================================
   sub Reset_Conditions_1_sub {
        ($fileout)=@_;
#        open OUT, ">>$fileout.m";
        $n=1;
        foreach $mol (@molecules) { foreach $rs (@react_species) {
                if ($mol eq $rs) {  $go="y"; foreach $t (@transitions) { if ($t eq $mol) { $go="n"; };};
                        if ($go eq "y") { print OUT "   P$mol=rP$mol; "; $n++; };};};}; print OUT "\n";
        foreach $rs (@react_species) { $do='yes'; foreach $mol (@molecules) { if ($mol eq $rs) { $do='no'; };};
                                                  foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') { print OUT "   coverage$rs=rcoverage$rs; "; $n++; };}; print OUT "\n";
        foreach $sur (@surfaces) { if ($y{$sur}) { print OUT "   coverage$sur=rcoverage$sur; "; $n++; };}; print OUT "\n";
#        close OUT;
    return();
   }; #--> sub Reset_Conditions
#=============================================================================================================================
   sub Reverse_Reactions_sub {
        ()=@_;
	@reverse=();
        for ($i=1; $i<=$#process; $i++ ) { $go="yes";
	       	$j=$i+1;
		while (($process[$j]) and (!$reverse[$i])) { @rev=();
			@PR1=split(/\s+/,@ProcessReactants[$i]);    @PP2=split(/\s+/,@ProcessProducts[$j]);
			for ($l=0; $l<=$#PR1; $l++) { $rev[$l]=0; foreach $P2 (@PP2) { if ($PR1[$l] eq $P2) { $rev[$l]=1; };};};
			foreach $r (@rev) { if ($r eq 0) { $go="no"; };};
			if ($go eq "yes") { @rev=();
				@PTS1=split(/\s+/,@ProcessTS[$i]);    @PTS2=split(/\s+/,@ProcessTS[$j]);
				for ($l=0; $l<=$#PTS1; $l++) { $rev[$l]=0; foreach $TS2 (@PTS2) { if ($PTS1[$l] eq $TS2) { $rev[$l]=1; };};};
				foreach $r (@rev) { if ($r eq 0) { $go="no"; };};
			};
			if ($go eq "yes") { @rev=();
                                @PP1=split(/\s+/,@ProcessProducts[$i]);     @PR2=split(/\s+/,@ProcessReactants[$j]);
                                for ($l=0; $l<=$#PP1; $l++) { $rev[$l]=0; foreach $R2 (@PR2) { if ($PP1[$l] eq $R2) { $rev[$l]=1; };};};
                                foreach $r (@rev) { if ($r eq 0) { $go="no"; };};
			};
                        if ($go eq "yes") {  $reverse[$i]=$j; };
		$j++; }; # while j
	}; # for process
    return(\@reverse);
   }; #--> sub Reset_Conditions
#=============================================================================================================================
   sub DRC_sub {
          ($fileout,@EQ)=@_;
	  $tincrease=$ftime*1E-3;

# 	Originally by Cambell :::		DOI: 10.1006/jcat.2001.3396
#	   					DOI: 10.1021/ja9000097
# 	              Kozuch and Shaik ::: 	DOI: 10.1021/ja0559146
#						DOI: 10.1021/jp8004772 	
#
#	"First-Principles-Based Microkinetics Simulations of Synthesis Gas Conversion on a Stepped Rhodium Surface"
#					ACS Catal. 2015, 5, 5453−5467  		DOI: 10.1021/acscatal.5b01391
#	>>> Degree of Rate Control
#	The DRC of a chemical reaction is defined as the relative change of the rate as a result of the 
#	relative change in the rate constant of a particular elementary reaction step while keeping the
#	equilibrium constant the same 
#
#					X_c,i = |(dR/R)/(dk_i/k_i)|K_i
# 					X_c,i = |dln(r_c)/dln(k_i)|K_i
#
#	R and r_c overall rate of species C	i elementary reaction step
#	k_i forward rate constant		K_i equilibrium reaction constant of step i
#
#					dR/R=(R'(t+dt)-R(t+dt))/R(t+dt)
#					dK_i/k_i= energy increment in the barrier ($BARRIER=1E-3)
#				
#		if X_c,i > 0 --> a decrease in Ea increases the overall rate
#		if X_c,i < 0 --> rate-inhibiting elementary step, if Ea decrease so does the overall rate
#		if X_c,i ~ 0 --> the elementary step does not affect at all the overall kinetics
#
     open OUT, ">>$fileout.m";
       	print OUT "\n";
        foreach $mol (@molecules) {  $go="y";
                foreach $t (@transitions) { if ($t eq $mol) { $go="n"; };};
                if ($go eq "y") {
                        if (((@pressure{$mol}) and (@pressure{$mol} != 0)) or ($ipressure{$mol})) { $mreactant=$mol; };};};
        foreach $rs (@react_species) { $do='yes';
	        foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
        	foreach $cat (@catalysts) { if ($cat eq $rs) { $do='no'; };};
	        foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
        	if (($do eq 'yes') and ($rs ne $mreactant)) {
                	if (($pressure{$rs}) or (@ipressure{$rs})) {
                        	print OUT "S$rs=(concentrations(end,$y{$rs})-concentrations(1,$y{$rs}))/(concentrations(1,$y{$mreactant})-concentrations(end,$y{$mreactant}));\n";
	};};}; # rs
        print OUT "\n";
	()=&Reset_Conditions_0_sub($fileout); open OUT, ">>$fileout.m";
	printf OUT "\n";
	for ($i=1; $i<=$#process; $i++ ) { print OUT "     K0rate$i=Krate$i;"; }; print OUT "\n";
	print OUT "\n";

	printf OUT "\n%%%%%%%%%%%%%%%%%%%%%%%%%%    R0    %%%%%%%%%%\n";
    	printf OUT "\n   tspan=[%.6f %.6f];\n   IC=[@Rspecies];\n",$ftime,$ftime+$tincrease;
        printf OUT "   options=odeset('Refine',5,'NonNegative',(1:%d),'RelTol',1e-7,'AbsTol',1e-5);\n",$#Rspecies+1;
        print OUT "solution=ode15s(\@(t,y) myode(t,y,";
 	foreach $trn (@transfer) { 
		if ($trn ne @transfer[$#transfer]) { print OUT "$trn,"; }else{ print OUT "$trn),tspan,IC,options);\n\n"; };
	};
        printf OUT "   concentrations=transpose(deval(solution,transpose(linspace(%.6f,%.6f,1))));\n",$ftime,$ftime+$tincrease;
        $n=1;
        foreach $mol (@molecules) { foreach $rs (@react_species) {
                if ($mol eq $rs) {  $go="y"; foreach $t (@transitions) { if ($t eq $mol) { $go="n"; };};
                        if ($go eq "y") { print OUT "   P$mol=concentrations(end,$n); "; $n++; };};};}; print OUT "\n";
        foreach $rs (@react_species) { $do='yes'; foreach $mol (@molecules) { if ($mol eq $rs) { $do='no'; };};
                                                  foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') { print OUT "   coverage$rs=concentrations(end,$n); "; $n++; };}; print OUT "\n";
        foreach $sur (@surfaces) { if ($y{$sur}) { print OUT "   coverage$sur=concentrations(end,$n); "; $n++; };}; print OUT "\n";
# R0
	foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') { print OUT "R0$rs=eval($EQ[$y{$rs}]);\n"; };};
# Ri Ki
	()=&Reverse_Reactions_sub();
	()=&DRC_call_sub($fileout); 
	()=&DSC_sub($fileout);
       	print OUT "\n";	
	()=&Reset_Conditions_1_sub($fileout); 
	print OUT "\n";
	()=&DRC_print_sub();
        ()=&DSC_print_sub();
	()=&ODE_closingloops_sub($fileout);
	close OUT;
    return();
   }; #--> sub DRC
#=============================================================================================================================
   sub DRC_call_sub {
        ($fileout)=@_;
	@done=();

	$BARRIER=1E-3;

	for ($pr=1; $pr<=$#process; $pr++ ) { if (!$done[$pr]) {
                printf OUT "\n%%%%%%%%%%%%%%%%%%%%%%%%%%    R$pr   @process[$pr]    %%%%%%%%%%\n";
                printf OUT "%%%%%%%%%%%%%%%%%%%%%%%%%%    R$reverse[$pr]   @process[$reverse[$pr]]    %%%%%%%%%%\n\n";

		for ($i=1; $i<=$#process; $i++ ) { if (($i eq $pr) or ($i eq $reverse[$pr])) {
			if (($typeproc[$i] eq "A") or ($typeproc[$i] eq "a")) {
        	       		print OUT "clear Krate$i;  Krate$i=process$i\{j,2}*process$i\{j,4}*exp(-((AEv$i+AEph$i+$BARRIER)*toeV)/(kb*T));\n";
	       		}else{  print OUT "clear Krate$i;  Krate$i=process$i\{j,3}*exp(-((AEv$i+AEph$i+$BARRIER)*toeV)/(kb*T));\n"; };
			$done[$pr]=1; $done[$reverse[$pr]]=1;
		}else{
			if (($typeproc[$i] eq "A") or ($typeproc[$reverse[$i]] eq "a")) {
        	               	print OUT "clear Krate$i;  Krate$i=process$i\{j,2}*process$i\{j,4}*exp(-((AEv$i+AEph$i)*toeV)/(kb*T));\n";
	                }else{  print OUT "clear Krate$i;  Krate$i=process$i\{j,3}*exp(-((AEv$i+AEph$i)*toeV)/(kb*T));\n"; };
		};};
                ()=&Reset_Conditions_1_sub($fileout); open OUT, ">>$fileout.m";
	        printf OUT "\n   tspan=[%.6f %.6f];\n   IC=[@Rspecies];\n",$ftime,$ftime+$tincrease;
        	printf OUT "   options=odeset('Refine',5,'NonNegative',(1:%d),'RelTol',1e-7,'AbsTol',1e-5);\n",$#Rspecies+1;
        	print OUT "solution=ode15s(\@(t,y) myode(t,y,";
        	foreach $trn (@transfer) {
                if ($trn ne @transfer[$#transfer]) { print OUT "$trn,"; }else{ print OUT "$trn),tspan,IC,options);\n\n"; };
        	};
	        printf OUT "   concentrations=transpose(deval(solution,transpose(linspace(%.6f,%.6f,1))));\n",$ftime,$ftime+$tincrease;
        	$n=1;
	        foreach $mol (@molecules) { foreach $rs (@react_species) {
        	        if ($mol eq $rs) {  $go="y"; foreach $t (@transitions) { if ($t eq $mol) { $go="n"; };};
                        if ($go eq "y") { print OUT "   P$mol=concentrations(end,$n); "; $n++; };};};}; print OUT "\n";
	        foreach $rs (@react_species) { $do='yes'; foreach $mol (@molecules) { if ($mol eq $rs) { $do='no'; };};
        	                                          foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                	                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
	                if ($do eq 'yes') { print OUT "   coverage$rs=concentrations(end,$n); "; $n++; };}; print OUT "\n";
        	foreach $sur (@surfaces) { if ($y{$sur}) { print OUT "   coverage$sur=concentrations(end,$n); "; $n++; };}; print OUT "\n";

                foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                          foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                        if ($do eq 'yes') { print OUT "R$pr$rs=eval($EQ[$y{$rs}]);\tR$reverse[$pr]$rs=eval($EQ[$y{$rs}]);\n"; };};
                print OUT "  Kirate$pr=Krate$pr;\tKirate$reverse[$pr]=Krate$reverse[$pr];\n";
		foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                          foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
			if ($do eq 'yes') { $X{$rs}->[$pr]="((R$pr$rs-R0$rs)/R0$rs)/((Kirate$pr-K0rate$pr)/K0rate$pr)";
                                            $X{$rs}->[$reverse[$pr]]="((R$reverse[$pr]$rs-R0$rs)/R0$rs)/((Kirate$reverse[$pr]-K0rate$reverse[$pr])/K0rate$reverse[$pr])"; };};
	}; # if done
	}; # for process
    return();
   }; #--> sub DRC_call
#=============================================================================================================================
   sub DRC_print_sub {
          ($fileout)=@_;
        print OUT "fprintf(DRCfile,\'# conditions :: V= %.2f pH= %.2f T= %.2f  "; printf OUT "t0= %.6f  t= %.6f",$ftime,$ftime+$tincrease;
        foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') {     if (($pressure{$rs}) or (@ipressure{$rs})) { print OUT "  $rs= %.3G"; };
                                        if (((@coverage{$rs}) and (@coverage{$rs} != 0)) or (@icoverage{$rs})) {
                                                print OUT "  $rs= %.3G"; };};};
         printf OUT "\\n\',V,pH,T";
        foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') {     if ($pressure{$rs}) { print OUT ",$pressure{$rs}";
                                        }elsif ($ipressure{$rs}) { print OUT ",0.5"; };
                                        if ((@coverage{$rs}) and (@coverage{$rs} != 0)) { print OUT ",$coverage{$rs}";
                                        }elsif (@icoverage{$rs}) { print OUT ",0.5"; };};};  print OUT ");\n";
        print OUT "fprintf(DRCfile,\'# %50s  %5s\\t";
        foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};

                foreach $mol (@molecules) { if ($do eq 'yes') { if ($mol eq $rs) { 
					if (($pressure{$rs}) or ($ipressure{$rs})) { print OUT "  %7s"; };};};};
#                                       foreach $cat (@catalysts) { if ($cat eq $rs) { print OUT "  %7s"; };};
                                };
         print OUT "\\n\',\'Reaction\',\'process\'";
        foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                                                  foreach $cat (@catalysts) { if ($cat eq $rs) { $do='no'; };};
                if ($do eq 'yes') { if (($pressure{$rs}) or ($ipressure{$rs})) { print OUT ",\'X$rs\'"; };};};
         print OUT ");\n";
        for ($pr=1; $pr<=$#process; $pr++ ) { $reaction=();
                @PR=split(/\s+/,@ProcessReactants[$pr]); @PTS=split(/\s+/,@ProcessTS[$pr]); @PP=split(/\s+/,@ProcessProducts[$pr]);
                foreach $r (@PR) { if ($r eq @PR[$#PR]) { $reaction.="$stoichio{$r}->[$pr].$r";
                        }else{ $reaction.="$stoichio{$r}->[$pr].$r + ";};};
                $reaction.=" --> ";
                foreach $p (@PP) { if ($p eq @PP[$#PP]) { $reaction.="$stoichio{$p}->[$pr].$p";
                        }else{ $reaction.="$stoichio{$p}->[$pr].$p + ";};};
                print OUT "fprintf(DRCfile,\'  %50s     %d\\t";
        foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                foreach $mol (@molecules) { if ($do eq 'yes') { if ($mol eq $rs) {
				       if (($pressure{$rs}) or ($ipressure{$rs})) { print OUT "  %+.4f"; };};};};
#                                        foreach $cat (@catalysts) { if ($cat eq $rs) { print OUT "  %+.4f"; };};
                                };
# Alberto 04/2019
                 print OUT "\\n\',\'$reaction\',$pr";
        foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                                                  foreach $cat (@catalysts) { if ($cat eq $rs) { $do='no'; };};
                if ($do eq 'yes') { if (($pressure{$rs}) or ($ipressure{$rs})) { print OUT ",$X{$rs}->[$pr]"; };};};
        print OUT ");\n";
        }; # for process
        print OUT "fclose(DRCfile);\n";
    return();
   }; #--> sub DRC_print

#=============================================================================================================================   
   sub DSC_sub {
          ($fileout)=@_;

#       Originally by Cambell :::               DOI: 10.1006/jcat.2001.3396
#
#       "First-Principles-Based Microkinetics Simulations of Synthesis Gas Conversion on a Stepped Rhodium Surface"
#                                       ACS Catal. 2015, 5, 5453−5467           DOI: 10.1021/acscatal.5b01391
#
#       >>> Degree of Selectivty Control
#       The DSC quantifies the extent to which a particular elementary reaction
#          step influences the selectivity to certain products
#       The sensitivity of the absolute change in selectivity as a result of the
#          relative change in the rate constant of a particular elementary reaction step.
#       while keeping the equilibrium constant the same 
#
#                                        S_c = [c]/([A0]-[A])
#                                       E_c,i = |dS_c/dln(k_i)|K_i
#                                       E_c,i = S_c (X_c,i - X_reactant,i)
#
#       S_c selectivity of species C over reactant A
#       k_i forward elemetary rate constant               K_i equilibrium reaction constant of step i
#     
#       if S_c1,i > 0 and S_c2,i < 0 --> inverse correltation -> when c1 increases c2 decreases
#       if S_c1,i > 0 and S_c2,i > 0 --> direct correltation -> when c1 increases c2 increases
#
#       correlation coefficient ρ_c1,c2,i based on the DSC parameters 
#
#                                       ρ_c1,c2,i = - E_c1,i * E_c2,i
#
#       if ρ_c1,c2,i < 0 --> i controls the formation of c1 and c2 in a competitive manner
#
#
        @done=(); 
	$mreactant=(); $creactant=();
        foreach $mol (@molecules) {  $go="y"; foreach $t (@transitions) { if ($t eq $mol) { $go="n"; };};
                        if ($go eq "y") { if (((@pressure{$mol}) and (@pressure{$mol} != 0)) or ($ipressure{$mol})) { $mreactant=$mol; };};};
        foreach $cat (@catalysts) { if (((@coverage{$cat}) and (@coverage{$cat} != 0)) or (@icoverage{$cat})) { $creactant=$cat; };};
        for ($pr=1; $pr<=$#process; $pr++ ) { if (!$done[$pr]) {
		foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                          foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                	foreach $mol (@molecules) { if ($do eq 'yes') { if ($mol ne $mreactant) {
				if (($pressure{$rs}) or ($ipressure{$rs})) {
# Alberto 11/2019
                        $DSC{$mol}->[$pr]="S$mol * ( $X{$mol}->[$pr] - $X{$mreactant}->[$pr] )";
			$DSC{$mol}->[$reverse[$pr]]="S$mol * ( $X{$mol}->[$reverse[$pr]] - $X{$mreactant}->[$reverse[$pr]] )";
	 };};};};};
#        foreach $cat (@catalysts) { if ($cat ne $creactant) {
#                $S{$cat}->[$pr]="concentrations(end,$y{$cat})/(concentrations(1,$y{$creactant})-concentrations(end,$y{$creactant}))";
#                $S{$cat}->[$reverse[$pr]]="concentrations(end,$y{$cat})/(concentrations(1,$y{$creactant})-concentrations(end,$y{$creactant}))";
#       		$DSC{$cat}->[$pr]="$S{$cat}->[$pr]/(Kirate$pr/K0rate$pr)";
#                $DSC{$cat}->[$reverse[$pr]]="$S{$cat}->[$reverse[$pr]]/(Kirate$reverse[$pr]/K0rate$reverse[$pr])"; };};
	$done[$pr]=1; $done[$reverse[$pr]]=1;
	};}; # for pr
    return();
   }; #--> sub DSC
#=============================================================================================================================   
   sub DSC_print_sub {
          ($fileout)=@_;
#    open OUT, ">>$fileout.m";
        print OUT "\n"; $mreactant=(); $creactant=();
        foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
						  foreach $cat (@catalysts) { if ($cat eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') { if (((@pressure{$rs}) and (@pressure{$rs} != 0)) or ($ipressure{$rs})) { $mreactant=$rs; };};};

	print OUT "fprintf(DSCfile,\'# conditions :: V= %.2f pH= %.2f T= %.2f  "; printf OUT "t0= %.6f  t= %.6f",$ftime,$ftime+$tincrease;
        foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') {     if (($pressure{$rs}) or (@ipressure{$rs})) { print OUT "  $rs= %.3G"; };
                                        if (((@coverage{$rs}) and (@coverage{$rs} != 0)) or (@icoverage{$rs})) {
                                                print OUT "  $rs= %.3G"; };};};
         printf OUT "\\n\',V,pH,T";
        foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') {     if ($pressure{$rs}) { print OUT ",$pressure{$rs}";
                                        }elsif ($ipressure{$rs}) { print OUT ",0.5"; };
                                        if ((@coverage{$rs}) and (@coverage{$rs} != 0)) { print OUT ",$coverage{$rs}";
                                        }elsif (@icoverage{$rs}) { print OUT ",0.5"; };};};  print OUT ");\n";

        print OUT "fprintf(DSCfile,\'# %50s  %5s\\t"; 
        foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') { if ($rs ne $mreactant) { if (($pressure{$rs}) or ($ipressure{$mol})) { 
			print OUT " %7s"; };};};};  # S
			print OUT "    ";
        foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') { if ($rs ne $mreactant) { if (($pressure{$rs}) or ($ipressure{$rs})) { 
			print OUT " %7s"; };};};};   # DSC
        		print OUT "    "; %done=();
        foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') { if ($rs ne $mreactant) { if (($pressure{$rs}) or ($ipressure{$rs})) { 
			foreach $mol2 (@molecules) { $go="yes"; foreach $t (@transitions) { if ($t eq $mol2) { $go='no'; };};
                               	if ($go eq 'yes') { if (($mol2 ne $mreactant) and ($rs ne $mol2) and (!$done{$rs}{$mol2})) {
                                       	if (($pressure{$mol2}) or ($ipressure{$mol2})) { $done{$mol2}{$rs}=1; 
			print OUT " %7s"; };};};};};};};};  # DSC_i * DSC_j
        		print OUT "\\n\',\'Reaction\',\'process\'";
    	foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') { if ($rs ne $mreactant) { if (($pressure{$rs}) or ($ipressure{$mol})) {  
                        print OUT ", \'S$rs\'"; };};};};  # S
        foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') { if ($rs ne $mreactant) { if (($pressure{$rs}) or ($ipressure{$rs})) { 
                        print OUT ", \'DSC$rs\'"; };};};};   # DSC   
       	%done=();
        foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                if ($do eq 'yes') { if ($rs ne $mreactant) { if (($pressure{$rs}) or ($ipressure{$rs})) {                       
                        foreach $mol2 (@molecules) { $go="yes"; foreach $t (@transitions) { if ($t eq $mol2) { $go='no'; };};
                                if ($go eq 'yes') { if (($mol2 ne $mreactant) and ($rs ne $mol2) and (!$done{$rs}{$mol2})) {
                                        if (($pressure{$mol2}) or ($ipressure{$mol2})) { $done{$mol2}{$rs}=1;
                        print OUT ", \'-(DSC$rs*DSC$mol2)\'"; };};};};};};};};  # DSC_i * DSC_j
         		print OUT ");\n";

        for ($pr=1; $pr<=$#process; $pr++ ) { $reaction=(); 
                @PR=split(/\s+/,@ProcessReactants[$pr]); @PTS=split(/\s+/,@ProcessTS[$pr]); @PP=split(/\s+/,@ProcessProducts[$pr]);
                foreach $r (@PR) { if ($r eq @PR[$#PR]) { $reaction.="$stoichio{$r}->[$pr].$r";
                        }else{ $reaction.="$stoichio{$r}->[$pr].$r + ";};};
                $reaction.=" --> ";
                foreach $p (@PP) { if ($p eq @PP[$#PP]) { $reaction.="$stoichio{$p}->[$pr].$p";
                        }else{ $reaction.="$stoichio{$p}->[$pr].$p + ";};};
# alberto 04/2019 + 11/2019
		print OUT "fprintf(DSCfile,\'  %50s   %d\\t";
    		foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                	                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
	                if ($do eq 'yes') { if ($rs ne $mreactant) { if (($pressure{$rs}) or ($ipressure{$mol})) {
        	                print OUT " %+.4f"; };};};};  # S
			        print OUT "    ";
        	foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                	                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
	                if ($do eq 'yes') { if ($rs ne $mreactant) { if (($pressure{$rs}) or ($ipressure{$rs})) {
        	                print OUT " %+.4f"; };};};};   # DSC   
			        print OUT "    "; %done=();
        	foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
	                                                  foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
        	        if ($do eq 'yes') { if ($rs ne $mreactant) { if (($pressure{$rs}) or ($ipressure{$rs})) {
                	        foreach $mol2 (@molecules) { $go="yes"; foreach $t (@transitions) { if ($t eq $mol2) { $go='no'; };};
                        	        if ($go eq 'yes') { if (($mol2 ne $mreactant) and ($rs ne $mol2) and (!$done{$rs}{$mol2})) {
                                	        if (($pressure{$mol2}) or ($ipressure{$mol2})) { $done{$mol2}{$rs}=1;
	                        print OUT " %+.4f"; };};};};};};};};  # DSC_i * DSC_j
		        	print OUT "\\n\',\'$reaction\',    $pr";
                foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                          foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                        if ($do eq 'yes') { if ($rs ne $mreactant) { if (($pressure{$rs}) or ($ipressure{$mol})) { 
                                print OUT ", S$rs"; };};};};  # S
                foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                          foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                        if ($do eq 'yes') { if ($rs ne $mreactant) { if (($pressure{$rs}) or ($ipressure{$rs})) {
                                print OUT ", $DSC{$rs}->[$pr]"; };};};};   # DSC   
                %done=();
                foreach $rs (@react_species) { $do='yes'; foreach $t (@transitions) { if ($t eq $rs) { $do='no'; };};
                                                          foreach $sur (@surfaces) { if (($sur eq $rs) and (!@freq{$sur})) {$do='no'; };};
                        if ($do eq 'yes') { if ($rs ne $mreactant) { if (($pressure{$rs}) or ($ipressure{$rs})) {
                                foreach $mol2 (@molecules) { $go="yes"; foreach $t (@transitions) { if ($t eq $mol2) { $go='no'; };};
                                        if ($go eq 'yes') { if (($mol2 ne $mreactant) and ($rs ne $mol2) and (!$done{$rs}{$mol2})) {
                                                if (($pressure{$mol2}) or ($ipressure{$mol2})) { $done{$mol2}{$rs}=1;
                                print OUT ", -( $DSC{$rs}->[$pr]*$DSC{$mol2}->[$pr])" };};};};};};};};  # DSC_i * DSC_j
         			print OUT ");\n";
 	}; # for process
	print OUT "fclose(DSCfile);\n";
#        if ($ttemp) { print OUT " end % Temperature\n"; };
#       foreach $mol (@molecules) { if (($ipressure{$mol}) and ($fpressure{$mol})) { print OUT " end % P$mol\n"; }; };
#       foreach $cat (@catalysts) { if (($icov{$cat}) and ($fcov{$cat})) { print OUT " end % coverage$cat\n"; }; };
#        if (($ipH) and ($fpH) and ($ipH != $fpH)) { print OUT " end % pH\n"; };
#        if (($nVext) and ($pVext) and ($nVext != $pVext)) { print OUT " end % Vext\n"; };
    return();
   }; #--> sub DSC
#==============================================================================================================================
#==============================================================================================================================



#==============================================================================================================================







   sub example_sub {
print "#    ---- Alberto Roldan ----\n#    v1 --> 06-2014\n#    coment line\n\n";
print "KINET  = TRUE\nRTEMP  = 50 1000 10\nRTIME  = 1 0.01\n##RVEXT  = -1 1 0.5\n#RpH    = 0 0 1\n\n";
print "#--------------------------------------------------------------------------------------- R1\n";
print "PROCESS = A N2H4 + Cu > N2H4_Cu\nSTOICHIOMETRY = 1 + 3 > 1\n#INEEXCH = 0\n";
print "PROCESS = D N2H4_Cu > N2H4 + Cu\nSTOICHIOMETRY = 1 > 1 + 3\n";
print "#--------------------------------------------------------------------------------------- R2\n";
print "PROCESS = R N2H4_Cu > TS_N2H4_N2H3_Cu >> N2H3_H_Cu\nSTOICHIOMETRY = 1 > 1 \nIREDOX = R 2\nIPH = H\n";
print "PROCESS = R  N2H3_H_Cu > TS_N2H4_N2H3_Cu >> N2H4_Cu\nSTOICHIOMETRY = 1 > 1\nIREDOX = O 2\nIPH = DH\n";
print "#===========================================================================================\n";
print "#---------------------------------------------------------------------- gas phase molecules\n";
print "SYSTEM = N2H4\nSYSPATH = '/home/alberto/ALBERTO/saeedeh-kinetics/OUTCARs/N2H4'\n";
print "MOLECULE = TRUE\nRPRESSURE = 50650 1013250 50650       # from 1/2*atm to 10*atm\n";
print "ISITES = 3 A\nEND\n";
print "#---------------------------------------------------------------------- Naked Cu(111)\n";
print "SYSTEM = Cu\nSYSPATH = '/home/alberto/ALBERTO/saeedeh-kinetics/OUTCARs/Cu'\n";
print "IACAT =   9.819E-20   # m^2\n#DEGENERATION = 1\nSSURFACE = A\nend\n";
print "#---------------------------------------------------------------------- On surface\n";
print "SYSTEM = N2H4_Cu\nSYSPATH = 'path../N2H4_Cu'\nIACAT =    9.819E-20  # m^2\n";
print "#DEGENERATION = 1\nICOVERAGE = 0.0\n##RCOVERAGE = 0.1 0.9 0.1\nTSITES = A\nend\n";
print "SYSTEM = N2H3_H_Cu\nSYSPATH = 'path../N2H3_H_Cu'\nIACAT =    9.819E-20  # m^2\n";
print "#DEGENERATION = 1\nICOVERAGE = 0.0\n##RCOVERAGE = 0.1 0.9 0.1\nTSITES = A\nend\n";
print "SYSTEM = TS_N2H4_N2H3_Cu\nSYSPATH = '../TS_N2H4_N2H3_H_Cu'\nIACAT =    9.819E-20  # m^2\n";
print "#DEGENERATION = 1\nICOVERAGE = 0.0\n##RCOVERAGE = 0.1 0.9 0.1\nTSITES = A\nend\n\n";
print "#---------------------------------------------- WHAT AND HOW TO DO --------------------------------------------------------------\n\n";
print "#============================================== GLOBAL PARAMETERS ===============================================================\n";
print "#THERMO = TRUE			# thermodynamic analysis\n";
print "#KINET = TRUE			# microkintic analysis\n";
print "#RTEMP = 10 500 50		# temperature ramp in Kelvins from 10 to 500 with steps of 5, for microkinetic and TPD analysis; avoid non natural values, Ti>=Tstep.\n";
print "#RTIME = 1 0.1			# time ramp in arbitrary units from 0 to 1 by steps of 0.1, for microkinetic and TPD analysis\n";
print "#end                            # each SYSTEM should have an end\n";
print "#___________________________________________________ MOLECULES __________________________________________________________________\n";
print "#IPRESSURE = 101325             # initial pressure for the molcule in Pa ; DEFAULT 101325 Pa = 1 atm\n";
print "#RPRESSURE = iP fP stepP        # ramp pressure                   \n";
print "#ISITES = 1    			# sites_per_molecule__1site, same value than stoichiometry of surf in Adsorption process and type of sites, one type per molecule-description\n";
print "#SYMMETRY = 2			# in case it is not proper in the OUCTAR\n";
print "#____________________________________________________ SURFACES ___________________________________________________________________\n";
print "#IACAT = 7e-20     		# Each catalyst site area in m^2\n";
print "#ICOVERAGE = 0.1		# coverage for this system (molecule on catalyst)\n";
print "#RCOVERAGE iCOV fCOV stepCOV    # ramp coverage\n";
print "#SSURFACE = A B			# Sites Naked Surface, when the surface's sites are the free sites\n";
print "#TSITES = A                     # Occupied sites\n#.\n#.\n\n";
print "#=============================================== PROCESS PARAMETERS ==============================================================\n";
print "#PROCESS = A a + b > c		# avoid numbers, the different processes must be ordered as the reaction profile\n";
print "#STOICHIOMETRY = 1 + 2 > 1	# same order and format than PROCESS, excluding PROCESS-TYPE( A, D, R)\n";
print "#INEEXCH = 1			# number of electrons in this process strongly related with the number of H.\n";
   }; #--> sub example
#==============================================================================================================================
