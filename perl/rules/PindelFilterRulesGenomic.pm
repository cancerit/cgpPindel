####################################################
# Copyright (c) 2013 Genome Research Ltd.
# Author: Cancer Genome Project, cgpit@sanger.ac.uk
# See LICENCE.TXT for details
####################################################
{
	tag  => 'INFO/LEN',
	name => 'F004',
	desc => 'Tum medium read depth strand bias check: Calls In 8% Reads Bt Depth 10 And 200 (inclusive)',
	test => sub {
#//		elsif($input->{'D MT'} >= 10 && $input->{'D MT'} < 200) {
#		if(pindelVariant.getRdMt() >= 10  && pindelVariant.getRdMt() < 200 ){
#//		    if($input->{'P MT +'} && $input->{'P MT -'}) {
#			if(pindelVariant.getpMtPos() > 0  && pindelVariant.getpMtNeg() > 0){
#//			    $pass = 1 if($input->{'P MT'} >= $input->{'D MT'} * 0.05);
#			//	System.out.println(pindelVariant.getpMt() + " " + pindelVariant.getRdMt()+ " " +(pindelVariant.getRdMt() * 0.05f));
#				if(pindelVariant.getpMt() >= (pindelVariant.getRdMt() * 0.05f)) ret = Boolean.TRUE;
#			}
#//		    elsif($input->{'P MT +'} && !$input->{'D MT +'}) {
#			else if(pindelVariant.getpMtPos() > 0 && pindelVariant.getRdMtPos() <= 0){
#				ret = Boolean.TRUE;
#			}
#//		    elsif($input->{'P MT +'}) {
#			else if(pindelVariant.getpMtPos() != null && pindelVariant.getpMtPos() > 0 ){
#//			    $pass = 1 if($input->{'P MT +'} >= $input->{'D MT +'} * 0.08);
#				if(pindelVariant.getpMtPos() >= (pindelVariant.getRdMtPos() * 0.08f)) ret = Boolean.TRUE;
#			}
#//		    elsif($input->{'P MT -'} && !$input->{'D MT -'}) {
#			else if(pindelVariant.getpMtNeg() > 0 && pindelVariant.getRdMtNeg() <= 0){
#				ret = Boolean.TRUE;
#			}
#//		    elsif($input->{'P MT -'}) {
#			else if(pindelVariant.getpMtNeg() != null && pindelVariant.getpMtNeg() > 0 ){
#//			    $pass = 1 if($input->{'P MT -'} >= $input->{'D MT -'} * 0.08);
#				if(pindelVariant.getpMtNeg() >= (pindelVariant.getRdMtNeg() * 0.08f)) ret = Boolean.TRUE;
#			}
#		}else{
#			ret = Boolean.TRUE;
#		}

		### HACK Dirty dirty dirty...... done to try and cut down the number of times I have to parse the FORMAT string I am storing it as a global variable.
		if($$RECORD[8] ne $main::previous_format_string){
			my @geno_formats = split(':',$$RECORD[8]);
			my $i = 0;
			map {$main::previous_format_hash->{$_} = $i++} split(':',$$RECORD[8]);
			$main::previous_format_string = $$RECORD[8];
		}

		my $ret = $FAIL;
		my @tum_geno = split(':',$$RECORD[10]);

		my $prd = $tum_geno[$main::previous_format_hash->{'PR'}];
		my $nrd = $tum_geno[$main::previous_format_hash->{'NR'}];
		my $rd = $prd + $nrd;
		my $pp = $tum_geno[$main::previous_format_hash->{'PU'}];
		my $np = $tum_geno[$main::previous_format_hash->{'NU'}];
		my $p = $pp + $np;

		if($rd >= 10 && $rd < 200){

# re written below to make easier to read....
#			if($pp > 0 && $np > 0){
#				$ret = $PASS if($p >= ($rd * 0.05));
#			}elsif($pp && !$prd){
#				$ret = $PASS;
#			}elsif($pp){
#				$ret = $PASS if($pp >= ($prd * 0.08));
#			}elsif($np && !$nrd){
#				$ret = $PASS;
#			}elsif($np){
#				$ret = $PASS if($np >= ($nrd * 0.08));
#			}


			if($pp > 0 && $np > 0){
				$ret = $PASS if($p >= ($rd * 0.05));
			}elsif($pp > 0){

				if($prd){
					$ret = $PASS if($pp >= ($prd * 0.08));
				}else{
					$ret = $PASS;
				}

			}elsif($np > 0){

				if($nrd){
					$ret = $PASS if($np >= ($nrd * 0.08));
				}else{
					$ret = $PASS;
				}

			}

		}else{
			$ret = $PASS;
		}
		return $ret;
	}
},

{
	tag  => 'INFO/LEN',
	name => 'F005',
	desc => 'Tum high read depth strand bias check: Calls In 4% Reads > Depth 200',
	test => sub {
#//		 elsif($input->{'D MT'} >= 200) {
#//			    if($input->{'P MT +'} && $input->{'P MT -'}) {
#//			      $pass = 1 if($input->{'P MT'} >= $input->{'D MT'} * 0.04);
#//			    }
#//			    elsif($input->{'P MT +'} && !$input->{'D MT +'}) {
#//			      $pass = 1;
#//			    }
#//			    elsif($input->{'P MT +'}) {
#//			      $pass = 1 if($input->{'P MT +'} >= $input->{'D MT +'} * 0.04);
#//			    }
#//			    elsif($input->{'P MT -'} && !$input->{'D MT -'}) {
#//			      $pass = 1;
#//			    }
#//			    elsif($input->{'P MT -'}) {
#//			      $pass = 1 if($input->{'P MT -'} >= $input->{'D MT -'} * 0.04);
#//			    }
#//			  }
#//			  else {
#//			    $pass = 1;
#//			  }
#		PindelVariant pindelVariant = getPindelVariant();
#		Boolean ret = Boolean.FALSE;
#
#		if(pindelVariant.getRdMt() >= 200 ){
#			if(pindelVariant.getpMtPos() > 0  && pindelVariant.getpMtNeg() > 0){
#				if(pindelVariant.getpMt() >= (pindelVariant.getRdMt() * 0.04f)) ret = Boolean.TRUE;
#			}
#			else if( pindelVariant.getpMtPos() > 0 && pindelVariant.getRdMtPos() <= 0){
#				ret = Boolean.TRUE;
#			}
#			else if(pindelVariant.getpMtPos() != null && pindelVariant.getpMtPos() > 0 ){
#				if(pindelVariant.getpMtPos() >= (pindelVariant.getRdMtPos() * 0.04f)) ret = Boolean.TRUE;
#			}
#			else if(&& pindelVariant.getpMtNeg() > 0 && pindelVariant.getRdMtNeg() <= 0){
#				ret = Boolean.TRUE;
#			}
#			else if(pindelVariant.getpMtNeg() != null && pindelVariant.getpMtNeg() > 0 ){
#				if(pindelVariant.getpMtNeg() >= (pindelVariant.getRdMtNeg() * 0.04f)) ret = Boolean.TRUE;
#			}
#		}else{
#			ret = Boolean.TRUE;
#		}

		### HACK Dirty dirty dirty...... done to try and cut down the number of times I have to parse the FORMAT string I am storing it as a global variable.
		if($$RECORD[8] ne $main::previous_format_string){
			my @geno_formats = split(':',$$RECORD[8]);
			my $i = 0;
			map {$main::previous_format_hash->{$_} = $i++} split(':',$$RECORD[8]);
			$main::previous_format_string = $$RECORD[8];
		}

		my $ret = $FAIL;
		my @tum_geno = split(':',$$RECORD[10]);

		my $prd = $tum_geno[$main::previous_format_hash->{'PR'}];
		my $nrd = $tum_geno[$main::previous_format_hash->{'NR'}];
		my $rd = $prd + $nrd;
		my $pp = $tum_geno[$main::previous_format_hash->{'PU'}];
		my $np = $tum_geno[$main::previous_format_hash->{'NU'}];
		my $p = $pp + $np;

		if($rd >= 200){

			if($pp > 0 && $np > 0){
				$ret = $PASS if($p >= ($rd * 0.04));
			}elsif($pp > 0){

				if($prd){
					$ret = $PASS if($pp >= ($prd * 0.04));
				}else{
					$ret = $PASS;
				}

			}elsif($np > 0){

				if($nrd){
					$ret = $PASS if($np >= ($nrd * 0.04));
				}else{
					$ret = $PASS;
				}
			}

		}else{
			$ret = $PASS;
		}
		return $ret;
	}

},

{
	tag  => 'INFO/LEN',
	name => 'F006',
	desc => 'Small call excessive repeat check: Fail if Length <= 4 and Repeats > 9',
	test => sub {
#if(pindelVariant.getLength() <= 4){
#			if(pindelVariant.getRepeats() > 9){
#				ret = Boolean.FALSE;
#			}
#		}

		if($MATCH <= 4){
			my ($rep) = $$RECORD[7] =~ /REP=(\d+)/;
			if($rep > 9) {
				return $FAIL;
			}
		}
		return $PASS;
	}
},

{
	tag  => 'INFO/LEN',                       # The VCF tag to apply this filter on
	name => 'F010',                       # The filter ID
	desc => "Variant must not exist within the Unmatched Normal Panel",  # Description for the VCF header
	test => sub {

#update PINDEL_VARIANT
#set FLAG_BIT_CODE = FLAG_BIT_CODE - ?
#where id_pindel_run = ?
#and id_variant in (
#	select distinct pv.id_variant
#	from pindel_variant pv, PINDEL_EXCL_MAT PE
#	where pe.data_type = ?
#	and pv.id_pindel_run = ?
#	and pv.VARIANT_TYPE = pe.type
#	and pv.id_sequence = pe.id_sequence
#	and pv.length <= 2
#	and pv.length between pe.length and pe.length_ext
#	and pv.MIN_POSITION between pe.lhs_pos_1 and pe.rhs_pos_1
#	union all
#	select distinct pv.id_variant
#	from pindel_variant pv, PINDEL_EXCL_MAT PE
#	where pe.data_type = ?
#	and pv.id_pindel_run = ?
#	and pv.VARIANT_TYPE = pe.type
#	and pv.id_sequence = pe.id_sequence
#	and pv.length > 2
#	and pv.length between pe.length and pe.length_ext
#	and pv.MIN_POSITION between pe.lhs_pos_20 and pe.rhs_pos_20)
#

		#22	16404839	.	GA	G	.	.	PC=D;RS=16404838;RE=16404857;LEN=1;SM=138;S1=10;S2=203.791;REP=18	PP:NP:PB:NB:PD:ND:PR:NR:PU:NU	1:0:1:0:4:5:4:6:1:1	3:1:4:2:25:20:32:21:18:8

		### HACK Dirty dirty dirty......
		unless($main::VCF_FLAGGING_UNMATCHED_NORMALS_TABIX){
			use Tabix;
			$main::VCF_FLAGGING_UNMATCHED_NORMALS_TABIX = new Tabix(-data => $ENV{VCF_FLAGGING_UNMATCHED_NORMALS},-index => $ENV{VCF_FLAGGING_UNMATCHED_NORMALS}.'.tbi');
		}

		my $length_off = ($MATCH <= 2) ? 1 : 20;

		my ($from) = ";$$RECORD[7]" =~ m/;RS=([^;]+)/;
		my ($to) = ";$$RECORD[7]" =~ m/;RE=([^;]+)/;
		#Range is the bases surrounding a position so need to bump back to the actual repetitive tract (or single pos)
		# but only when REP is > 0
		if(";$$RECORD[7];" !~ m/;REP=0;/) {
			$from++;
			$to--;
		}
		# then apply the range fudging
		$from -= ($length_off + 1); # additional -1 to switch to 0 based
		$to += $length_off;

		my $ret = eval{
			my $res = $main::VCF_FLAGGING_UNMATCHED_NORMALS_TABIX->query($CHROM,$from,$to);

			return $PASS if(!defined $res->get); # no valid entries (chromosome not in index) so must pass

			while(my $line = $main::VCF_FLAGGING_UNMATCHED_NORMALS_TABIX->read($res)){
				return $FAIL;
			}
			return $PASS;
		};
		if($@) {
	    die $@;
		}
		return $ret;
	},
},

{
	tag  => 'INFO/LEN',
	name => 'F012',
	desc => 'Germline: When length < 11 and depth > 9, fail if the variant is seen in both 20% of normal reads AND 20% of tumour reads in either pindel or bwa',
	test => sub {

#		if(pindelVariant.getLength() < 11 &&
#		pindelVariant.getdMt() > 9
#		&& pindelVariant.getdWt() > 9
#		&&(((pindelVariant.getpWt() >= (pindelVariant.getdWt().doubleValue() * 0.2d) && pindelVariant.getpMt() >= (pindelVariant.getdMt().doubleValue() * 0.2d))
#			||(pindelVariant.getbWt() >= (pindelVariant.getdWt().doubleValue() * 0.2d) && pindelVariant.getbMt() >= (pindelVariant.getdMt().doubleValue() * 0.2d)))
#		)

		if($MATCH < 11){

			### HACK Dirty dirty dirty...... to try and cut down the number of times I have to parse the FORMAT string I am storing it as a global variable.
			if($$RECORD[8] ne $main::previous_format_string){
				my @geno_formats = split(':',$$RECORD[8]);
				my $i = 0;
				map {$main::previous_format_hash->{$_} = $i++} split(':',$$RECORD[8]);
				$main::previous_format_string = $$RECORD[8];
			}

			my @nor_geno = split(':',$$RECORD[9]);
			my @tum_geno = split(':',$$RECORD[10]);

			if(($tum_geno[$main::previous_format_hash->{'ND'}] + $tum_geno[$main::previous_format_hash->{'PD'}]) > 9 &&
				 ($nor_geno[$main::previous_format_hash->{'ND'}] + $nor_geno[$main::previous_format_hash->{'PD'}]) > 9){

			 	my $tum_target_depth = ($tum_geno[$main::previous_format_hash->{'ND'}] + $tum_geno[$main::previous_format_hash->{'PD'}]) * 0.2;
			 	my $nor_target_depth = ($nor_geno[$main::previous_format_hash->{'ND'}] + $nor_geno[$main::previous_format_hash->{'PD'}]) * 0.2;

				if((($nor_geno[$main::previous_format_hash->{'NP'}] + $nor_geno[$main::previous_format_hash->{'PP'}] >= $nor_target_depth) &&
					($tum_geno[$main::previous_format_hash->{'NP'}] + $tum_geno[$main::previous_format_hash->{'PP'}] >= $tum_target_depth))
					||
					(($nor_geno[$main::previous_format_hash->{'NB'}] + $nor_geno[$main::previous_format_hash->{'PB'}] >= $nor_target_depth) &&
					($tum_geno[$main::previous_format_hash->{'NB'}] + $tum_geno[$main::previous_format_hash->{'PB'}] >= $tum_target_depth))
					){
					return $FAIL;
				}
			}

		}
		return $PASS;
	}
},

{
	tag  => 'INFO/LEN',
	name => 'F018',
	desc => 'Sufficient Depth: Pass if depth > 10',
	test => sub {
#		if(pindelVariant.getLength() < 50 && pindelVariant.getRdMt() >= 20 && pindelVariant.getRdWt() >= 20){
#			ret = Boolean.TRUE;
#		}

		### HACK Dirty dirty dirty...... done to try and cut down the number of times I have to parse the FORMAT string I am storing it as a global variable.
		if($$RECORD[8] ne $main::previous_format_string){
			my @geno_formats = split(':',$$RECORD[8]);
			my $i = 0;
			map {$main::previous_format_hash->{$_} = $i++} split(':',$$RECORD[8]);
			$main::previous_format_string = $$RECORD[8];
		}

		my @nor_geno = split(':',$$RECORD[9]);
		my @tum_geno = split(':',$$RECORD[10]);
		if(($nor_geno[$main::previous_format_hash->{'PR'}] + $nor_geno[$main::previous_format_hash->{'NR'}] >= 10) &&
			($tum_geno[$main::previous_format_hash->{'PR'}] + $tum_geno[$main::previous_format_hash->{'NR'}] >= 10)){
			return $PASS;
		}

		return $FAIL;
	}
},


{
	tag  => 'INFO/LEN',
	name => 'F015',
	desc => 'No normal calls',
	test => sub {
#if(pindelVariant.getpWt() == 0 && pindelVariant.getbWt() == 0){
#			ret = Boolean.TRUE;
#		}

		### HACK Dirty dirty dirty...... done to try and cut down the number of times I have to parse the FORMAT string I am storing it as a global variable.
		if($$RECORD[8] ne $main::previous_format_string){
			my @geno_formats = split(':',$$RECORD[8]);
			my $i = 0;
			map {$main::previous_format_hash->{$_} = $i++} split(':',$$RECORD[8]);
			$main::previous_format_string = $$RECORD[8];
		}

		my @nor_geno = split(':',$$RECORD[9]);

		if($nor_geno[$main::previous_format_hash->{'PU'}] + $nor_geno[$main::previous_format_hash->{'NU'}] == 0){
			return $PASS;
		}
		return $FAIL;
	}
},

{
	tag  => 'INFO/REP',
	name => 'F016',
	desc => 'Verify indel condition',
	test => sub {
#if(pindelVariant.getpMt() > 4){
#			if(pindelVariant.getbMt() > 0){
#				//small indel filter
#
#				ret = Boolean.TRUE;
#			}else if(pindelVariant.getbMt() == 0){
#				//large indel filter
#				if(pindelVariant.getbMt() == 0
#						&& pindelVariant.getRepeats() == 0
#						&& pindelVariant.getpMtNeg() > 0
#						&& pindelVariant.getpMtPos() > 0
#						){
#					ret = Boolean.TRUE;
#				}
#			}
#		}

		### HACK Dirty dirty dirty...... done to try and cut down the number of times I have to parse the FORMAT string I am storing it as a global variable.
		if($$RECORD[8] ne $main::previous_format_string){
			my @geno_formats = split(':',$$RECORD[8]);
			my $i = 0;
			map {$main::previous_format_hash->{$_} = $i++} split(':',$$RECORD[8]);
			$main::previous_format_string = $$RECORD[8];
		}

		my @nor_geno = split(':',$$RECORD[9]);
		my @tum_geno = split(':',$$RECORD[10]);
	
		if($tum_geno[$main::previous_format_hash->{'PP'}] + $tum_geno[$main::previous_format_hash->{'NP'}] > 4){
			if($tum_geno[$main::previous_format_hash->{'PB'}] + $tum_geno[$main::previous_format_hash->{'NB'}] > 0){
				return $PASS;
			}else{
				
				## RECORD[3]: ref
				## RECORD[4]: alt
				my $comp = length($$RECORD[3]) > 1 && length($$RECORD[4]) == 1 ? 1 : 0;
				
				if($MATCH == $comp){ ##repeats...
				
					if($tum_geno[$main::previous_format_hash->{'PP'}] > 0 &&
					$tum_geno[$main::previous_format_hash->{'NP'}] > 0){
						return $PASS;
					}
				}
			}
		}
		return $FAIL;
	}
},
