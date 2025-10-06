use strict;
use warnings;

our %configuration;

# Variable Type is two chars.
# The first char:
#  R for raw integrated variables
#  D for dgt integrated variables
#  S for raw step by step variables
#  M for digitized multi-hit variables
#  V for voltage(time) variables
#
# The second char:
# i for integers
# d for doubles

sub define_bank
{
	# uploading the hit definition
	my $bankId = 300;
	my $bankname = "solid_ec";
	
	insert_bank_variable(\%configuration, $bankname, "bankid", $bankId, "Di", "$bankname bank ID");
	insert_bank_variable(\%configuration, $bankname, "pid",     1,   "Di", "ID of the first particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "mpid",    2,   "Di", "ID of the mother of the first particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "tid",     3,   "Di", "Track ID of the first particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "mtid",    4,   "Di", "Track ID of the mother of the first particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "otid",    5,   "Di", "Track ID of the original track that generated the first particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "trackE",  6,   "Dd", "Energy of the track");
	insert_bank_variable(\%configuration, $bankname, "totEdep", 7,   "Dd", "Total Energy Deposited");
	insert_bank_variable(\%configuration, $bankname, "avg_x",   8,   "Dd", "Average X position in global reference system");
	insert_bank_variable(\%configuration, $bankname, "avg_y",   9,   "Dd", "Average Y position in global reference system");
	insert_bank_variable(\%configuration, $bankname, "avg_z",  10,   "Dd", "Average Z position in global reference system");
	insert_bank_variable(\%configuration, $bankname, "avg_lx", 11,   "Dd", "Average X position in local reference system");
	insert_bank_variable(\%configuration, $bankname, "avg_ly", 12,   "Dd", "Average Y position in local reference system");
	insert_bank_variable(\%configuration, $bankname, "avg_lz", 13,   "Dd", "Average Z position in local reference system");
	insert_bank_variable(\%configuration, $bankname, "px",     14,   "Dd", "x component of momentum of the particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "py",     15,   "Dd", "y component of momentum of the particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "pz",     16,   "Dd", "z component of momentum of the particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "vx",     17,   "Dd", "x component of primary vertex of the particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "vy",     18,   "Dd", "y component of primary vertex of the particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "vz",     19,   "Dd", "z component of primary vertex of the particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "mvx",    20,   "Dd", "x component of primary vertex of the mother of the particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "mvy",    21,   "Dd", "y component of primary vertex of the mother of the particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "mvz",    22,   "Dd", "z component of primary vertex of the mother of the particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "avg_t",  23,   "Dd", "Average time");
	insert_bank_variable(\%configuration, $bankname, "nsteps",  24,   "Di", "Number of geant4 steps");
	insert_bank_variable(\%configuration, $bankname, "procID",  25,   "Di", "Process that created the first particle");
	insert_bank_variable(\%configuration, $bankname, "totEdepB",  26,   "Dd", "");	
	insert_bank_variable(\%configuration, $bankname, "totEend",  27,   "Dd", "");	
	
	insert_bank_variable(\%configuration, $bankname, "Edep_seg0",  30,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Edep_seg1",  31,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Edep_seg2",  32,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Edep_seg3",  33,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Edep_seg4",  34,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Edep_seg5",  35,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Edep_seg6",  36,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Edep_seg7",  37,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Edep_seg8",  38,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Edep_seg9",  39,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "EdepB_seg0",  40,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "EdepB_seg1",  41,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "EdepB_seg2",  42,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "EdepB_seg3",  43,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "EdepB_seg4",  44,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "EdepB_seg5",  45,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "EdepB_seg6",  46,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "EdepB_seg7",  47,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "EdepB_seg8",  48,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "EdepB_seg9",  49,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Eend_seg0",  50,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Eend_seg1",  51,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Eend_seg2",  52,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Eend_seg3",  53,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Eend_seg4",  54,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Eend_seg5",  55,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Eend_seg6",  56,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Eend_seg7",  57,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Eend_seg8",  58,   "Dd", "");		
	insert_bank_variable(\%configuration, $bankname, "Eend_seg9",  59,   "Dd", "");		
	
	insert_bank_variable(\%configuration, $bankname, "id",     98,   "Di", "id number");		
	insert_bank_variable(\%configuration, $bankname, "hitn",   99,   "Di", "Hit Number");
}
define_bank();











