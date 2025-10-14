use strict;
use warnings;

our %configuration;

# Variable Type is two chars.
# The first char:
#  R for raw integrated variables  //outdated
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
	my $bankId = 100;
	my $bankname = "solid_gem";
	
	insert_bank_variable(\%configuration, $bankname, "bankid", $bankId, "Di", "$bankname bank ID");
	insert_bank_variable(\%configuration, $bankname, "hitn",        99, "Di", "hit number");
	
	insert_bank_variable(\%configuration, $bankname, "ETot",         1, "Dd", "Total Energy Deposited");
	insert_bank_variable(\%configuration, $bankname, "x",            2, "Dd", "Average global x position");
	insert_bank_variable(\%configuration, $bankname, "y",            3, "Dd", "Average global y position");
	insert_bank_variable(\%configuration, $bankname, "z",            4, "Dd", "Average global z position");
	insert_bank_variable(\%configuration, $bankname, "lxin",         5, "Dd", "local x entrance position");
	insert_bank_variable(\%configuration, $bankname, "lyin",         6, "Dd", "local y entrance position");
	insert_bank_variable(\%configuration, $bankname, "lzin",         7, "Dd", "local z entrance position");
	insert_bank_variable(\%configuration, $bankname, "tin",          8, "Dd", "entrance time");
	insert_bank_variable(\%configuration, $bankname, "lxout",        9, "Dd", "local x exit position");
	insert_bank_variable(\%configuration, $bankname, "lyout",       10, "Dd", "local y exit position");
	insert_bank_variable(\%configuration, $bankname, "lzout",       11, "Dd", "local z exit position");
	insert_bank_variable(\%configuration, $bankname, "tout",        12, "Dd", "exit time");
	insert_bank_variable(\%configuration, $bankname, "pid",         13, "Di", "Particle ID");
	insert_bank_variable(\%configuration, $bankname, "vx",          14, "Dd", "y coordinate of primary vertex");
	insert_bank_variable(\%configuration, $bankname, "vy",          15, "Dd", "z coordinate of primary vertex");
	insert_bank_variable(\%configuration, $bankname, "vz",          16, "Dd", "x coordinate of primary vertex");
	insert_bank_variable(\%configuration, $bankname, "trE",         17, "Dd", "Energy of the track at the entrance point");
	insert_bank_variable(\%configuration, $bankname, "trid",        18, "Di", "Track ID");
	insert_bank_variable(\%configuration, $bankname, "weight",      19, "Dd", "Event weight(not used)");
	insert_bank_variable(\%configuration, $bankname, "px",          20, "Dd", "Momentum x");
	insert_bank_variable(\%configuration, $bankname, "py",          21, "Dd", "Momentum y");
	insert_bank_variable(\%configuration, $bankname, "pz",          22, "Dd", "Momentum z");
	insert_bank_variable(\%configuration, $bankname, "id",          23, "Di", "Volume ID");
	insert_bank_variable(\%configuration, $bankname, "mpid",        24, "Di", "ID of the mother of the first particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "mtid",        25, "Di", "Track ID of the mother of the first particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "otid",        26, "Di", "Track ID of the original track that generated the first particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "mvx",         27, "Dd", "x component of primary vertex of the mother of the particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "mvy",         28, "Dd", "y component of primary vertex of the mother of the particle entering the sensitive volume");
	insert_bank_variable(\%configuration, $bankname, "mvz",         29, "Dd", "z component of primary vertex of the mother of the particle entering the sensitive volume");	
	insert_bank_variable(\%configuration, $bankname, "nsteps",      30, "Di", "Number of geant4 steps");			
	insert_bank_variable(\%configuration, $bankname, "procID",      31, "Di", "Process that created the first particle");		
}
define_bank();

