use strict;
use warnings;

our %configuration;
our %parameters;

sub define_hit
{
	# uploading the hit definition
	my %hit = init_hit();
	$hit{"name"}            = "solid_gem";
	$hit{"description"}     = "solid gem hit definition";
	$hit{"identifiers"}     = "id";
	$hit{"signalThreshold"} = "0*KeV";
	$hit{"timeWindow"}      = "0*ns";
	$hit{"prodThreshold"}   = "1*mm";
	$hit{"maxStep"}         = "1*mm";
	$hit{"delay"}           = "10*ns";
	$hit{"riseTime"}        = "1*ns";
	$hit{"fallTime"}        = "1*ns";
	$hit{"mvToMeV"}         = 100;
	$hit{"pedestal"}        = -20;
	print_hit(\%configuration, \%hit);
}
define_hit();
