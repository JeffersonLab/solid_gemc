use strict;
use warnings;

our %configuration;

sub solid_lmu_gem_hit
{
    # LMU-style urwell hit definition.
    # This matches geometry identifiers:
    #   region sector chamber layer component

    my %hit = init_hit();

    $hit{"name"}            = "lmu_gem";
    $hit{"description"}     = "LMU-style urwell/GEM strip hit definition";
    $hit{"identifiers"}     = "region sector chamber layer component";
    $hit{"signalThreshold"} = "0.1*KeV";
    $hit{"timeWindow"}      = "50*ns";
    $hit{"prodThreshold"}   = "1*mm";
    $hit{"maxStep"}         = "100*um";
    $hit{"delay"}           = "50*ns";
    $hit{"riseTime"}        = "1*ns";
    $hit{"fallTime"}        = "2*ns";
    $hit{"mvToMeV"}         = 100;
    $hit{"pedestal"}        = -20;

    print_hit(\%configuration, \%hit);
}

# Optional alias if your main driver expects define_hit()
sub define_hit
{
    solid_lmu_gem_hit();
}

1;
