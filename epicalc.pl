#!/usr/bin/perl
use warnings;
use strict;
 
#dmitrijs.lvovs@gmail.com, favorov@sensi.org
#calculates synergy factor as in Mario Cortina-Borja et al., “The Synergy Factor: a Statistic to Measure Interactions in Complex Diseases,” BMC Research Notes 2, no. 1 (2009): 105.

#also Fisher-like p-value is calculated, as in Douglas R White, Robert Pesner, and Karl P Reitz, “An Exact Significance Test for Three-Way Interaction Effects,” Cross-Cultural Research 18, no. 2 (May 1, 1983): 103–122.

#input as 8-pole contingency table:
#A	B	N	Y
#a	b	Nab	Yab
#A	b	NAb YAb
#a	B	NaB	YaB
#A	A	NAB	YAB


#program depends on normal approximaion by pjaclam:
#------------------------------------------------------------------------------------

sub ltqnorm ($) {
    #
    # Lower tail quantile for standard normal distribution function.
    #
    # This function returns an approximation of the inverse cumulative
    # standard normal distribution function.  I.e., given P, it returns
    # an approximation to the X satisfying P = Pr{Z <= X} where Z is a
    # random variable from the standard normal distribution.
    #
    # The algorithm uses a minimax approximation by rational functions
    # and the result has a relative error whose absolute value is less
    # than 1.15e-9.
    #
    # Author:      Peter John Acklam
    # Time-stamp:  2000-07-19 18:26:14
    # E-mail:      pjacklam@online.no
    # WWW URL:     http://home.online.no/~pjacklam
    
    my $p = shift;
    die "input argument must be in (0,1)\n" unless 0 < $p && $p < 1;
    
    # Coefficients in rational approximations.
    my @a = (-3.969683028665376e+01,  2.209460984245205e+02,
    -2.759285104469687e+02,  1.383577518672690e+02,
    -3.066479806614716e+01,  2.506628277459239e+00);
    my @b = (-5.447609879822406e+01,  1.615858368580409e+02,
    -1.556989798598866e+02,  6.680131188771972e+01,
    -1.328068155288572e+01 );
    my @c = (-7.784894002430293e-03, -3.223964580411365e-01,
    -2.400758277161838e+00, -2.549732539343734e+00,
    4.374664141464968e+00,  2.938163982698783e+00);
    my @d = ( 7.784695709041462e-03,  3.224671290700398e-01,
    2.445134137142996e+00,  3.754408661907416e+00);
    
    # Define break-points.
    my $plow  = 0.02425;
    my $phigh = 1 - $plow;
    
    # Rational approximation for lower region:
    if ( $p < $plow ) {
        my $q  = sqrt(-2*log($p));
        return ((((($c[0]*$q+$c[1])*$q+$c[2])*$q+$c[3])*$q+$c[4])*$q+$c[5]) /
        (((($d[0]*$q+$d[1])*$q+$d[2])*$q+$d[3])*$q+1);
    }
    
    # Rational approximation for upper region:
    if ( $phigh < $p ) {
        my $q  = sqrt(-2*log(1-$p));
        return -((((($c[0]*$q+$c[1])*$q+$c[2])*$q+$c[3])*$q+$c[4])*$q+$c[5]) /
        (((($d[0]*$q+$d[1])*$q+$d[2])*$q+$d[3])*$q+1);
    }
    
    # Rational approximation for central region:
    my $q = $p - 0.5;
    my $r = $q*$q;
    return ((((($a[0]*$r+$a[1])*$r+$a[2])*$r+$a[3])*$r+$a[4])*$r+$a[5])*$q /
    ((((($b[0]*$r+$b[1])*$r+$b[2])*$r+$b[3])*$r+$b[4])*$r+1);
}

#program depends on Factorials by favorov@sensi.org:
#------------------------------------------------------------------------------------

use constant PI => 3.141592653589793;

#counted as 0.5*log(PI), surprisingly :)
use constant log_sq_root_pi => 0.572364942924700;

#counted as 0.5*log(2*PI), surprisingly :)
use constant log_sq_root_two_pi => 0.918938533204673;

#from 0 to 40 including both
use constant first_40_factorials => (
1,                    1,
2,                    6,
24,                   120,
720,                  5040,
40320,                362880,
3628800,              39916800,
479001600,            6227020800,
87178291200,          1307674368000,
20922789888000,       355687428096000,
6402373705728000,     121645100408832000,
2432902008176640000,  5.10909421717094e+19,
1.12400072777761e+21, 2.5852016738885e+22,
6.20448401733239e+23, 1.5511210043331e+25,
4.03291461126606e+26, 1.08888694504184e+28,
3.04888344611714e+29, 8.8417619937397e+30,
2.65252859812191e+32, 8.22283865417792e+33,
2.63130836933694e+35, 8.68331761881189e+36,
2.95232799039604e+38, 1.03331479663861e+40,
3.71993326789901e+41, 1.37637530912263e+43,
5.23022617466601e+44, 2.03978820811974e+46,
8.15915283247898e+47
);

#from 0 to 40 including both
use constant first_40_factorial_logs => (
0,                0,                0.693147180559945, 1.79175946922805,
3.17805383034795, 4.78749174278205, 6.5792512120101,   8.52516136106541,
10.6046029027453, 12.8018274800815, 15.1044125730755,  17.5023078458739,
19.9872144956619, 22.5521638531234, 25.1912211827387,  27.8992713838409,
30.6718601060807, 33.5050734501369, 36.3954452080331,  39.3398841871995,
42.3356164607535, 45.3801388984769, 48.4711813518352,  51.6066755677644,
54.7847293981123, 58.0036052229805, 61.261701761002,   64.5575386270063,
67.8897431371815, 71.257038967168,  74.6582363488302,  78.0922235533153,
81.557959456115,  85.0544670175815, 88.5808275421977,  92.1361756036871,
95.7196945421432, 99.3306124547874, 102.968198614514,  106.631760260643,
110.320639714757
);

sub log_stirling_approx

# in scalar context returns the logarithm of Gosper's enrichment
# of Wells' variation of Striling's approximation for n!,
#
# in list context return 3 values : [0] is the same vith scalar,
# [1] and [2] are lower and upper boundaries of the Robbins-Feller
#	evaluation

# see # http://mathworld.wolfram.com/StirlingsApproximation.html
{
	my ($n)     = @_;
	my ($wells) =
    log_sq_root_pi + 0.5 * log( 2 * $n + 0.3333333333333 ) + $n * log($n) -
    $n;
	if (wantarray) {
        
		#List context;
		my ($robbins) =
        log_sq_root_two_pi + ( $n + 0.5 ) * log($n) - $n + 1 /
        ( 12 * $n + 1 );
		my ($feller) =
        log_sq_root_two_pi + ( $n + 0.5 ) * log($n) - $n + 1 / ( 12 * $n );
		( $wells, $robbins, $feller );
	}
    
	# False, but defined
	elsif ( defined wantarray ) {
        
		#Scalar context;
		$wells;
	}
    
	# False and undefined
	else {
        
		#void
		0;
	}
}

sub log_fact

#the first 0..40 are tabulated, others are from Striling's formulas
{
	my ($n) = @_;
	if ( $n <= 40 ) {
		my ($result) = (first_40_factorial_logs)[$n];
		$result;
	}
	else    #>40
	{
		my ($res) = log_stirling_approx($n);
		$res;
	}
}

sub fact

#the first 0..40 are tabulated, others are from Striling's formulas
{
	my ($n) = @_;
	if ( $n <= 40 ) {
		my ($result) = (first_40_factorials)[$n];
		$result;
	}
	else    #>40
	{
		my ($res) = exp( log_stirling_approx($n) );
		$res;
	}
}

#just a check that all input fields are numbers
#------------------------------------------------------------------------------------
sub notdigits(@){
#check that all supplied parameters are numbers
	my @m=@_;
	my @n=grep{/\d+/} @m;
	return $#m-$#n;
}

#this sub calculates FLINT weight as in White(1983)
#------------------------------------------------------------------------------------
sub flw_calc (@) {
    my @orcountarr=@_;
    
    #just a copy, in order not to calculate SF all the time
    my $Nab=$orcountarr[0];
    my $NAb=$orcountarr[1];
    my $NaB=$orcountarr[2];
    my $NAB=$orcountarr[3];
    my $Yab=$orcountarr[4];
    my $YAb=$orcountarr[5];
    my $YaB=$orcountarr[6];
    my $YAB=$orcountarr[7];
    
    my $YA=$YAb+$YAB;
    my $NA=$NAb+$NAB;
    my $Ya=$Yab+$YaB;
    my $Na=$Nab+$NaB;
    my $YB=$YaB+$YAB;
    my $NB=$NaB+$NAB;
    my $Yb=$Yab+$YAb;
    my $Nb=$Nab+$NAb;
    
    my $ab=$Yab+$Nab;
    my $Ab=$YAb+$NAb;
    my $aB=$YaB+$NaB;
    my $AB=$YAB+$NAB;
    
#we calculate and return log factorials, because otherwise we run into inf and resulting NaN
    
    my $flw=log_fact($YA) - (log_fact($YAB) + log_fact($YA - $YAB)) +
            log_fact($NA) - (log_fact($NAB) + log_fact($NA - $NAB)) +
            log_fact($Ya) - (log_fact($YaB) + log_fact($Ya - $YaB)) +
            log_fact($Na) - (log_fact($NaB) + log_fact($Na - $NaB));
    
    return $flw;
}

#this sub calculates the weights in a familily of tables
#------------------------------------------------------------------------------------
sub family_flw_calc(@){
	
#calculate the first table
    my @seed=@_;
    
    my $Yab=$seed[4];
    my $NAb=$seed[1];
    my $NaB=$seed[2];
    my $YAB=$seed[7];
    
#we need two-sided p-value, therefore going through all tables from max table
    my $upto=max(@seed);
    my $maxdelta=min($YAB,$Yab,$NaB,$NAb);
    
    
    $seed[0]=$seed[0]+$maxdelta;
    $seed[1]=$seed[1]-$maxdelta;
    $seed[2]=$seed[2]-$maxdelta;
    $seed[3]=$seed[3]+$maxdelta;
    $seed[4]=$seed[4]-$maxdelta;
    $seed[5]=$seed[5]+$maxdelta;
    $seed[6]=$seed[6]+$maxdelta;
    $seed[7]=$seed[7]-$maxdelta;

#calculate weights in family by walking downward to zero from max table
    
    my @allflw;
    my @allt;
    
    WALK: for (my $i=1; $i<=$upto; $i++) {
    
        if ($seed[0] lt 0 or $seed[1] lt 0 or $seed[2] lt 0 or $seed[3] lt 0 or $seed[4] lt 0 or $seed[5] lt 0 or $seed[6] lt 0 or $seed[7] lt 0) {
            last WALK}
        else {
            my $a=flw_calc(@seed);
            push(@allflw,$a);
        }
        
        $seed[0]=$seed[0]-1;
        $seed[1]=$seed[1]+1;
        $seed[2]=$seed[2]+1;
        $seed[3]=$seed[3]-1;
        $seed[4]=$seed[4]+1;
        $seed[5]=$seed[5]-1;
        $seed[6]=$seed[6]-1;
        $seed[7]=$seed[7]+1;
    }
    return @allflw;
}


#main
#------------------------------------------------------------------------------------
use List::Util 'max','min';
my $usage="Usage: perl epicalc.pl [Nab NAb NaB NAB Yab YAb YaB YAB] (--alpha n)\n";

#set defaults and check input
my @orcountarr;
my $alpha=0.05;

for (my $i=0; $i<=$#ARGV; $i++) {
 	unless ($ARGV[$i]=~/\D/) {push(@orcountarr,$ARGV[$i])}
 	if (($ARGV[$i] eq '--alpha') and (defined $ARGV[$i+1])) {$alpha=$ARGV[$i+1]}
 }
 
die "\"@orcountarr\" is not valid 8-pole contingency table\n", $usage if (notdigits(@orcountarr) ne 0) | $#orcountarr ne 7;

#calculate weights in families
my @allflw=family_flw_calc(@orcountarr);

#calculate my weight in families
my $myflw=flw_calc(@orcountarr);


#normalize by weight of our table in order to exponentiate smaller values
foreach my $w (@allflw) {
	$w=exp($w-$myflw);
}

#calculate p-value
my $allflwsum=0;
my $lteqflwsum=0;

for (my $i=0; $i<=$#allflw; $i++){

    $allflwsum=$allflwsum+$allflw[$i];
    
    if ($allflw[$i] <= 1) { #1 because $myflw-$myflw=0 and exp on that is 1
        $lteqflwsum=$lteqflwsum+$allflw[$i];
    }
}

#print STDERR "family is:", join("\n",(@allflw)),"\n";
#print STDERR "my table is:", join(" ",@orcountarr),"\n";
#print STDERR "my weight is:", $myflw,"\n";


#print STDERR "p-value is\n",$lteqflwsum, "\\", $allflwsum, "=",$lteqflwsum/$allflwsum,"\n";


#calculate SF as in (Cortina-Borja 2009)

#adjust if there are zeroes
my $string = join(" ",@orcountarr);
if ($string=~m/\b0\b/) {
	for (my $i=0; $i<=$#orcountarr; $i++){
		$orcountarr[$i]=$orcountarr[$i]+0.5;
	}
}

#contingency table..
my $Nab=$orcountarr[0];
my $Yab=$orcountarr[1];
my $NAb=$orcountarr[2];
my $YAb=$orcountarr[3];
my $NaB=$orcountarr[4];
my $YaB=$orcountarr[5];
my $NAB=$orcountarr[6];
my $YAB=$orcountarr[7];

my $YA=$YAb+$YAB;
my $NA=$NAb+$NAB;
my $Ya=$Yab+$YaB;
my $Na=$Nab+$NaB;
my $YB=$YaB+$YAB;
my $NB=$NaB+$NAB;
my $Yb=$Yab+$YAb;
my $Nb=$Nab+$NAb;

my $ab=$Yab+$Nab;
my $Ab=$YAb+$NAb;
my $aB=$YaB+$NaB;
my $AB=$YAB+$NAB;

my $ORab = ($Nab*$YAB)/($NAB*$Yab);
my $ORa = ($Nab*$YaB)/($NaB*$Yab);
my $ORb = ($Nab*$YAb)/($NAb*$Yab);

my $SF = $ORab/($ORa*$ORb);

my $SFse = sqrt(1/$YAB+1/$NAB+1/$YaB+1/$NaB+1/$YAb+1/$NAb+1/$Yab+1/$Nab);
my $SFlower = exp(log($SF)+$SFse*ltqnorm($alpha/2));
my $SFupper = exp(log($SF)+$SFse*ltqnorm(1-$alpha/2));


#printout results
print "SF=";
printf("%.2e",$SF);
print " CI(",(1-$alpha)*100,"%)","=[";
printf ("%.2e",$SFlower);
print "..";
printf("%.2e",$SFupper,);
print "]"," p-value=";
printf("%.2e",$lteqflwsum/$allflwsum);
print "\n";


