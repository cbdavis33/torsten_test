$PROB ABS ONE COMP PK MODEL

$PARAM @annotated
CL 		  : 5   : Clearance (volume/time)
VC 		  : 25  : Central Volume (volume)
KA 		  : 0.5 : Absorption Rate Constant (1/time)
F1 		  : 0.8 : Bioavailability (unitless)
LAGTIME :  0  : Lag Time (time)
  
$OMEGA @annotated @block
ETA_CL      : 0 :         ETA on CL
ETA_VC      : 0 0 :       ETA on VC
ETA_KA      : 0 0 0 :     ETA on KA
ETA_F1      : 0 0 0 0 :   ETA on F1
ETA_LAGTIME : 0 0 0 0 0 : ETA on LAGTIME

$SIGMA @annotated @block
SIGMA : 0 : PK lognormal Error


$CMT @annotated
GUT  : Extravascular (mass)
CENT : Central (mass)
  
$MAIN
double CLi = CL * exp(ETA_CL);
double VCi = VC * exp(ETA_VC);
double KAi = KA * exp(ETA_KA);
double F1i = F1 * exp(ETA_F1);
double LAGTIMEi = LAGTIME * exp(ETA_LAGTIME);

_F(1) = F1i;
_ALAG(1) = LAGTIMEi;

$ODE
dxdt_GUT = -KAi*GUT;
dxdt_CENT = KAi*GUT - (CLi/VCi)*CENT;

$TABLE
double IPRED = CENT/VCi;
double DV = IPRED*exp(SIGMA);

$CAPTURE @annotated
IPRED : Plasma Concentration (mass/volume)
DV    : Plasma Concentration with error (mass/volume)