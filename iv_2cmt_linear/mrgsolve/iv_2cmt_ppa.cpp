$PROB IV TWO COMP PK MODEL

$PARAM @annotated
TVCL 		  : 5   : Clearance (volume/time)
TVVC 		  : 25  : Central Volume (volume)
TVQ  		  : 5   : Intercompartmental Clearance (volume/time)
TVVP 		  : 25  : Peripiheral Volume (volume)

$PKMODEL cmt="CENT PERIPH", depot = FALSE  
  
$OMEGA @annotated @block
ETA_CL      : 0 :         ETA on CL
ETA_VC      : 0 0 :       ETA on VC
ETA_Q       : 0 0 0 :     ETA on Q
ETA_VP      : 0 0 0 0 :   ETA on VP

$SIGMA @annotated @block
SIGMA_ADD  : 0   : PK Additive Error
SIGMA_PROP : 0 0 : PK Proportional Error
  
$MAIN
double CL = TVCL * exp(ETA_CL);
double V1 = TVVC * exp(ETA_VC);
double Q  = TVQ  * exp(ETA_Q);
double V2 = TVVP * exp(ETA_VP);

$TABLE
double IPRED = CENT/V1;
double DV = IPRED * (1 + SIGMA_PROP) + SIGMA_ADD;

int i = 1;
while(DV < 0 & i <= 50) {
  simeps();
  DV = IPRED * (1 + SIGMA_PROP) + SIGMA_ADD;
  i++;
}

if(i > 50) mrg::report("Positive concentrations could not be obtained.");

$CAPTURE @annotated
IPRED : Plasma Concentration (mass/volume)
DV    : Plasma Concentration with error (mass/volume)