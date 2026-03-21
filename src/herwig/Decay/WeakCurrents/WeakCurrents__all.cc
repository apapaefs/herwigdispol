#line 1 "./FourPionNovosibirskCurrent.cc"
// -*- C++ -*-
//
// FourPionNovosibirskCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FourPionNovosibirskCurrent class.
//

#include "FourPionNovosibirskCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include <functional>

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

DescribeClass<FourPionNovosibirskCurrent,WeakCurrent>
describeHerwigFourPionNovosibirskCurrent("Herwig::FourPionNovosibirskCurrent",
					 "HwWeakCurrents.so");

HERWIG_INTERPOLATOR_CLASSDESC(FourPionNovosibirskCurrent1,Energy,Energy2)

HERWIG_INTERPOLATOR_CLASSDESC(FourPionNovosibirskCurrent2,double,Energy)

HERWIG_INTERPOLATOR_CLASSDESC(FourPionNovosibirskCurrent3,double,Energy2)

IBPtr FourPionNovosibirskCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr FourPionNovosibirskCurrent::fullclone() const {
  return new_ptr(*this);
}

void FourPionNovosibirskCurrent::doupdate() {
  WeakCurrent::doupdate();
  // update running width if needed
  if ( !touched() ) return;
  if(_maxmass!=_maxcalc) inita1width(-1);
}

FourPionNovosibirskCurrent::FourPionNovosibirskCurrent() : _mpic(), _mpi0(),
							   _mpic2(), _mpi02(), _prho()
{
  // set the number of modes
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  setInitialModes(2);
  // masses of the particles used in the current
  _rhomass    = 0.7761*GeV;
  _a1mass     = 1.2300*GeV;
  _omegamass  = 0.7820*GeV;
  _sigmamass  = 0.8000*GeV;
  // widths of the particles used in the current
  _rhowidth   = 0.14450*GeV;
  _a1width    = 0.45000*GeV;
  _omegawidth = 0.00841*GeV;
  _sigmawidth = 0.80000*GeV;
  // parameters for the resonance used in the integration
  _intmass = 1.4*GeV;
  _intwidth = 0.5*GeV;
  // relative coupling of the sigma and rho in the a_1 decay
  _zmag   = 1.3998721;
  _zphase = 0.43585036;
  _zsigma=0.;
  // parameter for f_a1
  _lambda2 = 1.2*GeV2;
  _onedlam2 = 1./_lambda2;
  _a1massolam2 = _a1mass*_a1mass*_onedlam2;
  _hm2=ZERO; 
  _rhoD=ZERO;
  _dhdq2m2=0.;
  // use local values of the parameters
  _localparameters=true;
  // magic numbers from TAUOLA (N.B. conversion from GeV to MeV)
  _athreec = 76.565/GeV;
  _bthreec = 0.71709;
  _cthreec = 0.27505;
  _aomega  = 886.84/GeV;
  _bomega  = 0.70983;
  _comega  = 0.26689;
  _aonec   = 96.867/GeV;
  _bonec   = 0.70907;
  _conec   = 0.26413;
  //parameters for the running omega width
  _omegaparam.resize(10);
  _omegaparam[0] = 17.560;
  _omegaparam[1] = 141.110;
  _omegaparam[2] = 894.884;
  _omegaparam[3] = 4977.35;
  _omegaparam[4] = 7610.66;
  _omegaparam[5] =-42524.4;
  _omegaparam[6] =-1333.26;
  _omegaparam[7] = 4860.19;
  _omegaparam[8] =-6000.81;
  _omegaparam[9] = 2504.97;
  // values of the g omega function from hep-ph/0201149
  double Fomegainit[98]={ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 2.2867811,
			  2.9710648, 2.9344304, 2.6913538, 2.5471206, 2.3557470,
			  2.2448280, 2.1074708, 2.0504866, 1.9270257, 1.8669430,
			  1.7907301, 1.7184515, 1.6535717, 1.6039416, 1.5535343,
			  1.5065620, 1.4608675, 1.4215596, 1.3849826, 1.3480113,
			  1.3147917, 1.2793381, 1.2487282, 1.2184237, 1.1952927,
			  1.1683835, 1.1458827, 1.1145806, 1.0935820, 1.0608720,
			  1.0390474, 1.0164336, 0.9908721, 0.9585276, 0.9307971,
			  0.9017274, 0.8731154, 0.8452763, 0.8145532, 0.7817339,
			  0.7493086, 0.7199919, 0.6887290, 0.6568120, 0.6255773,
			  0.5944664, 0.5661956, 0.5391204, 0.5102391, 0.4786543,
			  0.4546428, 0.4316779, 0.4063754, 0.3769831, 0.3561141,
			  0.3333555, 0.3139160, 0.2949214, 0.2814728, 0.2602444,
			  0.2349602, 0.2269845, 0.2192318, 0.2286938, 0.2839763,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000};
  // values of the three charged pion G function from hep-ph/0201149
  double Fthreeinit[98]={ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  13.1664906,10.7234087, 8.8219614,10.7989664, 9.1883001,
			  7.8526378, 7.7481031, 8.2633696, 5.5042820, 4.9029269,
			  4.4794345, 3.9654009, 4.5254011, 3.6509495, 3.5005512,
			  3.2274280, 3.1808922, 2.9925177, 2.6886659, 2.5195024,
			  2.4678771, 2.3540580, 2.2123868, 2.1103525, 2.0106986,
			  1.8792295, 1.8250662, 1.7068460, 1.6442842, 1.5503920,
			  1.4814349, 1.4225838, 1.3627135, 1.3205355, 1.2784383,
			  1.2387408, 1.1975995, 1.1633024, 1.1318133, 1.1114354,
			  1.0951439, 1.0691465, 1.0602311, 1.0392803, 1.0220672,
			  1.0154786, 1.0010130, 0.9908018, 0.9710845, 0.9602382,
			  0.9488459, 0.9316537, 0.9118049, 0.8920435, 0.8719332,
			  0.8520256, 0.8280582, 0.8064085, 0.7767881, 0.7570597,
			  0.7382626, 0.7100251, 0.6846500, 0.6666913, 0.6372250,
			  0.6162248, 0.6007728, 0.5799103, 0.5674670, 0.5446148,
			  0.5352115, 0.5128809, 0.4932536, 0.5310397, 0.8566489,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000};
  double   Foneinit[98]={ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  1.4819183, 1.7086354, 1.6958492, 1.6172935, 1.6301320,
			  1.5719706, 1.5459771, 1.5377471, 1.5008864, 1.4924121,
			  1.4720882, 1.4371741, 1.3990080, 1.3879193, 1.4030601,
			  1.3768673, 1.3493533, 1.3547127, 1.3275831, 1.3167892,
			  1.3035913, 1.2968298, 1.2801558, 1.2650299, 1.2557997,
			  1.2325822, 1.2210644, 1.1935984, 1.1746194, 1.1510350,
			  1.1358515, 1.1205584, 1.1010553, 1.0903869, 1.0731295,
			  1.0578678, 1.0438409, 1.0377911, 1.0253277, 1.0103551,
			  1.0042409, 0.9937978, 0.9858117, 0.9770346, 0.9724492,
			  0.9656686, 0.9606671, 0.9525813, 0.9488522, 0.9417335,
			  0.9399430, 0.9323438, 0.9281269, 0.9244171, 0.9237418,
			  0.9174354, 0.9177181, 0.9120840, 0.9047825, 0.9065579,
			  0.9034142, 0.8992961, 0.9011586, 0.9036470, 0.8954964,
			  0.8898208, 0.8911991, 0.8854824, 0.8888282, 0.8868449,
			  0.9004632, 0.8981572, 0.9096183, 0.9046990, 1.7454215,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000};
  // eninit in GeV
  double     eninit[98]={ 0.6000000, 0.6131313, 0.6262626, 0.6393939, 0.6525252,
			  0.6656566, 0.6787879, 0.6919192, 0.7050505, 0.7181818,
			  0.7313131, 0.7444444, 0.7575758, 0.7707071, 0.7838384,
			  0.7969697, 0.8101010, 0.8232324, 0.8363636, 0.8494949,
			  0.8626263, 0.8757576, 0.8888889, 0.9020202, 0.9151515,
			  0.9282829, 0.9414141, 0.9545454, 0.9676768, 0.9808081,
			  0.9939394, 1.0070707, 1.0202020, 1.0333333, 1.0464647,
			  1.0595959, 1.0727273, 1.0858586, 1.0989898, 1.1121212,
			  1.1252525, 1.1383839, 1.1515151, 1.1646465, 1.1777778,
			  1.1909091, 1.2040404, 1.2171717, 1.2303030, 1.2434343,
			  1.2565657, 1.2696970, 1.2828283, 1.2959596, 1.3090909,
			  1.3222222, 1.3353535, 1.3484849, 1.3616161, 1.3747475,
			  1.3878788, 1.4010102, 1.4141414, 1.4272727, 1.4404041,
			  1.4535353, 1.4666667, 1.4797980, 1.4929293, 1.5060606,
			  1.5191919, 1.5323232, 1.5454545, 1.5585859, 1.5717171,
			  1.5848485, 1.5979798, 1.6111112, 1.6242424, 1.6373737,
			  1.6505051, 1.6636363, 1.6767677, 1.6898990, 1.7030303,
			  1.7161616, 1.7292930, 1.7424242, 1.7555555, 1.7686869,
			  1.7818182, 1.7949495, 1.8080808, 1.8212122, 1.8343434,
			  1.8474747, 1.8606061, 1.8737373};          
  // ensigma in GeV2
  double   ensigma[100]={ 0.2916000, 0.3206586, 0.3497172, 0.3787757, 0.4078344,
			  0.4368929, 0.4659515, 0.4950101, 0.5240687, 0.5531273,
			  0.5821859, 0.6112444, 0.6403030, 0.6693616, 0.6984202,
			  0.7274788, 0.7565374, 0.7855960, 0.8146545, 0.8437131,
			  0.8727717, 0.9018303, 0.9308889, 0.9599475, 0.9890060,
			  1.0180646, 1.0471232, 1.0761818, 1.1052403, 1.1342990,
			  1.1633576, 1.1924162, 1.2214748, 1.2505333, 1.2795919,
			  1.3086505, 1.3377091, 1.3667676, 1.3958262, 1.4248848,
			  1.4539435, 1.4830021, 1.5120606, 1.5411192, 1.5701778,
			  1.5992364, 1.6282949, 1.6573535, 1.6864121, 1.7154707,
			  1.7445292, 1.7735878, 1.8026465, 1.8317051, 1.8607637,
			  1.8898222, 1.9188808, 1.9479394, 1.9769980, 2.0060565,
			  2.0351152, 2.0641737, 2.0932324, 2.1222908, 2.1513495,
			  2.1804080, 2.2094667, 2.2385252, 2.2675838, 2.2966425,
			  2.3257010, 2.3547597, 2.3838181, 2.4128768, 2.4419353,
			  2.4709940, 2.5000525, 2.5291111, 2.5581696, 2.5872283,
			  2.6162868, 2.6453454, 2.6744041, 2.7034626, 2.7325213,
			  2.7615798, 2.7906384, 2.8196969, 2.8487556, 2.8778141,
			  2.9068727, 2.9359312, 2.9649899, 2.9940486, 3.0231071,
			  3.0521657, 3.0812242, 3.1102829, 3.1393414, 3.1684000};
  double    Fsigma[100]={ 2.0261996, 2.2349865, 2.4839740, 2.7840748, 3.1488798,
			  3.5936222, 4.1301847, 4.7517977, 5.3984838, 5.9147439,
			  6.0864558, 5.8283591, 5.2841811, 4.6615186, 4.0839195,
			  3.5914702, 3.1841860, 2.8494759, 2.5732665, 2.3434010,
			  2.1502059, 1.9862038, 1.8456544, 1.7241427, 1.6182493,
			  1.5253036, 1.4432002, 1.3702650, 1.3051554, 1.2467849,
			  1.1942677, 1.1468738, 1.1039963, 1.0651271, 1.0298390,
			  0.9977714, 0.9686196, 0.9421255, 0.9180685, 0.8962603,
			  0.8765374, 0.8587573, 0.8427927, 0.8285285, 0.8158574,
			  0.8046767, 0.7948853, 0.7863811, 0.7790571, 0.7728010,
			  0.7674922, 0.7630011, 0.7591889, 0.7559078, 0.7530031,
			  0.7503147, 0.7476809, 0.7449428, 0.7419487, 0.7385587,
			  0.7346500, 0.7301207, 0.7248930, 0.7189151, 0.7121620,
			  0.7046344, 0.6963565, 0.6873729, 0.6777444, 0.6675445,
			  0.6568548, 0.6457604, 0.6343476, 0.6227004, 0.6108983,
			  0.5990148, 0.5871165, 0.5752623, 0.5635037, 0.5518846,
			  0.5404415, 0.5292045, 0.5181981, 0.5074410, 0.4969472,
			  0.4867267, 0.4767860, 0.4671288, 0.4577557, 0.4486661,
			  0.4398569, 0.4313242, 0.4230627, 0.4150662, 0.4073282,
			  0.3998415, 0.3925985, 0.3855914, 0.3788125, 0.3722538};
  // set up the interpolators
  _Fomega  = make_InterpolatorPtr( 98,Fomegainit,1.0,eninit, GeV, 1);
  _Fthreec = make_InterpolatorPtr( 98,Fthreeinit,1.0,eninit, GeV, 1);
  _Fonec   = make_InterpolatorPtr( 98,Foneinit  ,1.0,eninit, GeV, 1);
  _Fsigma  = make_InterpolatorPtr(100,Fsigma    ,1.0,ensigma,GeV2,1);
  // initialise the calculation of the a_1 width
  _initializea1=false;
  // in GeV2
  double  a1q2in[200]={0,15788.6,31577.3,47365.9,63154.6,78943.2,94731.9,110521,
		       126309,142098,157886,173675,189464,205252,221041,236830,252618,
		       268407,284196,299984,315773,331562,347350,363139,378927,394716,
		       410505,426293,442082,457871,473659,489448,505237,521025,536814,
		       552603,568391,584180,599969,615757,631546,647334,663123,678912,
		       694700,710489,726278,742066,757855,773644,789432,805221,821010,
		       836798,852587,868375,884164,899953,915741,931530,947319,963107,
		       978896,994685,1.01047e+06,1.02626e+06,1.04205e+06,1.05784e+06,
		       1.07363e+06,1.08942e+06,1.10521e+06,1.12099e+06,1.13678e+06,
		       1.15257e+06,1.16836e+06,1.18415e+06,1.19994e+06,1.21573e+06,
		       1.23151e+06,1.2473e+06,1.26309e+06,1.27888e+06,1.29467e+06,
		       1.31046e+06,1.32625e+06,1.34203e+06,1.35782e+06,1.37361e+06,
		       1.3894e+06,1.40519e+06,1.42098e+06,1.43677e+06,1.45256e+06,
		       1.46834e+06,1.48413e+06,1.49992e+06,1.51571e+06,1.5315e+06,
		       1.54729e+06,1.56308e+06,1.57886e+06,1.59465e+06,1.61044e+06,
		       1.62623e+06,1.64202e+06,1.65781e+06,1.6736e+06,1.68939e+06,
		       1.70517e+06,1.72096e+06,1.73675e+06,1.75254e+06,1.76833e+06,
		       1.78412e+06,1.79991e+06,1.81569e+06,1.83148e+06,1.84727e+06,
		       1.86306e+06,1.87885e+06,1.89464e+06,1.91043e+06,1.92621e+06,
		       1.942e+06,1.95779e+06,1.97358e+06,1.98937e+06,2.00516e+06,
		       2.02095e+06,2.03674e+06,2.05252e+06,2.06831e+06,2.0841e+06,
		       2.09989e+06,2.11568e+06,2.13147e+06,2.14726e+06,2.16304e+06,
		       2.17883e+06,2.19462e+06,2.21041e+06,2.2262e+06,2.24199e+06,
		       2.25778e+06,2.27356e+06,2.28935e+06,2.30514e+06,2.32093e+06,
		       2.33672e+06,2.35251e+06,2.3683e+06,2.38409e+06,2.39987e+06,
		       2.41566e+06,2.43145e+06,2.44724e+06,2.46303e+06,2.47882e+06,
		       2.49461e+06,2.51039e+06,2.52618e+06,2.54197e+06,2.55776e+06,
		       2.57355e+06,2.58934e+06,2.60513e+06,2.62092e+06,2.6367e+06,
		       2.65249e+06,2.66828e+06,2.68407e+06,2.69986e+06,2.71565e+06,
		       2.73144e+06,2.74722e+06,2.76301e+06,2.7788e+06,2.79459e+06,
		       2.81038e+06,2.82617e+06,2.84196e+06,2.85774e+06,2.87353e+06,
		       2.88932e+06,2.90511e+06,2.9209e+06,2.93669e+06,2.95248e+06,
		       2.96827e+06,2.98405e+06,2.99984e+06,3.01563e+06,3.03142e+06,
		       3.04721e+06,3.063e+06,3.07879e+06,3.09457e+06,3.11036e+06,
		       3.12615e+06,3.14194e+06};
  // in GeV
  double a1widthin[200]={0,0,0,0,0,0,0,0,
			 0,0,0,0,0.000634625,0.00686721,0.026178,0.066329,
			 0.134996,0.239698,0.387813,0.586641,0.843471,1.16567,
			 1.56076,2.03654,2.60115,3.26324,4.03202,4.91749,
			 5.93053,7.08313,8.38858,9.86176,11.5194,13.3805,
			 15.4667,17.8029,20.4175,23.3438,26.6202,30.2917,
			 34.4108,39.0384,44.2457,50.1143,56.7369,64.2147,
			 72.6566,82.1666,92.8329,104.708,117.786,131.981,
			 147.124,162.974,179.244,195.64,211.904,227.818,
			 243.223,257.991,272.06,285.392,297.971,309.8,
			 320.894,331.278,340.979,350.03,358.463,366.31,
			 373.608,380.386,386.677,392.511,397.945,402.935,
			 407.563,411.841,415.79,419.433,422.766,425.853,
			 428.695,431.302,433.715,435.883,437.887,439.716,
			 441.426,442.947,444.326,445.575,446.65,447.666,
			 448.578,449.395,450.123,450.821,451.343,451.847,
			 452.283,452.859,452.987,453.266,453.496,453.686,
			 453.839,453.958,454.05,454.113,454.149,454.16,
			 454.154,454.13,454.091,454.037,453.966,453.9,
			 453.814,453.724,453.628,453.528,453.417,453.314,
			 453.206,453.1,452.995,452.891,452.79,452.697,
			 452.598,452.509,452.423,452.343,452.269,452.201,
			 452.141,452.086,452.039,452.004,451.966,451.941,
			 451.926,451.888,451.919,451.928,451.945,451.971,
			 452.006,452.05,452.102,452.163,452.234,452.421,
			 452.401,452.498,452.605,452.718,452.84,452.971,
			 453.111,453.261,453.417,453.583,453.756,453.937,
			 454.126,454.324,454.529,455.023,454.964,455.719,
			 455.428,455.671,455.921,456.179,456.444,456.695,
			 456.996,457.276,457.57,457.867,458.171,458.478,
			 458.793,459.115,459.442,459.777,460.115,460.458,
			 460.809,461.161,461.52,461.884,462.253,462.626,
			 463.004,463.832,463.778,464.166};
  vector<double> tmp1(a1widthin,a1widthin+200);
  _a1runwidth.clear();
  std::transform(tmp1.begin(), tmp1.end(),
		 back_inserter(_a1runwidth),
		 [](double x){return x*GeV;});
  
  vector<double> tmp2(a1q2in,a1q2in+200);
  _a1runq2.clear();
  std::transform(tmp2.begin(), tmp2.end(),
		 back_inserter(_a1runq2),
		 [](double x){return x*GeV2;});

  _maxmass=ZERO;
  _maxcalc=ZERO;
}

void FourPionNovosibirskCurrent::doinit() {
  WeakCurrent::doinit();
  // pion masses
  _mpic=getParticleData(ParticleID::piplus)->mass();
  _mpic2=sqr(_mpic);
  _mpi0=getParticleData(ParticleID::pi0)->mass();
  _mpi02=sqr(_mpi0);
  if(!_localparameters) {
    _rhomass    = getParticleData(ParticleID::rhominus)->mass();
    _rhowidth   = getParticleData(ParticleID::rhominus)->width();
    _omegamass  = getParticleData(ParticleID::omega)->mass();
    _omegawidth = getParticleData(ParticleID::omega)->width();
    _sigmamass  = getParticleData(9000221)->mass();
    _sigmawidth = getParticleData(9000221)->width();
    _a1mass    = getParticleData(ParticleID::a_1minus)->mass();
    _a1width   = getParticleData(ParticleID::a_1minus)->width();
  }
  // calculate the constants for the a_1 form factor
  _onedlam2 = 1./_lambda2;
  _a1massolam2 = _a1mass*_a1mass*_onedlam2;
  // parameter for the sigma breit-wigner
  _psigma.push_back(Kinematics::pstarTwoBodyDecay(_sigmamass,_mpi0,_mpi0));
  _psigma.push_back(Kinematics::pstarTwoBodyDecay(_sigmamass,_mpic,_mpic));
  // parameters for the rho breit wigner
  _prho=Kinematics::pstarTwoBodyDecay(_rhomass,_mpic,_mpic);
  _hm2 = hFunction(_rhomass);
  _dhdq2m2=dhdq2Parameter();
  _rhoD=DParameter();
  // convert the magnitude and phase of z into a phase
  _zsigma = _zmag*(cos(_zphase)+Complex(0.,1.)*sin(_zphase));
  // initialize the a_1 width
  inita1width(-1);
}

void FourPionNovosibirskCurrent::doinitrun() {
  // set up the running a_1 width
  inita1width(0);
  WeakCurrent::doinitrun();
}

void FourPionNovosibirskCurrent::persistentOutput(PersistentOStream & os) const {
  os << _a1runinter << _Fomega << _Fthreec << _Fonec << _Fsigma
     << ounit(_rhomass,GeV) << ounit(_a1mass,GeV) << ounit(_omegamass,GeV) 
     << ounit(_sigmamass,GeV) << ounit(_rhowidth,GeV) << ounit(_a1width,GeV)
     << ounit(_omegawidth,GeV) << ounit(_sigmawidth,GeV) 
     << _zsigma << ounit(_lambda2,GeV2)
     << _initializea1 << _localparameters 
     << ounit(_a1runwidth,GeV) << ounit(_a1runq2,GeV2) << ounit(_onedlam2,1/GeV2) 
     << _a1massolam2 << ounit(_psigma,GeV) << ounit(_mpic,GeV) << ounit(_mpi0,GeV)
     << ounit(_aomega,1/GeV) << ounit(_athreec,1/GeV) << ounit(_aonec,1/GeV) 
     << _bomega << _bthreec << _bonec 
     << _comega << _cthreec <<_conec << _omegaparam 
     << ounit(_intwidth,GeV) << ounit(_intmass,GeV)
     << ounit(_mpic2,GeV2) << ounit(_mpi02,GeV2) << ounit(_hm2,GeV2) << _dhdq2m2 
     << ounit(_prho,GeV) << ounit(_rhoD,GeV2) << _zmag << _zphase
     << ounit(_maxmass,GeV) << ounit(_maxcalc,GeV);
}

void FourPionNovosibirskCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _a1runinter >> _Fomega >> _Fthreec >> _Fonec >> _Fsigma
     >> iunit(_rhomass,GeV) >> iunit(_a1mass,GeV) >> iunit(_omegamass,GeV) 
     >> iunit(_sigmamass,GeV) >> iunit(_rhowidth,GeV) >> iunit(_a1width,GeV)
     >> iunit(_omegawidth,GeV) >> iunit(_sigmawidth,GeV) 
     >> _zsigma >> iunit(_lambda2,GeV2)
     >> _initializea1 >> _localparameters 
     >> iunit(_a1runwidth,GeV) >> iunit(_a1runq2,GeV2) >> iunit(_onedlam2,1/GeV2)
     >> _a1massolam2 >> iunit(_psigma,GeV) >> iunit(_mpic,GeV) >> iunit(_mpi0,GeV)
     >> iunit(_aomega,1/GeV) >> iunit(_athreec,1/GeV) >> iunit(_aonec,1/GeV) 
     >> _bomega >> _bthreec >> _bonec 
     >> _comega >> _cthreec >>_conec >> _omegaparam 
     >> iunit(_intwidth,GeV) >> iunit(_intmass,GeV)
     >> iunit(_mpic2,GeV2) >> iunit(_mpi02,GeV2)>> iunit(_hm2,GeV2) >> _dhdq2m2
     >> iunit(_prho,GeV) >> iunit(_rhoD,GeV2) >> _zmag >> _zphase
     >> iunit(_maxmass,GeV) >> iunit(_maxcalc,GeV);
}

// Definition of the static class description member.

void FourPionNovosibirskCurrent::Init() {

  static ClassDocumentation<FourPionNovosibirskCurrent> documentation
    ("The FourPionNovosibirskCurrent class performs the decay"
     " of the tau to four pions using currents based on the the"
     " Novosibirsk e+e- data",
     "The decay of the tau to four pions uses currents based on \\cite{Bondar:2002mw}.",
     "%\\cite{Bondar:2002mw}\n"
     "\\bibitem{Bondar:2002mw}\n"
     "  A.~E.~Bondar, S.~I.~Eidelman, A.~I.~Milstein, T.~Pierzchala, N.~I.~Root, Z.~Was and M.~Worek,\n"
     "   ``Novosibirsk hadronic currents for tau --> 4pi channels of tau decay\n"
     "  %library TAUOLA,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 146}, 139 (2002)\n"
     "  [arXiv:hep-ph/0201149].\n"
     "  %%CITATION = CPHCB,146,139;%%\n"
     );

  static Parameter<FourPionNovosibirskCurrent,Energy> interfacerhoMass
    ("rhoMass",
     "The local value of the rho mass",
     &FourPionNovosibirskCurrent::_rhomass, GeV,0.7761*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfacea1mass
    ("a1Mass",
     "The local value of the square of the a_1 mass",
     &FourPionNovosibirskCurrent::_a1mass, GeV, 1.2300*GeV, 0.5*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceSigmaMass
    ("sigmaMass",
     "The local value of the sigma mass",
     &FourPionNovosibirskCurrent::_sigmamass, GeV, 0.8*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceOmegaMass
    ("omegaMass",
     "The local value of the omega mass",
     &FourPionNovosibirskCurrent::_omegamass, GeV, 0.7820*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfacerhoWidth
    ("rhoWidth",
     "The local value of the rho width",
     &FourPionNovosibirskCurrent::_rhowidth, GeV,0.1445*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfacea1width
    ("a1Width",
     "The local value of the square of the a_1 width",
     &FourPionNovosibirskCurrent::_a1width, GeV, 0.45*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceSigmaWidth
    ("sigmaWidth",
     "The local value of the sigma width",
     &FourPionNovosibirskCurrent::_sigmawidth, GeV, 0.8*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceOmegaWidth
    ("omegaWidth",
     "The local value of the omega width",
     &FourPionNovosibirskCurrent::_omegawidth, GeV, 0.00841*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceIntegrationMass
    ("IntegrationMass",
     "Mass of the pseudoresonance used to improve integration effciency",
     &FourPionNovosibirskCurrent::_intmass, GeV, 1.4*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceIntegrationWidth
    ("IntegrationWidth",
     "Width of the pseudoresonance used to improve integration effciency",
     &FourPionNovosibirskCurrent::_intwidth, GeV, 0.5*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,double> interfaceSigmaMagnitude
    ("SigmaMagnitude",
     "magnitude of the relative sigma coupling",
     &FourPionNovosibirskCurrent::_zmag, 1.3998721, 0.0, 10.0e20,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,double> interfaceSigmaPhase
    ("SigmaPhase",
     "phase of the relative sigma coupling",
     &FourPionNovosibirskCurrent::_zphase, 0.43585036, 0.0, Constants::twopi,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy2> interfaceLambda2
    ("Lambda2",
     "The value of the mass scale squared to use in the form-factor",
     &FourPionNovosibirskCurrent::_lambda2, GeV2, 1.2*GeV2, 0.0001*GeV2, 10.0*GeV2,
     false, false, true);

  static Switch<FourPionNovosibirskCurrent,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the intermediate resonances masses and widths",
     &FourPionNovosibirskCurrent::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use the local values",
     true);
  static SwitchOption interfaceLocalParametersDefault
    (interfaceLocalParameters,
     "ParticleData",
     "Use the values from the particleData objects",
     false);

  static Switch<FourPionNovosibirskCurrent,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &FourPionNovosibirskCurrent::_initializea1, false, false, false);
  static SwitchOption interfaceInitializea1Initialization
    (interfaceInitializea1,
     "Yes",
     "Initialize the calculation",
     true);
  static SwitchOption interfaceInitializea1NoInitialization
    (interfaceInitializea1,
     "No",
     "Use the default values",
     false);
  
  static ParVector<FourPionNovosibirskCurrent,Energy> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &FourPionNovosibirskCurrent::_a1runwidth, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<FourPionNovosibirskCurrent,Energy2> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &FourPionNovosibirskCurrent::_a1runq2, GeV2, -1, 1.0*GeV2, ZERO, 10.0*GeV2,
     false, false, true);

}

// initialisation of the a_1 running width 
void FourPionNovosibirskCurrent::inita1width(int iopt) {
  // set up the interpolator
  if(iopt==0||!_initializea1) {
    _a1runinter = make_InterpolatorPtr(_a1runwidth,_a1runq2,3);
    return;
  }
  _maxcalc=_maxmass;
  if(_maxmass==ZERO) return;
  // parameters for the table of values
  Energy2 step(sqr(_maxmass)/200.);
  // function to be integrated to give the matrix element
  // integrator to perform the integral
  // weights for the integration channels
  vector<double> inweights;
  inweights.push_back(0.3);inweights.push_back(0.3);inweights.push_back(0.3);
  vector<double> inpower(3, 0.0);
  // types of integration channels
  vector<int> intype;
  intype.push_back(2);intype.push_back(3);intype.push_back(1);
  // masses for the integration channels
  vector<Energy> inmass(2,_rhomass);inmass.push_back(_sigmamass);
  // widths for the integration channels
  vector<Energy> inwidth(2,_rhowidth);inwidth.push_back(_sigmawidth);
  ThreeBodyAllOnCalculator<FourPionNovosibirskCurrent> 
    widthgen1(inweights,intype,inmass,inwidth,inpower,*this,0,_mpi0,_mpic,_mpic); 
  ThreeBodyAllOnCalculator<FourPionNovosibirskCurrent>
    widthgen2(inweights,intype,inmass,inwidth,inpower,*this,1,_mpi0,_mpi0,_mpi0); 
  // normalisation constant to give physical width if on shell
  double a1const(_a1width/(widthgen1.partialWidth(sqr(_a1mass))+
			   widthgen2.partialWidth(sqr(_a1mass))));
  // loop to give the values
  Energy2 moff2(ZERO);
  _a1runwidth.clear();_a1runq2.clear();
  for(;moff2<=sqr(_maxmass);moff2+=step) {
    Energy total = a1const*(widthgen1.partialWidth(moff2)+widthgen2.partialWidth(moff2));
    _a1runwidth.push_back(total);
    _a1runq2.push_back(moff2);
  }
}

// complete the construction of the decay mode for integration
bool FourPionNovosibirskCurrent::createMode(int icharge, tcPDPtr resonance,
					    FlavourInfo flavour,
					    unsigned int imode,PhaseSpaceModePtr mode,
					    unsigned int iloc,int ires,
					    PhaseSpaceChannel phase, Energy upp ) {
  if(resonance) return false;
  // check the charge
  if(abs(icharge)!=3) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      return false;
      break;
    case IsoSpin::I3One:
      if(icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
  // check that the modes are kinematical allowed
  Energy min(ZERO);
  if(imode==0) {
    min=   getParticleData(ParticleID::piplus)->mass()
      +3.*getParticleData(ParticleID::pi0)->mass();
  }
  else {
    min=3.*getParticleData(ParticleID::piplus)->mass()
      +getParticleData(ParticleID::pi0)->mass();
  }
  if(min>upp) return false;
  _maxmass=max(upp,_maxmass);
  // intermediates for the channels
  tPDPtr omega(getParticleData(ParticleID::omega)),rhop,rhom,
    rho0(getParticleData(ParticleID::rho0)),a1m,a10(getParticleData(ParticleID::a_10)),
    sigma(getParticleData(9000221)),rhot;
  if(icharge==3) {
    rhop = getParticleData(ParticleID::rhominus);
    rhom = getParticleData(ParticleID::rhoplus);
    a1m  = getParticleData(ParticleID::a_1plus);
    rhot = getParticleData(24);
  }
  else {
    rhop = getParticleData(ParticleID::rhoplus);
    rhom = getParticleData(ParticleID::rhominus);
    a1m  = getParticleData(ParticleID::a_1minus);
    rhot = getParticleData(-24);
  }
  if(imode==1) {
    // the omega channels for the three charged pion mode
    // first  channel two channels with rho0
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,omega,ires+1,iloc+1,
		      ires+2,rho0,ires+2,iloc+4,ires+3,iloc+2,ires+3,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,omega,ires+1,iloc+2,
		      ires+2,rho0,ires+2,iloc+4,ires+3,iloc+1,ires+3,iloc+3));
    // second two channels with rho -
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,omega,ires+1,iloc+1,
		      ires+2,rhom,ires+2,iloc+3,ires+3,iloc+2,ires+3,iloc+4));
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,omega,ires+1,iloc+2,
		      ires+2,rhom,ires+2,iloc+3,ires+3,iloc+1,ires+3,iloc+4));
    // third two channels with rho +
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,omega,ires+1,iloc+1,
		      ires+2,rhop,ires+2,iloc+2,ires+3,iloc+3,ires+3,iloc+4));
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,omega,ires+1,iloc+2,
		      ires+2,rhop,ires+2,iloc+1,ires+3,iloc+3,ires+3,iloc+4));
    //  a_1 channels with rhos
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a1m,ires+1,iloc+4,
		      ires+2,rho0,ires+2,iloc+1,ires+3,iloc+2,ires+3,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a1m,ires+1,iloc+4,
		      ires+2,rho0,ires+2,iloc+2,ires+3,iloc+1,ires+3,iloc+3));
    // neutral a_1 channels with rhos
    // rho-
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a10,ires+1,iloc+1,
		      ires+2,rhom,ires+2,iloc+3,ires+3,iloc+2,ires+3,iloc+4));
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a10,ires+1,iloc+2,
		      ires+2,rhom,ires+2,iloc+3,ires+3,iloc+1,ires+3,iloc+4));
    // rho+
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a10,ires+1,iloc+1,
		      ires+2,rhop,ires+2,iloc+2,ires+3,iloc+3,ires+3,iloc+4));
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a10,ires+1,iloc+2,
		      ires+2,rhop,ires+2,iloc+1,ires+3,iloc+3,ires+3,iloc+4));
    //  a_1 channels with sigmas
    if(sigma) {
      // charged a_1 channels with sigma
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a1m,ires+1,iloc+4,
			ires+2,sigma,ires+2,iloc+1,ires+3,iloc+2,ires+3,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a1m,ires+1,iloc+4,
			ires+2,sigma,ires+2,iloc+2,ires+3,iloc+1,ires+3,iloc+3));
      // neutral a_1 channels with sigma
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a10,ires+1,iloc+1,
			ires+2,sigma,ires+2,iloc+4,ires+3,iloc+2,ires+3,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a10,ires+1,iloc+2,
			ires+2,sigma,ires+2,iloc+4,ires+3,iloc+1,ires+3,iloc+3));
    }
  }
  else {
  //   // channels with an a1- and a rho -
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a1m,ires+1,iloc+2,
		      ires+2,rhom,ires+2,iloc+3,ires+3,iloc+1,ires+3,iloc+4));
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a1m,ires+1,iloc+2,
		      ires+2,rhom,ires+2,iloc+4,ires+3,iloc+1,ires+3,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a1m,ires+1,iloc+3,
		      ires+2,rhom,ires+2,iloc+2,ires+3,iloc+1,ires+3,iloc+4));
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a1m,ires+1,iloc+3,
		      ires+2,rhom,ires+2,iloc+4,ires+3,iloc+1,ires+3,iloc+2));
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a1m,ires+1,iloc+4,
		      ires+2,rhom,ires+2,iloc+2,ires+3,iloc+1,ires+3,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a1m,ires+1,iloc+4,
		      ires+2,rhom,ires+2,iloc+3,ires+3,iloc+1,ires+3,iloc+2));
    // channels with a sigma and a10
    if(sigma ) {
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a10,ires+1,iloc+1,
			ires+2,sigma,ires+2,iloc+2,ires+3,iloc+3,ires+3,iloc+4));
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a10,ires+1,iloc+1,
			ires+2,sigma,ires+2,iloc+3,ires+3,iloc+2,ires+3,iloc+4));
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a10,ires+1,iloc+1,
			ires+2,sigma,ires+2,iloc+4,ires+3,iloc+2,ires+3,iloc+3));
      // channels with a1- and sigma
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a1m,ires+1,iloc+2,
			ires+2,sigma,ires+2,iloc+1,ires+3,iloc+3,ires+3,iloc+4));
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a1m,ires+1,iloc+3,
			ires+2,sigma,ires+2,iloc+1,ires+3,iloc+2,ires+3,iloc+4));
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhot,ires+1,a1m,ires+1,iloc+4,
			ires+2,sigma,ires+2,iloc+1,ires+3,iloc+2,ires+3,iloc+3));
    }
  }
  // reset the parameters of the dummy resonance used for integration
  mode->resetIntermediate(rhot,_intmass,_intwidth);
  // reset the parameters of the resonances if using local values
  if(_localparameters) {
    mode->resetIntermediate(rhom,_rhomass,_rhowidth);
    mode->resetIntermediate(rhop,_rhomass,_rhowidth);
    mode->resetIntermediate(rho0,_rhomass,_rhowidth);
    mode->resetIntermediate(omega,_omegamass,_omegawidth);
    if(sigma) mode->resetIntermediate(sigma,_sigmamass,_sigmawidth);
  }
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector FourPionNovosibirskCurrent::particles(int icharge, unsigned int imode,int,int) {
  if(abs(icharge)!=3) return tPDVector();
  tPDVector output(4);
  if(imode==1) {
    output[0]=getParticleData(ParticleID::piplus);
    output[1]=getParticleData(ParticleID::piplus);
    output[2]=getParticleData(ParticleID::piminus);
  }
  else {
    output[0]=getParticleData(ParticleID::piplus);
    output[1]=getParticleData(ParticleID::pi0);
    output[2]=getParticleData(ParticleID::pi0);
  }
  output[3]=getParticleData(ParticleID::pi0);
  if(icharge==-3) {
    for(unsigned int ix=0;ix<output.size();++ix) {
      if(output[ix]->CC()) output[ix]=output[ix]->CC();
    }
  }
  return output;
}

 
// the hadronic currents    
vector<LorentzPolarizationVectorE> 
FourPionNovosibirskCurrent::current(tcPDPtr resonance,
				    FlavourInfo flavour,
				    const int imode, const int ichan,Energy & scale,
				    const tPDVector & outgoing,
				    const vector<Lorentz5Momentum> & momenta,
				    DecayIntegrator::MEOption) const {
  assert(!resonance);
  int icharge(0);
  for(tPDPtr out : outgoing) icharge+=out->iCharge();
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return vector<LorentzPolarizationVectorE>();
  useMe();
  LorentzVector<complex<InvEnergy> > output;
  double fact(1.);
  // the momenta of the particles
  Lorentz5Momentum q1(momenta[0]),q2(momenta[2]),
    q3(momenta[1]),q4(momenta[3]);
  Lorentz5Momentum Q(q1+q2+q3+q4);Q.rescaleMass();
  scale = Q.mass();
  // decide which decay mode
  // three charged pions
  if(imode==1) {
    // momenta of the particles
    LorentzVector<complex<Energy5> > veca1rho,vecomega,veca1sig;
    if(ichan<0) {
      // a_1 rho current
      veca1rho = 
  	t1(q1,q2,q3,q4)+t1(q3,q2,q1,q4)+t1(q1,q3,q2,q4)
  	+t1(q3,q1,q2,q4)+t1(q4,q3,q1,q2)+t1(q4,q1,q3,q2);
      // a_1 sigma current
      veca1sig = 
  	t2(q4,q3,q1,q2,1)+t2(q4,q1,q3,q2,1)
  	-t2(q1,q4,q3,q2,1)-t2(q3,q4,q1,q2,1);
      // omega current
      vecomega = 
  	t3(q1,q2,q3,q4)+t3(q3,q2,q1,q4)-t3(q1,q3,q2,q4)
  	-t3(q3,q1,q2,q4)-t3(q1,q4,q3,q2)-t3(q3,q4,q1,q2);
    }
    else if(ichan== 0) vecomega = t3(q1,q4,q3,q2);
    else if(ichan== 1) vecomega = t3(q3,q4,q1,q2);
    else if(ichan== 2) vecomega = t3(q1,q2,q3,q4);
    else if(ichan== 3) vecomega = t3(q3,q2,q1,q4);
    else if(ichan== 4) vecomega = t3(q1,q3,q2,q4);
    else if(ichan== 5) vecomega = t3(q3,q1,q2,q4);
    else if(ichan== 6) veca1rho = t1(q4,q1,q3,q2);
    else if(ichan== 7) veca1rho = t1(q4,q3,q1,q2);
    else if(ichan== 8) veca1rho = t1(q1,q2,q3,q4);
    else if(ichan== 9) veca1rho = t1(q3,q2,q1,q4);
    else if(ichan==10) veca1rho = t1(q1,q3,q2,q4);
    else if(ichan==11) veca1rho = t1(q3,q1,q2,q4);
    else if(ichan==12) veca1sig = t2(q4,q1,q3,q2,1);
    else if(ichan==13) veca1sig = t2(q4,q3,q1,q2,1);
    else if(ichan==14) veca1sig = t2(q1,q4,q3,q2,1);
    else if(ichan==15) veca1sig = t2(q3,q4,q1,q2,1);
    // final manipulations
    veca1rho += veca1sig;
    LorentzVector<complex<InvEnergy> > 
      veca1rho1 = veca1rho * gFunction(Q.mass2(),1);
    LorentzVector<complex<InvEnergy> > 
      vecomega1 = vecomega * gFunction(Q.mass2(),2);
    output = vecomega1 + veca1rho1;
    // this is 1/sqrt(2) for identical particles
    fact *= 1./sqrt(2.);
  }
  else if(imode==0) {
    // momenta of the particles
    LorentzVector<complex<Energy5> > veca1rho,veca1sig;
    if(ichan<0) {
      // a_1 rho current
      veca1rho= t1(q2,q3,q1,q4)+t1(q2,q4,q1,q3)+t1(q3,q2,q1,q4)
  	+t1(q3,q4,q1,q2)+t1(q4,q2,q1,q3)+t1(q4,q3,q1,q2);
      // a_1 sigma current
      veca1sig=
  	t2(q2,q1,q3,q4,0)+t2(q3,q1,q2,q4,0)+t2(q4,q1,q3,q2,0)
  	-t2(q1,q2,q3,q4,0)-t2(q1,q3,q2,q4,0)-t2(q1,q4,q3,q2,0);
    }
    else if(ichan== 0) veca1rho = t1(q2,q3,q1,q4);
    else if(ichan== 1) veca1rho = t1(q2,q4,q1,q3);
    else if(ichan== 2) veca1rho = t1(q3,q2,q1,q4);
    else if(ichan== 3) veca1rho = t1(q3,q4,q1,q2);
    else if(ichan== 4) veca1rho = t1(q4,q2,q1,q3);
    else if(ichan== 5) veca1rho = t1(q4,q3,q1,q2);
    else if(ichan== 6) veca1sig = t2(q2,q1,q3,q4,0);
    else if(ichan== 7) veca1sig = t2(q3,q1,q2,q4,0);
    else if(ichan== 8) veca1sig = t2(q4,q1,q3,q2,0);
    else if(ichan== 9) veca1sig = t2(q1,q2,q3,q4,0);
    else if(ichan==10) veca1sig = t2(q1,q3,q2,q4,0);
    else if(ichan==11) veca1sig = t2(q1,q4,q3,q2,0);
    // add them up 
    output = (veca1rho + veca1sig) * gFunction(Q.mass2(),0);
    // this is sqrt(1/3!) for identical particles
    fact *= 1./sqrt(6.);
  }     
  else
    assert(false);
  return vector<LorentzPolarizationVectorE>(1, output * fact * Q.mass2()); 
}

bool FourPionNovosibirskCurrent::accept(vector<int> id) {
  bool allowed(false);
  // check four products
  if(id.size()!=4){return false;}
  int npiminus=0,npiplus=0,npi0=0;
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID:: piplus)      ++npiplus;
    else if(id[ix]==ParticleID::piminus) ++npiminus;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
  }
  if(npiminus==2&&npiplus==1&&npi0==1)      allowed=true;
  else if(npiminus==1&&npi0==3)             allowed=true;
  else if(npiplus==2&&npiminus==1&&npi0==1) allowed=true;
  else if(npiplus==1&&npi0==3)              allowed=true;
  return allowed;
}

// the decay mode
unsigned int FourPionNovosibirskCurrent::decayMode(vector<int> idout) {
  unsigned int npi(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::piplus) ++npi;
  }
  if(npi==3) return 1;
  return 0;
}


// output the information for the database
void FourPionNovosibirskCurrent::dataBaseOutput(ofstream & output,bool header,
						bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::FourPionNovosibirskCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  output << "newdef " << name() << ":rhoMass "    << _rhomass/GeV << "\n";
  output << "newdef " << name() << ":a1Mass  "    << _a1mass/GeV  << "\n";
  output << "newdef " << name() << ":sigmaMass  " << _sigmamass/GeV  << "\n";
  output << "newdef " << name() << ":omegaMass  " << _omegamass/GeV  << "\n";
  output << "newdef " << name() << ":rhoWidth "    << _rhowidth/GeV << "\n";
  output << "newdef " << name() << ":a1Width  "    << _a1width/GeV  << "\n";
  output << "newdef " << name() << ":sigmaWidth  " << _sigmawidth/GeV  << "\n";
  output << "newdef " << name() << ":omegaWidth  " << _omegawidth/GeV  << "\n";
  output << "newdef " << name() << ":IntegrationMass "  << _intmass/GeV  << "\n";
  output << "newdef " << name() << ":IntegrationWidth " << _intwidth/GeV  << "\n";
  output << "newdef " << name() << ":SigmaMagnitude "  <<  _zmag << "\n";
  output << "newdef " << name() << ":SigmaPhase " << _zphase  << "\n";
  output << "newdef " << name() << ":Lambda2 "  <<  _lambda2/GeV2 << "\n";
  output << "newdef " << name() << ":LocalParameters " <<  _localparameters << "\n";
  output << "newdef " << name() << ":Initializea1 " <<  _initializea1 << "\n";
  for(unsigned int ix=0;ix<_a1runwidth.size();++ix) {
    if(ix<200) output << "newdef ";
    else       output << "insert ";
    output << name() << ":a1RunningWidth " << ix << " " 
	   << _a1runwidth[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_a1runq2.size();++ix) {
    if(ix<200) output << "newdef ";
    else       output << "insert ";
    output << name() << ":a1RunningQ2 " << ix << " " << _a1runq2[ix]/GeV2 << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

double FourPionNovosibirskCurrent::
threeBodyMatrixElement(const int iopt, const Energy2 q2,
		       const Energy2 s3, const Energy2 s2, 
		       const Energy2 s1, const Energy,
		       const Energy, const Energy) const {
  unsigned int ix;
  // construct the momenta of the decay products
  Energy p1[5],p2[5],p3[5];
  Energy2 p1sq, p2sq, p3sq;
  Energy q(sqrt(q2));
  if(iopt==0) {
    p1[0] = 0.5*(q2+_mpi02-s1)/q; p1sq=p1[0]*p1[0]; p1[4]=sqrt(p1sq-_mpi02);
    p2[0] = 0.5*(q2+_mpic2-s2)/q; p2sq=p2[0]*p2[0]; p2[4]=sqrt(p2sq-_mpic2);
    p3[0] = 0.5*(q2+_mpic2-s3)/q; p3sq=p3[0]*p3[0]; p3[4]=sqrt(p3sq-_mpic2);
  }
  else {
    p1[0] = 0.5*(q2+_mpi02-s1)/q; p1sq=p1[0]*p1[0]; p1[4]=sqrt(p1sq-_mpi02);
    p2[0] = 0.5*(q2+_mpi02-s2)/q; p2sq=p2[0]*p2[0]; p2[4]=sqrt(p2sq-_mpi02);
    p3[0] = 0.5*(q2+_mpi02-s3)/q; p3sq=p3[0]*p3[0]; p3[4]=sqrt(p3sq-_mpi02);
  }
  // take momentum of 1 parallel to z axis
  p1[1]=ZERO;p1[2]=ZERO;p1[3]=p1[4];
  // construct 2 
  double cos2(0.5*(sqr(p1[4])+sqr(p2[4])-sqr(p3[4]))/p1[4]/p2[4]);
  p2[1] = p2[4]*sqrt(1.-sqr(cos2)); p2[2]=ZERO; p2[3]=-p2[4]*cos2;
  // construct 3
  double cos3(0.5*(sqr(p1[4])-sqr(p2[4])+sqr(p3[4]))/p1[4]/p3[4]);
  p3[1] =-p3[4]*sqrt(1.-sqr(cos3)); p3[2]=ZERO; p3[3]=-p3[4]*cos3;
  // pi+pi-pi0 term
  complex<Energy4> output(0.*sqr(MeV2));
  if(iopt==0) {
    // values for the different Breit-Wigner terms
    Complex rho1(2.365*rhoBreitWigner(s2)),
      rho2(2.365*rhoBreitWigner(s3)),
      sig1(sigmaBreitWigner(s1,1));
    // compute the vector
    complex<Energy2> term;
    for(ix=1;ix<4;++ix) { 
      term = (p1[0]*p2[ix]-p2[0]*p1[ix])*rho2+(p1[0]*p3[ix]-p3[0]*p1[ix])*rho1
	+_zsigma*q*p1[ix]*sig1;
      output+=term*conj(term);
    }
  }
  // pi0pi0pi0 term
  else if(iopt==1) {
    // values for the different Breit-Wigner terms
    Complex sig1(sigmaBreitWigner(s1,0)),
      sig2(sigmaBreitWigner(s2,0)),
      sig3(sigmaBreitWigner(s3,0));
    // compute the vector
    complex<Energy2> term;
    for(ix=1;ix<4;++ix) {
      term = _zsigma * q * (p1[ix]*sig1 + p2[ix]*sig2 + p3[ix]*sig3);
      output += term*conj(term);
    }
    output/=6.;
  }
  output *= a1FormFactor(q2);
  return output.real() / pow<4,1>(_rhomass);
}

Complex FourPionNovosibirskCurrent::sigmaBreitWigner(Energy2 q2,
						     unsigned int iopt) const {
  Energy q(sqrt(q2));
  Energy pcm = iopt==0 ? 
    Kinematics::pstarTwoBodyDecay(q,_mpi0,_mpi0) :
    Kinematics::pstarTwoBodyDecay(q,_mpic,_mpic);
  if(pcm<ZERO) pcm=ZERO;
  Energy  width(_sigmawidth*pcm/_psigma[iopt]);
  Energy2 msigma2 = sqr(_sigmamass);
  return msigma2/(q2-msigma2+Complex(0.,1.)*msigma2*width/q);
}

Complex FourPionNovosibirskCurrent::a1BreitWigner(Energy2 q2) const {
  Complex ii(0.,1.);
  Energy2 m2 = sqr(_a1mass);
  Energy q = sqrt(q2);
  return (m2/complex<Energy2>(q2 - m2 + ii*q*a1width(q2)));
}

Complex FourPionNovosibirskCurrent::omegaBreitWigner(Energy2 q2) const {
  Energy q(sqrt(q2));
  // calcluate the running width
  double diff((q-_omegamass)/GeV),temp(diff);
  double gomega(1.);
  Complex ii(0.,1.);
  if(q<=1.*GeV) {
    for(unsigned int ix=0;ix<6;++ix) {
      gomega +=temp*_omegaparam[ix];
      temp*=diff;
    }
  }
  else {
    gomega=_omegaparam[6]+q/GeV*(_omegaparam[7]+q/GeV*_omegaparam[8]
				 +q2/GeV2*_omegaparam[9]);
  }
  if(gomega<0.){gomega=0.;}
  Energy2 numer=_omegamass*_omegamass;
  complex<Energy2> denom=q2-_omegamass*_omegamass+ii*_omegamass*_omegawidth*gomega;
  return numer/denom;
}

Complex FourPionNovosibirskCurrent::rhoBreitWigner(Energy2 q2) const {
  Energy q(sqrt(q2));
  Energy2 grhom(8.*_prho*_prho*_prho/_rhomass);
  complex<Energy2> denom;
  Complex ii(0.,1.);
  if(q2<4.*_mpic2) {
    denom=q2-_rhomass*_rhomass
      -_rhowidth*_rhomass*(hFunction(q)-_hm2-(q2-_rhomass*_rhomass)*_dhdq2m2)/grhom;
  }
  else {
    Energy pcm(2.*Kinematics::pstarTwoBodyDecay(q,_mpic,_mpic));
    Energy2 grho(pcm*pcm*pcm/q);
    denom=q2-_rhomass*_rhomass
      -_rhowidth*_rhomass*(hFunction(q)-_hm2-(q2-_rhomass*_rhomass)*_dhdq2m2)/grhom
      +ii*_rhomass*_rhowidth*grho/grhom;
  }
  return _rhoD/denom;
}

LorentzVector<complex<Energy5> > 
FourPionNovosibirskCurrent::t1(Lorentz5Momentum & q1,Lorentz5Momentum & q2,
			       Lorentz5Momentum & q3,Lorentz5Momentum & q4) const {
  // momentum of the whole system
  Lorentz5Momentum Q(q1+q2+q3+q4);Q.rescaleMass();
  // compute the virtuality of the a_1
  Lorentz5Momentum a1(q2+q3+q4);a1.rescaleMass();
  // compute the virtuality of the  rho
  Lorentz5Momentum rho(q3+q4);rho.rescaleMass();
  // compute the prefactor    
  Complex pre(-a1FormFactor(a1.mass2())*a1BreitWigner(a1.mass2())*
	      rhoBreitWigner(rho.mass2()));
  // dot products we need
  Energy2 QdQmq1(Q*a1);
  complex<Energy4> consta(QdQmq1*(a1*q3)), constb(QdQmq1*(a1*q4)),
    constc(((Q*q4)*(q1*q3)-(Q*q3)*(q1*q4)));
  // compute the current
  return pre*(consta*q4-constb*q3+constc*a1);
}

LorentzVector<complex<Energy5> > 
FourPionNovosibirskCurrent::t2(Lorentz5Momentum & q1,Lorentz5Momentum & q2,
			       Lorentz5Momentum & q3,Lorentz5Momentum & q4,
			       unsigned int iopt) const {
  // momentum of the whole system
  Lorentz5Momentum Q(q1+q2+q3+q4);Q.rescaleMass();
  // compute the virtuality of the a_1
  Lorentz5Momentum a1(q2+q3+q4);a1.rescaleMass();
  // compute the virtuality of the  sigma
  Lorentz5Momentum sigma(q3+q4);sigma.rescaleMass();
  // compute the prefactor
  Complex pre(_zsigma*a1FormFactor(a1.mass2())
	      *a1BreitWigner(a1.mass2())*
	      sigmaBreitWigner(sigma.mass2(),iopt));
  // dot products we need
  complex<Energy4> consta((Q*a1)*a1.mass2()),constb((Q*q2)*a1.mass2());
  // compute the current
  return pre*(consta*q2-constb*a1);
}

LorentzVector<complex<Energy5> > 
FourPionNovosibirskCurrent::t3(Lorentz5Momentum & q1,Lorentz5Momentum & q2,
			       Lorentz5Momentum & q3,Lorentz5Momentum & q4) const {
  // momentum of the whole sysytem
  Lorentz5Momentum Q(q1+q2+q3+q4);Q.rescaleMass();
  // compute the virtuality of the omega
  Lorentz5Momentum omega(q2+q3+q4);omega.rescaleMass();
  // compute the virtuality of the  rho
  Lorentz5Momentum rho(q3+q4);rho.rescaleMass();
  // compute the prefactor
  Complex pre(omegaBreitWigner(omega.mass2())*rhoBreitWigner(rho.mass2()));
  // dot products we need
  complex<Energy4> consta((Q*q3)*(q1*q4)-(Q*q4)*(q1*q3)),
    constb(-(Q*q2)*(q1*q4)+(q1*q2)*(Q*q4)),
    constc((Q*q2)*(q1*q3)-(q1*q2)*(Q*q3));
  // compute the current
  return pre*(consta*q2+constb*q3+constc*q4);
}

InvEnergy6 FourPionNovosibirskCurrent::gFunction(Energy2 q2, int ichan) const {
  Energy q(sqrt(q2));
  InvEnergy4 invmrho4 = 1/sqr(sqr(_rhomass));
  // the one charged pion G function
  if(ichan==0) {
    return (*_Fonec)(q) * _aonec * (*_Fsigma)(q2) * sqrt(_bonec*q/GeV-_conec) *
      invmrho4/q;
  }
  // the three charged pion G function
  else if(ichan==1) {
    return (*_Fthreec)(q)*_athreec*sqrt(_bthreec*q/GeV-_cthreec)*invmrho4/q; 
  }
  // the omega G function
  else if(ichan==2) {
    return(*_Fomega)(q)*_aomega*sqrt(_bomega*q/GeV-_comega)*invmrho4/q;
  }
  assert(false);
  return InvEnergy6();
}

Energy2 FourPionNovosibirskCurrent::DParameter() const {
  Energy2 grhom(8.*_prho*_prho*_prho/_rhomass);
  return _rhomass*_rhomass+_rhowidth*_rhomass*
    (hFunction(ZERO)-_hm2+_rhomass*_rhomass*_dhdq2m2)/grhom;
}

double FourPionNovosibirskCurrent::dhdq2Parameter() const {
  Energy2 mrho2(_rhomass*_rhomass);
  double root(sqrt(1.-4.*_mpic2/mrho2));
  return root/Constants::pi*(root+(1.+2*_mpic2/mrho2)*log((1+root)/(1-root)));
}

Energy2 FourPionNovosibirskCurrent::hFunction(const Energy q) const {
  using Constants::pi;
  static const Energy2 eps(0.01*MeV2);

  Energy2 q2(q*q), output;
  if (q2 > 4*_mpic2) {
    double root = sqrt(1.-4.*_mpic2/q2);
    output = root*log((1.+root)/(1.-root))*(q2-4*_mpic2)/pi;
  }
  else if (q2 > eps) output = ZERO;
  else               output = -8.*_mpic2/pi;
  return output;
}

#line 1 "./ScalarMesonCurrent.cc"
// -*- C++ -*-
//
// ScalarMesonCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarMesonCurrent class.
//

#include "ScalarMesonCurrent.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void ScalarMesonCurrent::doinit() {
  unsigned int isize=numberOfModes();
  if(_id.size()!=isize||_decay_constant.size()!=isize)
    {throw InitException() << "Inconsistent parameters in ScalarMesonCurrent::doinit()"
			   << Exception::abortnow;}
  WeakCurrent::doinit();
}

ScalarMesonCurrent::ScalarMesonCurrent() {
  // the eta/eta' mixing angle
  _thetaeta=-0.194;
  // the decay constants for the different modes
  _id = {211,111,111,221,221,221,
	 331,331,331,
	 311,321,411,421,431,10431};
  _decay_constant = {130.7*MeV,130.7*MeV,130.7*MeV,130.7*MeV,130.7*MeV,130.7*MeV,
		     130.7*MeV,130.7*MeV,130.7*MeV,
		     159.8*MeV,159.8*MeV,200.0*MeV,200.0*MeV,241.0*MeV,73.7*MeV};
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(3,-3);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(3,-3);
  addDecayMode(1,-3);
  addDecayMode(2,-3);
  addDecayMode(4,-1);
  addDecayMode(4,-2);
  addDecayMode(4,-3);
  addDecayMode(4,-3);
  // initial size of the arrays
  _initsize = _id.size();
  setInitialModes(_initsize);
}

void ScalarMesonCurrent::persistentOutput(PersistentOStream & os) const {
  os << _id << ounit(_decay_constant,GeV) << _thetaeta;
}

void ScalarMesonCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _id >> iunit(_decay_constant,GeV) >> _thetaeta;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarMesonCurrent,WeakCurrent>
describeHerwigScalarMesonCurrent("Herwig::ScalarMesonCurrent", "HwWeakCurrents.so");

void ScalarMesonCurrent::Init() {

  static ClassDocumentation<ScalarMesonCurrent> documentation
    ("The ScalarMesonCurrent class implements the current"
     " for the decay of the weak current into a pseudoscalar meson.");

  static ParVector<ScalarMesonCurrent,long> interfaceID
    ("ID",
     "The PDG code for the outgoing meson.",
     &ScalarMesonCurrent::_id,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<ScalarMesonCurrent,Energy> interfaceDecay_Constant
    ("Decay_Constant",
     "The decay constant for the meson.",
     &ScalarMesonCurrent::_decay_constant, MeV, -1, 100.*MeV,-1000.0*MeV, 1000.0*MeV,
     false, false, true);

  static Parameter<ScalarMesonCurrent,double> interfaceThetaEtaEtaPrime
    ("ThetaEtaEtaPrime",
     "The eta-eta' mixing angle",
     &ScalarMesonCurrent::_thetaeta, -0.194, -Constants::pi, Constants::pi,
     false, false, true);

}

// create the decay phase space mode
bool ScalarMesonCurrent::createMode(int icharge, tcPDPtr resonance,
				    FlavourInfo flavour,
				    unsigned int imode,PhaseSpaceModePtr mode,
				    unsigned int iloc,int ires,
				    PhaseSpaceChannel phase, Energy upp ) {
  assert(!resonance);
  assert(flavour.I==IsoSpin::IUnknown && flavour.I3==IsoSpin::I3Unknown);
  // check the mode has the correct charge
  if(abs(icharge)!=abs(int(getParticleData(_id[imode])->iCharge()))) return false;
  // check if the particle is kinematically allowed
  tPDPtr part(getParticleData(_id[imode]));
  Energy min=part->massMin();
  if(min>upp) return false;
  // construct the mode
  mode->addChannel((PhaseSpaceChannel(phase),ires,iloc+1));
  return true;
}

// outgoing particles 
tPDVector ScalarMesonCurrent::particles(int icharge, unsigned int imode, int iq, int ia) {
  tPDPtr part(getParticleData(_id[imode]));
  tPDVector output;
  if(icharge==int(part->iCharge())) {
    if(icharge==0) {
      int iqb,iab; 
      decayModeInfo(imode,iqb,iab);
      if(iq==iqb&&ia==iab) output.push_back(part);
      else                 output.push_back(part->CC());
    }
    else {
      output.push_back(part);
    }
  }
  else if(icharge==-int(part->iCharge())) {
    output.push_back(part->CC());
  }
  return output;
}

vector<LorentzPolarizationVectorE> 
ScalarMesonCurrent::current(tcPDPtr resonance,
			    FlavourInfo flavour,
			    const int imode, const int , Energy & scale, 
			    const tPDVector & outgoing,
			    const vector<Lorentz5Momentum> & momenta,
			    DecayIntegrator::MEOption) const {
  assert(!resonance);
  assert(flavour.I==IsoSpin::IUnknown && flavour.I3==IsoSpin::I3Unknown);
  static const Complex ii(0.,1.);
  scale =momenta[0].mass();
  Complex pre(-ii*_decay_constant[imode]/scale);
  // quarks in the current
  int iq,ia;
  decayModeInfo(imode,iq,ia);
  if(abs(iq)==abs(ia)) {
    int id(outgoing[0]->id());
    if(id==ParticleID::eta) {
      if(abs(iq)==3) pre*=-2.*cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);
      else           pre*=cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);
    }
    else if(id==ParticleID::etaprime) {
      if(abs(iq)==3) pre*=-2.*sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);
      else           pre*=sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);
    }
    else if(id==ParticleID::pi0&&abs(iq)==1) {
      pre*=-sqrt(0.5);
    }
    else {
      pre*= sqrt(0.5);
    }
  }
  // return the answer
  return vector<LorentzPolarizationVectorE>(1,pre*momenta[0]);
}
  
bool ScalarMesonCurrent::accept(vector<int> id) {
  if(id.size()!=1){return false;}
  int idtemp(abs(id[0]));
  for(unsigned int ix=0;ix<_id.size();++ix) {
    if(abs(_id[ix])==idtemp) return true;
  }
  return false;
}

unsigned int ScalarMesonCurrent::decayMode(vector<int> idout) {
  int idtemp(abs(idout[0])); unsigned int ix(0);
  bool found(false);
  do {
    if(idtemp==abs(_id[ix])) found=true;
    else                     ++ix;
  }
  while(!found);
  return ix;
}

void ScalarMesonCurrent::dataBaseOutput(ofstream & output,
					bool header,bool create) const {
  if(header) {
    output << "update decayers set parameters=\"";
  }
  if(create) {
    output << "create Herwig::ScalarMesonCurrent " << name() 
	   << " HwWeakCurrents.so\n";
  }
  output << "newdef " << name() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  unsigned int ix;
  for(ix=0;ix<_id.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":ID " << ix 
	     << " " << _id[ix] << "\n";
      output << "newdef " << name() << ":Decay_Constant " << ix 
	     << " " << _decay_constant[ix]/MeV << "\n";
    }
    else {
      output << "insert " << name() << ":ID " << ix 
	     << " " << _id[ix] << "\n";
      output << "insert " << name() << ":Decay_Constant " << ix 
	     << " " << _decay_constant[ix]/MeV << "\n";
    }
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";\n";
  }
}
#line 1 "./ThreePionDefaultCurrent.cc"
// -*- C++ -*-
//
// ThreePionDefaultCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreePionDefaultCurrent class.
//

#include "ThreePionDefaultCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;
using namespace ThePEG;

DescribeClass<ThreePionDefaultCurrent,WeakCurrent>
describeHerwigThreePionDefaultCurrent("Herwig::ThreePionDefaultCurrent",
				       "HwWeakCurrents.so");
HERWIG_INTERPOLATOR_CLASSDESC(ThreePionDefaultCurrent,Energy,Energy2)


ThreePionDefaultCurrent::ThreePionDefaultCurrent() {
  // the quarks for the different modes
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(4);
  // the pion decay constant
  _fpi=130.7*MeV/sqrt(2.);
  _mpi=ZERO;
  // set the initial weights for the resonances
  // the rho weights
  _rhoF123wgts = {1.0,-0.145,0.};
  // local values of the a_1 parameters
  _a1mass  = 1.251*GeV;
  _a1width = 0.599*GeV;
  _a1opt   = true;
  // local values of the rho parameters
  _rhoF123masses = {0.773*GeV,1.370*GeV,1.750*GeV};
  _rhoF123widths = {0.145*GeV,0.510*GeV,0.120*GeV};
  // initialization of the a_1 running width
  _initializea1=false;
  double a1q2in[200]={0,15788.6,31577.3,47365.9,63154.6,78943.2,94731.9,110521,
		       126309,142098,157886,173675,189464,205252,221041,236830,
		       252618,268407,284196,299984,315773,331562,347350,363139,
		       378927,394716,410505,426293,442082,457871,473659,489448,
		       505237,521025,536814,552603,568391,584180,599969,615757,
		       631546,647334,663123,678912,694700,710489,726278,742066,
		       757855,773644,789432,805221,821010,836798,852587,868375,
		       884164,899953,915741,931530,947319,963107,978896,994685,
		       1.01047e+06,1.02626e+06,1.04205e+06,1.05784e+06,1.07363e+06,
		       1.08942e+06,1.10521e+06,1.12099e+06,1.13678e+06,1.15257e+06,
		       1.16836e+06,1.18415e+06,1.19994e+06,1.21573e+06,1.23151e+06,
		       1.2473e+06,1.26309e+06,1.27888e+06,1.29467e+06,1.31046e+06,
		       1.32625e+06,1.34203e+06,1.35782e+06,1.37361e+06,1.3894e+06,
		       1.40519e+06,1.42098e+06,1.43677e+06,1.45256e+06,1.46834e+06
		       ,1.48413e+06,1.49992e+06,1.51571e+06,1.5315e+06,1.54729e+06,
		       1.56308e+06,1.57886e+06,1.59465e+06,1.61044e+06,1.62623e+06,
		       1.64202e+06,1.65781e+06,1.6736e+06,1.68939e+06,1.70517e+06,
		       1.72096e+06,1.73675e+06,1.75254e+06,1.76833e+06,1.78412e+06,
		       1.79991e+06,1.81569e+06,1.83148e+06,1.84727e+06,1.86306e+06,
		       1.87885e+06,1.89464e+06,1.91043e+06,1.92621e+06,1.942e+06,
		       1.95779e+06,1.97358e+06,1.98937e+06,2.00516e+06,2.02095e+06,
		       2.03674e+06,2.05252e+06,2.06831e+06,2.0841e+06,2.09989e+06,
		       2.11568e+06,2.13147e+06,2.14726e+06,2.16304e+06,2.17883e+06,
		       2.19462e+06,2.21041e+06,2.2262e+06,2.24199e+06,2.25778e+06,
		       2.27356e+06,2.28935e+06,2.30514e+06,2.32093e+06,2.33672e+06,
		       2.35251e+06,2.3683e+06,2.38409e+06,2.39987e+06,2.41566e+06,
		       2.43145e+06,2.44724e+06,2.46303e+06,2.47882e+06,2.49461e+06,
		       2.51039e+06,2.52618e+06,2.54197e+06,2.55776e+06,2.57355e+06,
		       2.58934e+06,2.60513e+06,2.62092e+06,2.6367e+06,2.65249e+06,
		       2.66828e+06,2.68407e+06,2.69986e+06,2.71565e+06,2.73144e+06,
		       2.74722e+06,2.76301e+06,2.7788e+06,2.79459e+06,2.81038e+06,
		       2.82617e+06,2.84196e+06,2.85774e+06,2.87353e+06,2.88932e+06,
		       2.90511e+06,2.9209e+06,2.93669e+06,2.95248e+06,2.96827e+06,
		       2.98405e+06,2.99984e+06,3.01563e+06,3.03142e+06,3.04721e+06,
		       3.063e+06,3.07879e+06,3.09457e+06,3.11036e+06,3.12615e+06,
		       3.14194e+06};
  double a1widthin[200]={0,0,0,0,0,0,0,0,0,0,0,0,0.00153933,0.0136382,0.0457614,
			 0.105567,0.199612,0.333825,0.513831,0.745192,1.0336,1.38501,
			 1.80581,2.30295,2.88403,3.5575,4.33278,5.22045,6.23243,
			 7.38223,8.68521,10.1589,11.8234,13.7018,15.8206,18.2107,
			 20.9078,23.9533,27.3954,31.2905,35.7038,40.7106,46.3984,
			 52.8654,60.2207,68.581,78.0637,88.7754,100.794,114.145,
			 128.783,144.574,161.299,178.683,196.426,214.248,231.908,
			 249.221,266.059,282.336,298.006,313.048,327.46,341.254,
			 354.448,367.066,379.133,390.677,401.726,412.304,422.439,
			   432.155,441.474,450.419,459.01,467.267,475.207,482.847,
			 490.203,497.29,504.121,510.71,517.068,523.207,529.138,
			 534.869,540.411,545.776,550.961,556.663,560.851,565.566,
			 570.137,574.569,578.869,583.041,587.091,591.023,594.843,
			 598.553,602.16,605.664,609.072,612.396,615.626,618.754,
			 621.796,624.766,627.656,630.47,633.21,635.878,638.5,
			 641.006,643.471,645.873,648.213,650.493,652.715,654.88,
			 656.99,659.047,661.052,663.007,664.963,666.771,668.6,
			 670.351,672.075,673.828,675.397,676.996,678.567,680.083,
			 681.589,683.023,684.457,685.825,687.18,688.499,689.789,
			 691.058,692.284,693.501,694.667,695.82,696.947,698.05,
			 699.129,700.186,701.221,702.234,703.226,704.198,705.158,
			 706.085,707.001,707.899,708.78,709.644,710.474,711.334,
			 712.145,712.943,713.727,714.505,715.266,716.015,716.751,
			 717.474,718.183,718.88,719.645,720.243,720.91,721.565,
			 722.211,722.851,723.473,724.094,724.697,725.296,725.886,
			 726.468,727.041,727.608,728.166,728.718,729.262,729.808,
			 730.337,730.856,731.374,731.883,732.386,732.884,733.373,
			 733.859,734.339,734.813};

  vector<double> tmp1(a1widthin,a1widthin+200);
  _a1runwidth.clear();
  std::transform(tmp1.begin(), tmp1.end(),
		 back_inserter(_a1runwidth),
		 [](double x){return x*MeV;});
  
  vector<double> tmp2(a1q2in,a1q2in+200);
  _a1runq2.clear();
  std::transform(tmp2.begin(), tmp2.end(),
		 back_inserter(_a1runq2),
		 [](double x){return x*MeV2;});

  _maxmass=ZERO;
  _maxcalc=ZERO;
}

void ThreePionDefaultCurrent::doinit() {
  WeakCurrent::doinit();
  // masses for the running widths
  _mpi = getParticleData(ParticleID::piplus)->mass();
  // initialise the a_1 running width calculation
  inita1Width(-1);
}

void ThreePionDefaultCurrent::persistentOutput(PersistentOStream & os) const {
  os << _rhoF123wgts <<  ounit(_a1runwidth,GeV)<< ounit(_a1runq2,GeV2)
     << _initializea1 << ounit(_a1mass,GeV) << ounit(_a1width,GeV)
     << ounit(_fpi,GeV) << ounit(_mpi,GeV)
     << ounit(_rhoF123masses,GeV) << ounit(_rhoF123widths,GeV)
     << _a1opt << ounit(_maxmass,GeV) << ounit(_maxcalc,GeV) << _a1runinter;
}

void ThreePionDefaultCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _rhoF123wgts >>  iunit(_a1runwidth,GeV) >> iunit(_a1runq2,GeV2) 
     >> _initializea1 >> iunit(_a1mass,GeV) >> iunit(_a1width,GeV)
     >> iunit(_fpi,GeV) >> iunit(_mpi,GeV)
     >> iunit(_rhoF123masses,GeV) >> iunit(_rhoF123widths,GeV)
     >> _a1opt >> iunit(_maxmass,GeV) >> iunit(_maxcalc,GeV) >> _a1runinter;
}

void ThreePionDefaultCurrent::Init() {
        
  static ClassDocumentation<ThreePionDefaultCurrent> documentation
    ("The ThreePionDefaultCurrent class is designed to implement "
     "the three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, " 
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "pi- Kbar0 pi0, pi- pi0 eta. It uses the same currents as those in TAUOLA.",
     "The three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, "
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "and pi- Kbar0 pi0, pi- pi0 eta "
     "use the same currents as \\cite{Jadach:1993hs,Kuhn:1990ad,Decker:1992kj}.",
     "%\\cite{Jadach:1993hs}\n"
     "\\bibitem{Jadach:1993hs}\n"
     "  S.~Jadach, Z.~Was, R.~Decker and J.~H.~Kuhn,\n"
     "  %``The Tau Decay Library Tauola: Version 2.4,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 76}, 361 (1993).\n"
     "  %%CITATION = CPHCB,76,361;%%\n"
     "%\\cite{Kuhn:1990ad}\n"
     "\\bibitem{Kuhn:1990ad}\n"
     "  J.~H.~Kuhn and A.~Santamaria,\n"
     "  %``Tau decays to pions,''\n"
     "  Z.\\ Phys.\\  C {\\bf 48}, 445 (1990).\n"
     "  %%CITATION = ZEPYA,C48,445;%%\n"
     "%\\cite{Decker:1992kj}\n"
     "\\bibitem{Decker:1992kj}\n"
     "  R.~Decker, E.~Mirkes, R.~Sauer and Z.~Was,\n"
     "  %``Tau decays into three pseudoscalar mesons,''\n"
     "  Z.\\ Phys.\\  C {\\bf 58}, 445 (1993).\n"
     "  %%CITATION = ZEPYA,C58,445;%%\n"
     );

  
  static ParVector<ThreePionDefaultCurrent,double> interfaceF123RhoWgt
    ("F123RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &ThreePionDefaultCurrent::_rhoF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static Switch<ThreePionDefaultCurrent,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &ThreePionDefaultCurrent::_initializea1, false, false, false);
  static SwitchOption interfaceInitializea1Initialization
    (interfaceInitializea1,
     "Yes",
     "Initialize the calculation",
     true);
  static SwitchOption interfaceInitializea1NoInitialization
    (interfaceInitializea1,
     "No",
     "Use the default values",
     false);
  
  static Switch<ThreePionDefaultCurrent,bool> interfacea1WidthOption
    ("a1WidthOption",
     "Option for the treatment of the a1 width",
     &ThreePionDefaultCurrent::_a1opt, true, false, false);
  static SwitchOption interfacea1WidthOptionLocal
    (interfacea1WidthOption,
     "Local",
     "Use a calculation of the running width based on the parameters as"
     " interpolation table.",
     true);
  static SwitchOption interfacea1WidthOptionParam
    (interfacea1WidthOption,
     "Kuhn",
     "Use the parameterization of Kuhn and Santamaria for default parameters."
     " This should only be used for testing vs TAUOLA",
     false);

  static ParVector<ThreePionDefaultCurrent,Energy> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &ThreePionDefaultCurrent::_a1runwidth, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<ThreePionDefaultCurrent,Energy2> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &ThreePionDefaultCurrent::_a1runq2, GeV2, -1, 1.0*GeV2, ZERO, 10.0*GeV2,
     false, false, true);
    
  static Parameter<ThreePionDefaultCurrent,Energy> interfaceA1Width
    ("A1Width",
     "The a_1 width if using local values.",
     &ThreePionDefaultCurrent::_a1width, GeV, 0.599*GeV, ZERO, 10.0*GeV,
     false, false, false);
  
  static Parameter<ThreePionDefaultCurrent,Energy> interfaceA1Mass
    ("A1Mass",
     "The a_1 mass if using local values.",
     &ThreePionDefaultCurrent::_a1mass, GeV, 1.251*GeV, ZERO, 10.0*GeV,
     false, false, false);
  
  static ParVector<ThreePionDefaultCurrent,Energy> interfacerhoF123masses
    ("rhoF123masses",
     "The masses for the rho resonances if used local values",
     &ThreePionDefaultCurrent::_rhoF123masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<ThreePionDefaultCurrent,Energy> interfacerhoF123widths
    ("rhoF123widths",
     "The widths for the rho resonances if used local values",
     &ThreePionDefaultCurrent::_rhoF123widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionDefaultCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &ThreePionDefaultCurrent::_fpi, MeV, 92.4*MeV, ZERO, 200.0*MeV,
     false, false, true);
}

// complete the construction of the decay mode for integration
bool ThreePionDefaultCurrent::createMode(int icharge, tcPDPtr resonance,
					  FlavourInfo flavour,
					  unsigned int imode,PhaseSpaceModePtr mode,
					  unsigned int iloc,int ires,
					  PhaseSpaceChannel phase, Energy upp ) {
  // check the charge and resonance
  if(imode==2||imode==3) {
    if(icharge!=0) return false;
    if(resonance && resonance->id()!=ParticleID::a_10) return false;
  }
  else if(imode<2) {
    if(abs(icharge)!=3) return false;
    if(resonance && abs(resonance->id())!=ParticleID::a_1plus) return false;
  }
  else
    assert(false);
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode<=1) return false;
      break;
    case IsoSpin::I3One:
      if( imode>1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if( imode>1 || icharge == 3) return false;
      break;
    default:
      return false;
    }
  }
  // and other flavour
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
  // get the particles and check the masses
  int iq(0),ia(0);
  tPDVector extpart(particles(1,imode,iq,ia));
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  // the particles we will use a lot
  tPDPtr a1;
  if(icharge==-3) {
    a1=getParticleData(ParticleID::a_1minus);
  }
  else if(icharge==3) {
    a1=getParticleData(ParticleID::a_1plus);
  }
  else {
    a1=getParticleData(ParticleID::a_10);
  }
  _maxmass=max(_maxmass,upp);
  // the rho0 resonances
  tPDPtr rho0[3]   = { getParticleData(113), getParticleData(100113), getParticleData(30113)};
  tPDPtr rhoc[3]   = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  if(icharge==3)
    for(unsigned int ix=0;ix<3;++ix) rhoc  [ix] =   rhoc[ix]->CC();
  // create the mode
  if(imode==0) {
    // channels for pi- pi- pi+
    for(unsigned int ix=0;ix<3;++ix) {
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+1,ires+1,rho0[ix],
			ires+2,iloc+2,ires+2,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+2,ires+1,rho0[ix],
			ires+2,iloc+1,ires+2,iloc+3));
    }
  }
  else if(imode==1) {
    // channels for pi0 pi0 pi-
    for(unsigned int ix=0;ix<3;++ix) {
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+1,ires+1,rhoc[ix],
			ires+2,iloc+2,ires+2,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+2,ires+1,rhoc[ix],
			ires+2,iloc+1,ires+2,iloc+3));
    }
  }
  else if(imode>=2) {
    // channels  for pi+ pi- pi0
    for(unsigned int ix=0;ix<3;++ix) {
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+2,ires+1,rhoc[ix]->CC(),
			ires+2,iloc+1,ires+2,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+1,ires+1,rhoc[ix],
			ires+2,iloc+2,ires+2,iloc+3));
    }
  }
  // reset the parameters for the resonances in the integration
  mode->resetIntermediate(a1,_a1mass,_a1width);
  for(unsigned int ix=0;ix<_rhoF123masses.size();++ix) {
    mode->resetIntermediate(rhoc[ix],_rhoF123masses[ix],
			    _rhoF123widths[ix]);
    mode->resetIntermediate(rho0[ix],_rhoF123masses[ix],
			    _rhoF123widths[ix]);
  }
  return true;
}

// initialisation of the a_1 width
// (iopt=-1 initialises, iopt=0 starts the interpolation)
void ThreePionDefaultCurrent::inita1Width(int iopt) {
  if(iopt==-1) {
    _maxcalc=_maxmass;
    if(!_initializea1||_maxmass==ZERO) return;
    // parameters for the table of values
    Energy2 step(sqr(_maxcalc)/199.);
    // integrator to perform the integral
    vector<double> inweights;inweights.push_back(0.5);inweights.push_back(0.5);
    vector<int> intype;intype.push_back(2);intype.push_back(3);
    Energy mrho(getParticleData(ParticleID::rhoplus)->mass()),
      wrho(getParticleData(ParticleID::rhoplus)->width());
    vector<Energy> inmass(2,mrho),inwidth(2,wrho);
    vector<double> inpow(2,0.0);
    ThreeBodyAllOnCalculator<ThreePionDefaultCurrent> 
      widthgen(inweights,intype,inmass,inwidth,inpow,*this,0,_mpi,_mpi,_mpi);
    // normalisation constant to give physical width if on shell
    double a1const(_a1width/(widthgen.partialWidth(sqr(_a1mass))));
    // loop to give the values
    _a1runq2.clear(); _a1runwidth.clear();
    for(Energy2 moff2(ZERO); moff2<=sqr(_maxcalc); moff2+=step) {
      _a1runwidth.push_back(widthgen.partialWidth(moff2)*a1const);
      _a1runq2.push_back(moff2);
    }
  }
  // set up the interpolator
  else if(iopt==0) {
    _a1runinter = make_InterpolatorPtr(_a1runwidth,_a1runq2,3);
  }
}

void ThreePionDefaultCurrent::dataBaseOutput(ofstream & output,bool header,
					      bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::ThreePionDefaultCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<_rhoF123wgts.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F123RhoWeight " << ix << " " << _rhoF123wgts[ix] << "\n";
  }
  output << "newdef " << name() << ":Initializea1 "    << _initializea1    << "\n";
  output << "newdef " << name() << ":a1WidthOption "   << _a1opt           << "\n";
  for(unsigned int ix=0;ix<_a1runwidth.size();++ix) {
    output << "newdef " << name() << ":a1RunningWidth " << ix 
	   << " " << _a1runwidth[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_a1runq2.size();++ix) {
    output << "newdef " << name() << ":a1RunningQ2 " << ix 
	   << " " << _a1runq2[ix]/GeV2 << "\n";
  }
  output << "newdef " << name() << ":A1Width " << _a1width/GeV << "\n";
  output << "newdef " << name() << ":A1Mass "  << _a1mass/GeV  << "\n";
  output << "newdef " << name() << ":FPi "     << _fpi/MeV     << "\n";
  for(unsigned int ix=0;ix<_rhoF123masses.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF123masses " << ix 
	   << " " << _rhoF123masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF123widths.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF123widths " << ix << " " 
	   << _rhoF123widths[ix]/GeV << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

void ThreePionDefaultCurrent::doinitrun() {
  // set up the running a_1 width
  inita1Width(0);
  WeakCurrent::doinitrun();
}

void ThreePionDefaultCurrent::doupdate() {
  WeakCurrent::doupdate();
  // update running width if needed
  if ( !touched() ) return;
  if(_maxmass!=_maxcalc) inita1Width(-1);
}

double ThreePionDefaultCurrent::
threeBodyMatrixElement(const int       , const Energy2 q2,
		       const Energy2 s3, const Energy2 s2, 
		       const Energy2 s1, const Energy    , 
		       const Energy    , const Energy    ) const {
  Energy2 mpi2(sqr(_mpi));
  Complex propb(BrhoF123(s1,-1)),propa(BrhoF123(s2,-1)); 
  // the matrix element
  Energy2 output(ZERO); 
  // first resonance
  output += ((s1-4.*mpi2) + 0.25*(s3-s2)*(s3-s2)/q2) * real(propb*conj(propb)); 
  // second resonance
  output += ((s2-4.*mpi2) + 0.25*(s3-s1)*(s3-s1)/q2) * real(propa*conj(propa)); 
  // the interference term 
  output += (0.5*q2-s3-0.5*mpi2+0.25*(s3-s2)*(s3-s1)/q2)*real(propa*conj(propb)+
							      propb*conj(propa)); 
  return output/sqr(_rhoF123masses[0]);
}


// the hadronic currents    
vector<LorentzPolarizationVectorE> 
ThreePionDefaultCurrent::current(tcPDPtr resonance,
			      FlavourInfo flavour,
			      const int imode, const int ichan, Energy & scale, 
			      const tPDVector & ,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IOne)
    return vector<LorentzPolarizationVectorE>();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode<=1) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode>=2) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode>=2) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return vector<LorentzPolarizationVectorE>();
  useMe();
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q;
  for(unsigned int ix=0;ix<momenta.size();++ix)
    q+=momenta[ix];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 s1 = (momenta[1]+momenta[2]).m2();
  Energy2 s2 = (momenta[0]+momenta[2]).m2();
  Complex F1(0.), F2(0.);
  Complex a1fact(2./3.);
  if(!resonance) a1fact *= a1BreitWigner(q2);
  if(ichan<0) {
    F1 = a1fact*BrhoF123(s1,-1);
    F2 =-a1fact*BrhoF123(s2,-1);
  }
  else if(ichan%2==0) F1 = a1fact*BrhoF123(s1,    ichan/2);
  else if(ichan%2==1) F2 =-a1fact*BrhoF123(s2,(ichan-1)/2);
  // the first three form-factors
  LorentzPolarizationVectorE vect = (F2-F1)*momenta[2] +F1*momenta[1] -F2*momenta[0];
  // multiply by the transverse projection operator
  Complex dot=(vect*q)/q2;
  // scalar and parity violating terms
  vect -= dot*q;
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,q.mass()/_fpi*vect);
}

bool ThreePionDefaultCurrent::accept(vector<int> id) {
  if(id.size()!=3) return false;
  int npip(0),npim(0),npi0(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
  }
  if(      (npip==2&&npim==1) || (npim==2&&npip==1) ) return true;
  else if( (npip==1&&npi0==2) || (npim==1&&npi0==2) ) return true;
  else if( npip==1 && npim==1 && npi0 ==1 )           return true;
  return false;
}

unsigned int ThreePionDefaultCurrent::decayMode(vector<int> id) {
  assert(id.size()==3);
  int npip(0),npim(0),npi0(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
  }
  if(      (npip==2&&npim==1) || (npim==2&&npip==1) ) return 0;
  else if( (npip==1&&npi0==2) || (npim==1&&npi0==2) ) return 1;
  else if( npip==1 && npim==1 && npi0 ==1 )           return 2;
  else assert(false);
}

tPDVector ThreePionDefaultCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart(3);
  if(imode==0) {
    extpart[0]=getParticleData(ParticleID::piminus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::piplus);
  }
  else if(imode==1) {
    extpart[0]=getParticleData(ParticleID::pi0);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::piminus);
  }
  else if(imode==2||imode==3) {
    extpart[0]=getParticleData(ParticleID::piplus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::pi0);
  }
  else
    assert(false);
  // conjugate the particles if needed
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(extpart[ix]->CC()) extpart[ix]=extpart[ix]->CC();
    }
  }
  // return the answer
  return extpart;
}
#line 1 "./OneKaonTwoPionDefaultCurrent.cc"
// -*- C++ -*-
//
// OneKaonTwoPionDefaultCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneKaonTwoPionDefaultCurrent class.
//

#include "OneKaonTwoPionDefaultCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace Herwig;
using namespace ThePEG;

DescribeClass<OneKaonTwoPionDefaultCurrent,WeakCurrent>
describeHerwigOneKaonTwoPionDefaultCurrent("Herwig::OneKaonTwoPionDefaultCurrent",
				       "HwWeakCurrents.so");
HERWIG_INTERPOLATOR_CLASSDESC(OneKaonTwoPionDefaultCurrent,Energy,Energy2)


OneKaonTwoPionDefaultCurrent::OneKaonTwoPionDefaultCurrent() {
  // the quarks for the different modes
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  setInitialModes(3);
  // the pion decay constant
  _fpi=130.7*MeV/sqrt(2.);
  _mpi = ZERO;
  _mK  = ZERO;
  // set the initial weights for the resonances
  // the rho weights
  _rhoF123wgts = {1.0,-0.145,0.};
  // the Kstar weights
  _kstarF123wgts = {1.};
  _kstarF5wgts   = {1.};
  // relative rho/Kstar weights
  _rhoKstarwgt=-0.2;
  // local values of the K_1 parameters
  _k1mass  = 1.402*GeV;
  _k1width = 0.174*GeV;
  // local values of the rho parameters
  _rhoF123masses = {0.773*GeV,1.370*GeV,1.750*GeV};
  _rhoF123widths = {0.145*GeV,0.510*GeV,0.120*GeV};
  // local values for the Kstar parameters
  _kstarF123masses = {0.8921*GeV};
  _kstarF123widths = {0.0513*GeV};
  _kstarF5masses   = {0.8921*GeV};
  _kstarF5widths   = {0.0513*GeV};
}

void OneKaonTwoPionDefaultCurrent::doinit() {
  WeakCurrent::doinit();
  _mpi = getParticleData(ParticleID::piplus)->mass();
  _mK  = getParticleData(ParticleID::Kminus)->mass();

}

void OneKaonTwoPionDefaultCurrent::persistentOutput(PersistentOStream & os) const {
  os << _rhoF123wgts << _kstarF123wgts << _kstarF5wgts
     << _rhoKstarwgt << ounit(_k1mass,GeV)
     << ounit(_k1width,GeV) << ounit(_fpi,GeV) << ounit(_mpi,GeV) << ounit(_mK,GeV)
     << ounit(_rhoF123masses,GeV) << ounit(_rhoF123widths,GeV) 
     << ounit(_kstarF123masses,GeV) << ounit(_kstarF5masses,GeV)
     << ounit(_kstarF123widths,GeV) << ounit(_kstarF5widths,GeV);
}

void OneKaonTwoPionDefaultCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _rhoF123wgts >> _kstarF123wgts >> _kstarF5wgts
     >> _rhoKstarwgt >> iunit(_k1mass,GeV) 
     >> iunit(_k1width,GeV) >> iunit(_fpi,GeV) >> iunit(_mpi,GeV) >> iunit(_mK,GeV)
     >> iunit(_rhoF123masses,GeV) >> iunit(_rhoF123widths,GeV) 
     >> iunit(_kstarF123masses,GeV) >> iunit(_kstarF5masses,GeV)
     >> iunit(_kstarF123widths,GeV) >> iunit(_kstarF5widths,GeV);
}

void OneKaonTwoPionDefaultCurrent::Init() {
        
  static ClassDocumentation<OneKaonTwoPionDefaultCurrent> documentation
    ("The OneKaonTwoPionDefaultCurrent class is designed to implement "
     "the three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, " 
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "pi- Kbar0 pi0, pi- pi0 eta. It uses the same currents as those in TAUOLA.",
     "The three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, "
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "and pi- Kbar0 pi0, pi- pi0 eta "
     "use the same currents as \\cite{Jadach:1993hs,Kuhn:1990ad,Decker:1992kj}.",
     "%\\cite{Jadach:1993hs}\n"
     "\\bibitem{Jadach:1993hs}\n"
     "  S.~Jadach, Z.~Was, R.~Decker and J.~H.~Kuhn,\n"
     "  %``The Tau Decay Library Tauola: Version 2.4,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 76}, 361 (1993).\n"
     "  %%CITATION = CPHCB,76,361;%%\n"
     "%\\cite{Kuhn:1990ad}\n"
     "\\bibitem{Kuhn:1990ad}\n"
     "  J.~H.~Kuhn and A.~Santamaria,\n"
     "  %``Tau decays to pions,''\n"
     "  Z.\\ Phys.\\  C {\\bf 48}, 445 (1990).\n"
     "  %%CITATION = ZEPYA,C48,445;%%\n"
     "%\\cite{Decker:1992kj}\n"
     "\\bibitem{Decker:1992kj}\n"
     "  R.~Decker, E.~Mirkes, R.~Sauer and Z.~Was,\n"
     "  %``Tau decays into three pseudoscalar mesons,''\n"
     "  Z.\\ Phys.\\  C {\\bf 58}, 445 (1993).\n"
     "  %%CITATION = ZEPYA,C58,445;%%\n"
     );

  
  static ParVector<OneKaonTwoPionDefaultCurrent,double> interfaceF123RhoWgt
    ("F123RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &OneKaonTwoPionDefaultCurrent::_rhoF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<OneKaonTwoPionDefaultCurrent,double> interfaceF123KstarWgt
    ("F123KstarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &OneKaonTwoPionDefaultCurrent::_kstarF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<OneKaonTwoPionDefaultCurrent,double> interfaceF5KstarWgt
    ("F5KstarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &OneKaonTwoPionDefaultCurrent::_kstarF5wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static Parameter<OneKaonTwoPionDefaultCurrent,double> interfaceRhoKstarWgt
    ("RhoKstarWgt",
     "The relative weights of the rho and K* in the F5 form factor",
     &OneKaonTwoPionDefaultCurrent::_rhoKstarwgt, -0.2, -10., 10.,
     false, false, false);
  
  static Parameter<OneKaonTwoPionDefaultCurrent,Energy> interfaceK1Width
    ("K1Width",
     "The K_1 width if using local values.",
     &OneKaonTwoPionDefaultCurrent::_k1width, GeV, 0.174*GeV, ZERO, 10.0*GeV,
     false, false, false);
  
  static Parameter<OneKaonTwoPionDefaultCurrent,Energy> interfaceK1Mass
    ("K1Mass",
     "The K_1 mass if using local values.",
     &OneKaonTwoPionDefaultCurrent::_k1mass, GeV, 1.402*GeV, ZERO, 10.0*GeV,
     false, false, false);
  
  static ParVector<OneKaonTwoPionDefaultCurrent,Energy> interfacerhoF123masses
    ("rhoF123masses",
     "The masses for the rho resonances if used local values",
     &OneKaonTwoPionDefaultCurrent::_rhoF123masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionDefaultCurrent,Energy> interfacerhoF123widths
    ("rhoF123widths",
     "The widths for the rho resonances if used local values",
     &OneKaonTwoPionDefaultCurrent::_rhoF123widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<OneKaonTwoPionDefaultCurrent,Energy> interfaceKstarF123masses
    ("KstarF123masses",
     "The masses for the Kstar resonances if used local values",
     &OneKaonTwoPionDefaultCurrent::_kstarF123masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionDefaultCurrent,Energy> interfaceKstarF123widths
    ("KstarF123widths",
     "The widths for the Kstar resonances if used local values",
     &OneKaonTwoPionDefaultCurrent::_kstarF123widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionDefaultCurrent,Energy> interfaceKstarF5masses
    ("KstarF5masses",
     "The masses for the Kstar resonances if used local values",
     &OneKaonTwoPionDefaultCurrent::_kstarF5masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionDefaultCurrent,Energy> interfaceKstarF5widths
    ("KstarF5widths",
     "The widths for the Kstar resonances if used local values",
     &OneKaonTwoPionDefaultCurrent::_kstarF5widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<OneKaonTwoPionDefaultCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &OneKaonTwoPionDefaultCurrent::_fpi, MeV, 92.4*MeV, ZERO, 200.0*MeV,
     false, false, true);
}

// complete the construction of the decay mode for integration
bool OneKaonTwoPionDefaultCurrent::createMode(int icharge, tcPDPtr resonance,
					      FlavourInfo flavour,
					      unsigned int imode,PhaseSpaceModePtr mode,
					      unsigned int iloc,int ires,
					      PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(abs(icharge)!=3) return false; 
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IHalf) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Half:
      if(icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusHalf:
      if(icharge == 3) return false;
      break;
    default:
      return false;
    }
  }
  // strangeness
  if(flavour.strange != Strangeness::Unknown) {
    if(icharge== 3 and flavour.strange != Strangeness::PlusOne ) return false;
    if(icharge==-3 and flavour.strange != Strangeness::MinusOne) return false;
  }
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
  // get external particles and check mass
  int iq(0),ia(0);
  tPDVector extpart(particles(1,imode,iq,ia));
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  // the particles we will use a lot
  tPDPtr k1 = getParticleData(ParticleID::Kstar_1minus);
  if(icharge==3) k1 = k1->CC();
  // the rho0 resonances
  tPDPtr rho0[3]   = { getParticleData(113), getParticleData(100113), getParticleData(30113)};
  tPDPtr rhoc[3]   = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  tPDPtr Kstar0[3] = { getParticleData(313), getParticleData(100313), getParticleData(30313)};
  tPDPtr Kstarc[3] = {getParticleData(-323),getParticleData(-100323),getParticleData(-30323)};
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) {
      rhoc  [ix] =   rhoc[ix]->CC();
      Kstar0[ix] = Kstar0[ix]->CC();
      Kstarc[ix] = Kstarc[ix]->CC();
    }
  }
  if(imode==0) {
    if(resonance && resonance != k1) return false;
    // channels for pi0 pi0 K-
    for(unsigned int ix=0;ix<3;++ix) {
      mode->addChannel((PhaseSpaceChannel(phase),ires,k1,ires+1,iloc+1,ires+1,Kstarc[ix],
			ires+2,iloc+2,ires+2,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,k1,ires+1,iloc+2,ires+1,Kstarc[ix],
			ires+2,iloc+1,ires+2,iloc+3));
    }
  }
  else if(imode==1) {
    // channels for K- pi- pi+
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==k1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1,ires+1,iloc+1,ires+1,rho0[ix],
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1,ires+1,iloc+2,ires+1,Kstar0[ix],
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=Kstarc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,iloc+1,ires+1,rho0[iy],
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,iloc+2,ires+1,Kstar0[iy],
			  ires+2,iloc+1,ires+2,iloc+3));
      }
    }
  }
  else if(imode==2) {
    // channels for pi- kbar0 pi0
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==k1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1,ires+1,iloc+2,ires+1,rhoc[ix],
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=Kstarc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,iloc+1,ires+1,Kstar0[iy],
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,iloc+2,ires+1,rhoc[iy],
			  ires+2,iloc+1,ires+2,iloc+3));
      }
    }
  }
  for(unsigned int ix=0;ix<_rhoF123masses.size();++ix) {
    mode->resetIntermediate(rhoc[ix],_rhoF123masses[ix],_rhoF123widths[ix]);
    mode->resetIntermediate(rho0[ix],_rhoF123masses[ix],_rhoF123widths[ix]);
  }
  for(unsigned int ix=0;ix<_kstarF123masses.size();++ix) {
    mode->resetIntermediate(Kstarc[ix],_kstarF123masses[ix],_kstarF123widths[ix]);
    mode->resetIntermediate(Kstar0[ix],_kstarF123masses[ix],_kstarF123widths[ix]);
  }
  return true;
}

void OneKaonTwoPionDefaultCurrent::dataBaseOutput(ofstream & output,bool header,
					      bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::OneKaonTwoPionDefaultCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<_rhoF123wgts.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F123RhoWeight " << ix << " " << _rhoF123wgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF123wgts.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F123KstarWeight " << ix << " " 
	   << _kstarF123wgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF5wgts.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F5KstarWeight " << ix << " " << _kstarF5wgts[ix] << "\n";
  }
  output << "newdef " << name() << ":RhoKstarWgt "     << _rhoKstarwgt     << "\n";
  output << "newdef " << name() << ":K1Width " << _k1width/GeV << "\n";
  output << "newdef " << name() << ":K1Mass "  << _k1mass/GeV  << "\n";
  output << "newdef " << name() << ":FPi "     << _fpi/MeV     << "\n";
  for(unsigned int ix=0;ix<_rhoF123masses.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF123masses " << ix 
	   << " " << _rhoF123masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF123widths.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF123widths " << ix << " " 
	   << _rhoF123widths[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF123masses.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarF123masses " << ix << " " 
	   << _kstarF123masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF123widths.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarF123widths " << ix << " " 
	   << _kstarF123widths[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF5masses.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarF5masses " << ix << " " 
	   << _kstarF5masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF5widths.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarF5widths " << ix << " " 
	   << _kstarF5widths[ix]/GeV << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
OneKaonTwoPionDefaultCurrent::current(tcPDPtr resonance,
			      FlavourInfo flavour,
			      const int imode, const int ichan, Energy & scale, 
			      const tPDVector & outgoing,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  // charge
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge()+outgoing[2]->iCharge();
  // check the isospin
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IHalf) return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Half:
      if(icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusHalf:
      if(icharge == 3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  // strangeness
  if(flavour.strange != Strangeness::Unknown) {
    if(icharge== 3 and flavour.strange != Strangeness::PlusOne ) return vector<LorentzPolarizationVectorE>();
    if(icharge==-3 and flavour.strange != Strangeness::MinusOne) return vector<LorentzPolarizationVectorE>();
  }
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return vector<LorentzPolarizationVectorE>();
  // check the resonance
  int ires1=-1;
  if(resonance) {
    switch(abs(resonance->id())/1000) {
    case 0:
      ires1=0; break;
    case 100:
      ires1=1; break;
    case  30:
      ires1=2; break;
    case  20:
      ires1=3; break;
    default:
      assert(false);
    }
  }
  useMe();
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q;
  for(unsigned int ix=0;ix<momenta.size();++ix)
    q+=momenta[ix];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 s1 = (momenta[1]+momenta[2]).m2();
  Energy2 s2 = (momenta[0]+momenta[2]).m2();
  // calculatebthe form factors
  Complex F1(0.), F2(0.), F3(0.), F5(0.);
  // calculate the K- pi0 k0
  // calculate the pi0 pi0 K-
  Complex K1fact = ires1<0 || ires1==3 ? Resonance::BreitWignerFW_GN(q2,_k1mass,_k1width) : 0.;
  if(imode==0) {
    K1fact /=6.;
    if(ichan<0) {
      F1 = K1fact*BKstarF123(s1,-1);
      F2 =-K1fact*BKstarF123(s2,-1);
    }
    else if(ichan%2==0) F1 = K1fact*BKstarF123(s1,ichan/2);
    else                F2 =-K1fact*BKstarF123(s2,(ichan-1)/2);
  }
  // calculate the K- pi- pi+
  else if(imode==1) {
    K1fact *= sqrt(2.)/3.;
    if(ichan<0) {
      F1 =-K1fact*  BrhoF123(s1,-1);
      F2 = K1fact*BKstarF123(s2,-1);
      if(ires1<0)
	F5 = -BKstarF123(q2,   -1)*FKrho(s2,s1,-1)*sqrt(2.);
      else if(ires1<3)
	F5 = -BKstarF123(q2,ires1)*FKrho(s2,s1,-1)*sqrt(2.);
      else
	F5 = 0.;
    }
    else if(ichan%8==0) F1 =-K1fact*BrhoF123(s1,ichan/8);
    else if(ichan%8==1) F2 = K1fact*BKstarF123(s2,(ichan-1)/8);
    else                F5 = -BKstarF123(q2,ichan/8)*FKrho(s2,s1,(ichan-2)%8)*sqrt(2.);
  }
  // calculate the pi- K0bar pi0
  else if(imode==2) {
    if(ichan<0) {
      F2 =-K1fact*BrhoF123(s2,-1);
      if(ires1<0)
	F5 =-2.*BKstarF123(q2,   -1)*FKrho(s1,s2,-1);
      else if(ires1<3)
	F5 =-2.*BKstarF123(q2,ires1)*FKrho(s1,s2,-1);
      else
	F5  =0.;
    }
    else if(ichan%7==0) F2 =-K1fact*BrhoF123(s2,ichan/7);
    else                F5 =-2.*BKstarF123(q2,ichan/7)*FKrho(s1,s2,(ichan-1)%7);
  }
  // the first three form-factors
  LorentzPolarizationVectorE vect;
  vect = (F2-F1)*momenta[2]
        +(F1-F3)*momenta[1]
        +(F3-F2)*momenta[0];
  // multiply by the transverse projection operator
  Complex dot=(vect*q)/q2;
  // scalar and parity violating terms
  vect -= dot*q;
  if(F5!=0.) {
    using Constants::twopi;
    vect -= Complex(0.,1.)*F5/sqr(twopi)/sqr(_fpi)*
      Helicity::epsilon(momenta[0],momenta[1],momenta[2]);
  }
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,q.mass()/_fpi*vect);
}

bool OneKaonTwoPionDefaultCurrent::accept(vector<int> id) {
  if(id.size()!=3) return false;
  int npip(0),npim(0),nkp(0),nkm(0);
  int npi0(0),nk0(0),nk0bar(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
  }
  if( (nkp==1&&npi0==2) || (npi0==2&&nkm==1) )   return 0;
  else if( (npip==1&&npim==1&&nkp==1) ||
	   (nkm==1&&npim==1&&npip==1) )          return 1;
  else if( (nk0==1&&npip==1&&npi0==1)  ||
	   (npim==1&&nk0bar==1&&npi0==1))        return 2;
  return -1;
}

unsigned int OneKaonTwoPionDefaultCurrent::decayMode(vector<int> id) {
  int npip(0),npim(0),nkp(0),nkm(0);
  int  npi0(0),nk0(0),nk0bar(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
  }
  if( (nkp==1&&npi0==2) || (npi0==2&&nkm==1) )   return 0;
  else if( (npip==1&&npim==1&&nkp==1) ||
	   (nkm==1&&npim==1&&npip==1) )          return 1;
  else if( (nk0==1&&npip==1&&npi0==1)  ||
	   (npim==1&&nk0bar==1&&npi0==1))        return 2;
  assert(false);
}

tPDVector OneKaonTwoPionDefaultCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart(3);
  if(imode==0) {
    extpart[0]=getParticleData(ParticleID::pi0);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::Kminus);
  }
  else if(imode==1) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::piplus);
  }
  else if(imode==2) {
    extpart[0]=getParticleData(ParticleID::piminus);
    extpart[1]=getParticleData(ParticleID::Kbar0);
    extpart[2]=getParticleData(ParticleID::pi0);
  }
  // conjugate the particles if needed
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(extpart[ix]->CC()) extpart[ix]=extpart[ix]->CC();
    }
  }
  // return the answer
  return extpart;
}
#line 1 "./TwoKaonOnePionDefaultCurrent.cc"
// -*- C++ -*-
//
// TwoKaonOnePionDefaultCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoKaonOnePionDefaultCurrent class.
//

#include "TwoKaonOnePionDefaultCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace Herwig;
using namespace ThePEG;

DescribeClass<TwoKaonOnePionDefaultCurrent,WeakCurrent>
describeHerwigTwoKaonOnePionDefaultCurrent("Herwig::TwoKaonOnePionDefaultCurrent",
					   "HwWeakCurrents.so");
HERWIG_INTERPOLATOR_CLASSDESC(TwoKaonOnePionDefaultCurrent,Energy,Energy2)


TwoKaonOnePionDefaultCurrent::TwoKaonOnePionDefaultCurrent() {
  // the quarks for the different modes
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  setInitialModes(3);
  // the pion decay constant
  _fpi=130.7*MeV/sqrt(2.);
  _mpi = ZERO;
  _mK  = ZERO;
  // set the initial weights for the resonances
  // the rho weights
  _rhoF123wgts = { 1.0,-0.145,0.};
  _rhoF5wgts   = {-26.,   6.5,1.};
  // the Kstar weights
  _kstarF123wgts = {1.};
  // relative rho/Kstar weights
  _rhoKstarwgt=-0.2;
  // local values of the a_1 parameters
  _a1mass  = 1.251*GeV;
  _a1width = 0.599*GeV;
  _a1opt=true;
  // local values of the rho parameters
  _rhoF123masses = {0.773*GeV,1.370*GeV,1.750*GeV};
  _rhoF123widths = {0.145*GeV,0.510*GeV,0.120*GeV};
  _rhoF5masses   = {0.773*GeV,1.500*GeV,1.750*GeV};
  _rhoF5widths   = {0.145*GeV,0.220*GeV,0.120*GeV};
  // local values for the Kstar parameters
  _kstarF123masses = {0.8921*GeV};
  _kstarF123widths = {0.0513*GeV};
  // initialization of the a_1 running width
  _initializea1=false;
  double a1q2in[200]={0,15788.6,31577.3,47365.9,63154.6,78943.2,94731.9,110521,
		       126309,142098,157886,173675,189464,205252,221041,236830,
		       252618,268407,284196,299984,315773,331562,347350,363139,
		       378927,394716,410505,426293,442082,457871,473659,489448,
		       505237,521025,536814,552603,568391,584180,599969,615757,
		       631546,647334,663123,678912,694700,710489,726278,742066,
		       757855,773644,789432,805221,821010,836798,852587,868375,
		       884164,899953,915741,931530,947319,963107,978896,994685,
		       1.01047e+06,1.02626e+06,1.04205e+06,1.05784e+06,1.07363e+06,
		       1.08942e+06,1.10521e+06,1.12099e+06,1.13678e+06,1.15257e+06,
		       1.16836e+06,1.18415e+06,1.19994e+06,1.21573e+06,1.23151e+06,
		       1.2473e+06,1.26309e+06,1.27888e+06,1.29467e+06,1.31046e+06,
		       1.32625e+06,1.34203e+06,1.35782e+06,1.37361e+06,1.3894e+06,
		       1.40519e+06,1.42098e+06,1.43677e+06,1.45256e+06,1.46834e+06
		       ,1.48413e+06,1.49992e+06,1.51571e+06,1.5315e+06,1.54729e+06,
		       1.56308e+06,1.57886e+06,1.59465e+06,1.61044e+06,1.62623e+06,
		       1.64202e+06,1.65781e+06,1.6736e+06,1.68939e+06,1.70517e+06,
		       1.72096e+06,1.73675e+06,1.75254e+06,1.76833e+06,1.78412e+06,
		       1.79991e+06,1.81569e+06,1.83148e+06,1.84727e+06,1.86306e+06,
		       1.87885e+06,1.89464e+06,1.91043e+06,1.92621e+06,1.942e+06,
		       1.95779e+06,1.97358e+06,1.98937e+06,2.00516e+06,2.02095e+06,
		       2.03674e+06,2.05252e+06,2.06831e+06,2.0841e+06,2.09989e+06,
		       2.11568e+06,2.13147e+06,2.14726e+06,2.16304e+06,2.17883e+06,
		       2.19462e+06,2.21041e+06,2.2262e+06,2.24199e+06,2.25778e+06,
		       2.27356e+06,2.28935e+06,2.30514e+06,2.32093e+06,2.33672e+06,
		       2.35251e+06,2.3683e+06,2.38409e+06,2.39987e+06,2.41566e+06,
		       2.43145e+06,2.44724e+06,2.46303e+06,2.47882e+06,2.49461e+06,
		       2.51039e+06,2.52618e+06,2.54197e+06,2.55776e+06,2.57355e+06,
		       2.58934e+06,2.60513e+06,2.62092e+06,2.6367e+06,2.65249e+06,
		       2.66828e+06,2.68407e+06,2.69986e+06,2.71565e+06,2.73144e+06,
		       2.74722e+06,2.76301e+06,2.7788e+06,2.79459e+06,2.81038e+06,
		       2.82617e+06,2.84196e+06,2.85774e+06,2.87353e+06,2.88932e+06,
		       2.90511e+06,2.9209e+06,2.93669e+06,2.95248e+06,2.96827e+06,
		       2.98405e+06,2.99984e+06,3.01563e+06,3.03142e+06,3.04721e+06,
		       3.063e+06,3.07879e+06,3.09457e+06,3.11036e+06,3.12615e+06,
		       3.14194e+06};
  double a1widthin[200]={0,0,0,0,0,0,0,0,0,0,0,0,0.00153933,0.0136382,0.0457614,
			 0.105567,0.199612,0.333825,0.513831,0.745192,1.0336,1.38501,
			 1.80581,2.30295,2.88403,3.5575,4.33278,5.22045,6.23243,
			 7.38223,8.68521,10.1589,11.8234,13.7018,15.8206,18.2107,
			 20.9078,23.9533,27.3954,31.2905,35.7038,40.7106,46.3984,
			 52.8654,60.2207,68.581,78.0637,88.7754,100.794,114.145,
			 128.783,144.574,161.299,178.683,196.426,214.248,231.908,
			 249.221,266.059,282.336,298.006,313.048,327.46,341.254,
			 354.448,367.066,379.133,390.677,401.726,412.304,422.439,
			   432.155,441.474,450.419,459.01,467.267,475.207,482.847,
			 490.203,497.29,504.121,510.71,517.068,523.207,529.138,
			 534.869,540.411,545.776,550.961,556.663,560.851,565.566,
			 570.137,574.569,578.869,583.041,587.091,591.023,594.843,
			 598.553,602.16,605.664,609.072,612.396,615.626,618.754,
			 621.796,624.766,627.656,630.47,633.21,635.878,638.5,
			 641.006,643.471,645.873,648.213,650.493,652.715,654.88,
			 656.99,659.047,661.052,663.007,664.963,666.771,668.6,
			 670.351,672.075,673.828,675.397,676.996,678.567,680.083,
			 681.589,683.023,684.457,685.825,687.18,688.499,689.789,
			 691.058,692.284,693.501,694.667,695.82,696.947,698.05,
			 699.129,700.186,701.221,702.234,703.226,704.198,705.158,
			 706.085,707.001,707.899,708.78,709.644,710.474,711.334,
			 712.145,712.943,713.727,714.505,715.266,716.015,716.751,
			 717.474,718.183,718.88,719.645,720.243,720.91,721.565,
			 722.211,722.851,723.473,724.094,724.697,725.296,725.886,
			 726.468,727.041,727.608,728.166,728.718,729.262,729.808,
			 730.337,730.856,731.374,731.883,732.386,732.884,733.373,
			 733.859,734.339,734.813};

  vector<double> tmp1(a1widthin,a1widthin+200);
  _a1runwidth.clear();
  std::transform(tmp1.begin(), tmp1.end(),
		 back_inserter(_a1runwidth),
		 [](double x){return x*MeV;});
  
  vector<double> tmp2(a1q2in,a1q2in+200);
  _a1runq2.clear();
  std::transform(tmp2.begin(), tmp2.end(),
		 back_inserter(_a1runq2),
		 [](double x){return x*MeV2;});
  _maxmass=ZERO;
  _maxcalc=ZERO;
}

void TwoKaonOnePionDefaultCurrent::doinit() {
  WeakCurrent::doinit();
  // masses for the running widths
  _mpi = getParticleData(ParticleID::piplus)->mass();
  _mK  = getParticleData(ParticleID::Kminus)->mass();
  // initialise the a_1 running width calculation
  inita1Width(-1);
}

void TwoKaonOnePionDefaultCurrent::persistentOutput(PersistentOStream & os) const {
  os << _rhoF123wgts << _kstarF123wgts << _rhoF5wgts
     << _rhoKstarwgt <<  ounit(_a1runwidth,GeV)<< ounit(_a1runq2,GeV2)
     << _initializea1 << ounit(_a1mass,GeV)<< ounit(_a1width,GeV)
     << ounit(_fpi,GeV) << ounit(_mpi,GeV)<< ounit(_mK,GeV)
     << ounit(_rhoF123masses,GeV) << ounit(_rhoF5masses,GeV) 
     << ounit(_rhoF123widths,GeV) 
     << ounit(_rhoF5widths,GeV) << ounit(_kstarF123masses,GeV) 
     << ounit(_kstarF123widths,GeV)
     << _a1opt << ounit(_maxmass,GeV) << ounit(_maxcalc,GeV) << _a1runinter;
}

void TwoKaonOnePionDefaultCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _rhoF123wgts >> _kstarF123wgts >> _rhoF5wgts 
     >> _rhoKstarwgt >>  iunit(_a1runwidth,GeV) >> iunit(_a1runq2,GeV2) 
     >> _initializea1 >> iunit(_a1mass,GeV) >> iunit(_a1width,GeV)
     >> iunit(_fpi,GeV) >> iunit(_mpi,GeV) >> iunit(_mK,GeV)
     >> iunit(_rhoF123masses,GeV) >> iunit(_rhoF5masses,GeV) 
     >> iunit(_rhoF123widths,GeV) 
     >> iunit(_rhoF5widths,GeV) >> iunit(_kstarF123masses,GeV) 
     >> iunit(_kstarF123widths,GeV)
     >> _a1opt >> iunit(_maxmass,GeV) >> iunit(_maxcalc,GeV) >> _a1runinter;
}

void TwoKaonOnePionDefaultCurrent::Init() {
        
  static ClassDocumentation<TwoKaonOnePionDefaultCurrent> documentation
    ("The TwoKaonOnePionDefaultCurrent class is designed to implement "
     "the three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, " 
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "pi- Kbar0 pi0, pi- pi0 eta. It uses the same currents as those in TAUOLA.",
     "The three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, "
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "and pi- Kbar0 pi0, pi- pi0 eta "
     "use the same currents as \\cite{Jadach:1993hs,Kuhn:1990ad,Decker:1992kj}.",
     "%\\cite{Jadach:1993hs}\n"
     "\\bibitem{Jadach:1993hs}\n"
     "  S.~Jadach, Z.~Was, R.~Decker and J.~H.~Kuhn,\n"
     "  %``The Tau Decay Library Tauola: Version 2.4,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 76}, 361 (1993).\n"
     "  %%CITATION = CPHCB,76,361;%%\n"
     "%\\cite{Kuhn:1990ad}\n"
     "\\bibitem{Kuhn:1990ad}\n"
     "  J.~H.~Kuhn and A.~Santamaria,\n"
     "  %``Tau decays to pions,''\n"
     "  Z.\\ Phys.\\  C {\\bf 48}, 445 (1990).\n"
     "  %%CITATION = ZEPYA,C48,445;%%\n"
     "%\\cite{Decker:1992kj}\n"
     "\\bibitem{Decker:1992kj}\n"
     "  R.~Decker, E.~Mirkes, R.~Sauer and Z.~Was,\n"
     "  %``Tau decays into three pseudoscalar mesons,''\n"
     "  Z.\\ Phys.\\  C {\\bf 58}, 445 (1993).\n"
     "  %%CITATION = ZEPYA,C58,445;%%\n"
     );

  static ParVector<TwoKaonOnePionDefaultCurrent,double> interfaceF123RhoWgt
    ("F123RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &TwoKaonOnePionDefaultCurrent::_rhoF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<TwoKaonOnePionDefaultCurrent,double> interfaceF123KstarWgt
    ("F123KstarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &TwoKaonOnePionDefaultCurrent::_kstarF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<TwoKaonOnePionDefaultCurrent,double> interfaceF5RhoWgt
    ("F5RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &TwoKaonOnePionDefaultCurrent::_rhoF5wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static Parameter<TwoKaonOnePionDefaultCurrent,double> interfaceRhoKstarWgt
    ("RhoKstarWgt",
     "The relative weights of the rho and K* in the F5 form factor",
     &TwoKaonOnePionDefaultCurrent::_rhoKstarwgt, -0.2, -10., 10.,
     false, false, false);
  
  static Switch<TwoKaonOnePionDefaultCurrent,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &TwoKaonOnePionDefaultCurrent::_initializea1, false, false, false);
  static SwitchOption interfaceInitializea1Initialization
    (interfaceInitializea1,
     "Yes",
     "Initialize the calculation",
     true);
  static SwitchOption interfaceInitializea1NoInitialization
    (interfaceInitializea1,
     "No",
     "Use the default values",
     false);
  
  static Switch<TwoKaonOnePionDefaultCurrent,bool> interfacea1WidthOption
    ("a1WidthOption",
     "Option for the treatment of the a1 width",
     &TwoKaonOnePionDefaultCurrent::_a1opt, true, false, false);
  static SwitchOption interfacea1WidthOptionLocal
    (interfacea1WidthOption,
     "Local",
     "Use a calculation of the running width based on the parameters as"
     " interpolation table.",
     true);
  static SwitchOption interfacea1WidthOptionParam
    (interfacea1WidthOption,
     "Kuhn",
     "Use the parameterization of Kuhn and Santamaria for default parameters."
     " This should only be used for testing vs TAUOLA",
     false);

  static ParVector<TwoKaonOnePionDefaultCurrent,Energy> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &TwoKaonOnePionDefaultCurrent::_a1runwidth, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<TwoKaonOnePionDefaultCurrent,Energy2> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &TwoKaonOnePionDefaultCurrent::_a1runq2, GeV2, -1, 1.0*GeV2, ZERO, 10.0*GeV2,
     false, false, true);
    
  static Parameter<TwoKaonOnePionDefaultCurrent,Energy> interfaceA1Width
    ("A1Width",
     "The a_1 width if using local values.",
     &TwoKaonOnePionDefaultCurrent::_a1width, GeV, 0.599*GeV, ZERO, 10.0*GeV,
     false, false, false);
  
  static Parameter<TwoKaonOnePionDefaultCurrent,Energy> interfaceA1Mass
    ("A1Mass",
     "The a_1 mass if using local values.",
     &TwoKaonOnePionDefaultCurrent::_a1mass, GeV, 1.251*GeV, ZERO, 10.0*GeV,
     false, false, false);
  
  static ParVector<TwoKaonOnePionDefaultCurrent,Energy> interfacerhoF123masses
    ("rhoF123masses",
     "The masses for the rho resonances if used local values",
     &TwoKaonOnePionDefaultCurrent::_rhoF123masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionDefaultCurrent,Energy> interfacerhoF123widths
    ("rhoF123widths",
     "The widths for the rho resonances if used local values",
     &TwoKaonOnePionDefaultCurrent::_rhoF123widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionDefaultCurrent,Energy> interfacerhoF5masses
    ("rhoF5masses",
     "The masses for the rho resonances if used local values",
     &TwoKaonOnePionDefaultCurrent::_rhoF5masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionDefaultCurrent,Energy> interfacerhoF5widths
    ("rhoF5widths",
     "The widths for the rho resonances if used local values",
     &TwoKaonOnePionDefaultCurrent::_rhoF5widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<TwoKaonOnePionDefaultCurrent,Energy> interfaceKstarF123masses
    ("KstarF123masses",
     "The masses for the Kstar resonances if used local values",
     &TwoKaonOnePionDefaultCurrent::_kstarF123masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionDefaultCurrent,Energy> interfaceKstarF123widths
    ("KstarF123widths",
     "The widths for the Kstar resonances if used local values",
     &TwoKaonOnePionDefaultCurrent::_kstarF123widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<TwoKaonOnePionDefaultCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &TwoKaonOnePionDefaultCurrent::_fpi, MeV, 92.4*MeV, ZERO, 200.0*MeV,
     false, false, true);
}

// complete the construction of the decay mode for integration
bool TwoKaonOnePionDefaultCurrent::createMode(int icharge, tcPDPtr resonance,
					      FlavourInfo flavour,
					      unsigned int imode,PhaseSpaceModePtr mode,
					      unsigned int iloc,int ires,
					      PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(abs(icharge)!=3) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode<=1) return false;
      break;
    case IsoSpin::I3One:
      if( imode>1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if( imode>1 || icharge == 3) return false;
      break;
    default:
      return false;
    }
  }
  // other flavour stuff
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  != Beauty::Zero     ) return false;
  // get the particles and check the mass
  int iq(0),ia(0);
  tPDVector extpart(particles(1,imode,iq,ia));
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  // the particles we will use a lot
  tPDPtr a1=getParticleData(ParticleID::a_1minus);
  if(icharge==3) a1=a1->CC();
  _maxmass=max(_maxmass,upp);
  // the rho0 resonances
  tPDPtr rho0[3]   = { getParticleData(113), getParticleData(100113), getParticleData(30113)};
  tPDPtr rhoc[3]   = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  tPDPtr Kstar0[3] = { getParticleData(313), getParticleData(100313), getParticleData(30313)};
  tPDPtr Kstarc[3] = {getParticleData(-323),getParticleData(-100323),getParticleData(-30323)};
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) {
      rhoc  [ix] =   rhoc[ix]->CC();
      Kstar0[ix] = Kstar0[ix]->CC();
      Kstarc[ix] = Kstarc[ix]->CC();
    }
  }
  if(imode==0) {
    // channels for K- pi- K+
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==a1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+1,ires+1,Kstar0[ix],
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+2,ires+1,rho0[ix],
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=rhoc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,iloc+1,ires+1,Kstar0[iy],
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,iloc+2,ires+1,rho0[iy],
			  ires+2,iloc+1,ires+2,iloc+3));
      }
    }
  }
  else if(imode==1) {
    // channels for K0 pi- K0bar
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==a1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+1,ires+1,Kstarc[ix],
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+2,ires+1,rho0[ix],
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=rhoc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,iloc+1,ires+1,Kstarc[iy],
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,iloc+2,ires+1,rho0[iy],
			  ires+2,iloc+1,ires+2,iloc+3));
      }
    }
  }
  else if(imode==2) {
    // channels for K- pi0 K0
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==a1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+2,ires+1,rhoc[ix],
			  ires+2,iloc+1,ires+2,iloc+3));
      }
    }
  }
  for(unsigned int ix=0;ix<_rhoF123masses.size();++ix) {
    mode->resetIntermediate(rhoc[ix],_rhoF123masses[ix],
			    _rhoF123widths[ix]);
    mode->resetIntermediate(rho0[ix],_rhoF123masses[ix],
			    _rhoF123widths[ix]);
  }
  // K star parameters in the base class
  for(unsigned int ix=0;ix<_kstarF123masses.size();++ix) {
    mode->resetIntermediate(Kstarc[ix],_kstarF123masses[ix],
			    _kstarF123widths[ix]);
    mode->resetIntermediate(Kstar0[ix],_kstarF123masses[ix],
			    _kstarF123widths[ix]);
  }
  return true;
}

// initialisation of the a_1 width
// (iopt=-1 initialises, iopt=0 starts the interpolation)
void TwoKaonOnePionDefaultCurrent::inita1Width(int iopt) {
  if(iopt==-1) {
    _maxcalc=_maxmass;
    if(!_initializea1||_maxmass==ZERO) return;
    // parameters for the table of values
    Energy2 step(sqr(_maxcalc)/199.);
    // integrator to perform the integral
    vector<double> inweights;inweights.push_back(0.5);inweights.push_back(0.5);
    vector<int> intype;intype.push_back(2);intype.push_back(3);
    Energy mrho(getParticleData(ParticleID::rhoplus)->mass()),
      wrho(getParticleData(ParticleID::rhoplus)->width());
    vector<Energy> inmass(2,mrho),inwidth(2,wrho);
    vector<double> inpow(2,0.0);
    ThreeBodyAllOnCalculator<TwoKaonOnePionDefaultCurrent> 
      widthgen(inweights,intype,inmass,inwidth,inpow,*this,0,_mpi,_mpi,_mpi);
    // normalisation constant to give physical width if on shell
    double a1const(_a1width/(widthgen.partialWidth(sqr(_a1mass))));
    // loop to give the values
    _a1runq2.clear(); _a1runwidth.clear();
    for(Energy2 moff2(ZERO); moff2<=sqr(_maxcalc); moff2+=step) {
      _a1runwidth.push_back(widthgen.partialWidth(moff2)*a1const);
      _a1runq2.push_back(moff2);
    }
  }
  // set up the interpolator
  else if(iopt==0) {
    _a1runinter = make_InterpolatorPtr(_a1runwidth,_a1runq2,3);
  }
}

void TwoKaonOnePionDefaultCurrent::dataBaseOutput(ofstream & output,bool header,
					      bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::TwoKaonOnePionDefaultCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<_rhoF123wgts.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F123RhoWeight " << ix << " " << _rhoF123wgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF123wgts.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F123KstarWeight " << ix << " " 
	   << _kstarF123wgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF5wgts.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F5RhoWeight " << ix << " " << _rhoF5wgts[ix] << "\n";
  }
  output << "newdef " << name() << ":RhoKstarWgt "     << _rhoKstarwgt     << "\n";
  output << "newdef " << name() << ":Initializea1 "    << _initializea1    << "\n";
  output << "newdef " << name() << ":a1WidthOption "   << _a1opt           << "\n";
  for(unsigned int ix=0;ix<_a1runwidth.size();++ix) {
    output << "newdef " << name() << ":a1RunningWidth " << ix 
	   << " " << _a1runwidth[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_a1runq2.size();++ix) {
    output << "newdef " << name() << ":a1RunningQ2 " << ix 
	   << " " << _a1runq2[ix]/GeV2 << "\n";
  }
  output << "newdef " << name() << ":A1Width " << _a1width/GeV << "\n";
  output << "newdef " << name() << ":A1Mass "  << _a1mass/GeV  << "\n";
  output << "newdef " << name() << ":FPi "     << _fpi/MeV     << "\n";
  for(unsigned int ix=0;ix<_rhoF123masses.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF123masses " << ix 
	   << " " << _rhoF123masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF123widths.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF123widths " << ix << " " 
	   << _rhoF123widths[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF5masses.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF5masses " << ix << " " 
	   << _rhoF5masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF5widths.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF5widths " << ix << " " 
	   << _rhoF5widths[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF123masses.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarF123masses " << ix << " " 
	   << _kstarF123masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF123widths.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarF123widths " << ix << " " 
	   << _kstarF123widths[ix]/GeV << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

void TwoKaonOnePionDefaultCurrent::doinitrun() {
  // set up the running a_1 width
  inita1Width(0);
  WeakCurrent::doinitrun();
}

void TwoKaonOnePionDefaultCurrent::doupdate() {
  WeakCurrent::doupdate();
  // update running width if needed
  if ( !touched() ) return;
  if(_maxmass!=_maxcalc) inita1Width(-1);
}

double TwoKaonOnePionDefaultCurrent::
threeBodyMatrixElement(const int       , const Energy2 q2,
		       const Energy2 s3, const Energy2 s2, 
		       const Energy2 s1, const Energy    , 
		       const Energy    , const Energy    ) const {
  Energy2 mpi2(sqr(_mpi));
  Complex propb(BrhoF123(s1,-1)),propa(BrhoF123(s2,-1)); 
  // the matrix element
  Energy2 output(ZERO); 
  // first resonance
  output += ((s1-4.*mpi2) + 0.25*(s3-s2)*(s3-s2)/q2) * real(propb*conj(propb)); 
  // second resonance
  output += ((s2-4.*mpi2) + 0.25*(s3-s1)*(s3-s1)/q2) * real(propa*conj(propa)); 
  // the interference term 
  output += (0.5*q2-s3-0.5*mpi2+0.25*(s3-s2)*(s3-s1)/q2)*real(propa*conj(propb)+
							      propb*conj(propa)); 
  return output/sqr(_rhoF123masses[0]);
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
TwoKaonOnePionDefaultCurrent::current(tcPDPtr resonance,
				      FlavourInfo flavour,
				      const int imode, const int ichan, Energy & scale, 
				      const tPDVector & ,
				      const vector<Lorentz5Momentum> & momenta,
				      DecayIntegrator::MEOption) const {
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IOne)
    return vector<LorentzPolarizationVectorE>();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One: case IsoSpin::I3MinusOne:
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  // other flavour stuff
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  != Beauty::Zero     ) return vector<LorentzPolarizationVectorE>();
  // check the resonance
  int ires1=-1;
  if(resonance) {
    switch(abs(resonance->id())/1000) {
    case 0:
      ires1=0; break;
    case 100:
      ires1=1; break;
    case  30:
      ires1=2; break;
    case  10:
      ires1=3; break;
    default:
      assert(false);
    }
  }
  useMe();
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q;
  for(unsigned int ix=0;ix<momenta.size();++ix)
    q+=momenta[ix];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 s1 = (momenta[1]+momenta[2]).m2();
  Energy2 s2 = (momenta[0]+momenta[2]).m2();
  Complex F1(0.), F2(0.), F5(0.);
  Complex a1fact = ires1<0 || ires1==3 ? a1BreitWigner(q2) : 0.;
  // calculate the K- pi - K+ factor
  if(imode==0) {
    a1fact *= sqrt(2.)/3.;
    if(ichan<0) {
      F1 =-a1fact*BKstarF123(s1,-1); 
      F2 = a1fact*BrhoF123(s2,-1);
      if(ires1<0)
	F5 = BrhoF5(q2,-1)*FKrho(s1,s2,-1)*sqrt(2.);
      else if(ires1<3)
	F5 = BrhoF5(q2,ires1)*FKrho(s1,s2,-1)*sqrt(2.);
      else
	F5 = 0.;
    }
    else if(ichan%8==0) F1 =-a1fact*BKstarF123(s1,ichan/8);
    else if(ichan%8==1) F2 = a1fact*BrhoF123(s2,(ichan-1)/8);
    else if(ichan%8>=2) F5 = BrhoF5(q2,ichan/8)*FKrho(s1,s2,(ichan-2)%8)*sqrt(2.);
  }
  // calculate the K0 pi- K0bar
  else if(imode==1) {
    a1fact *= sqrt(2.)/3.;
    if(ichan<0) {
      F1 =-a1fact*BKstarF123(s1,-1);
      F2 = a1fact*BrhoF123(s2,-1);
      if(ires1<0)
	F5 =-BrhoF5(q2,-1)*FKrho(s1,s2,-1)*sqrt(2.);
      else if(ires1<3)
	F5 =-BrhoF5(q2,ires1)*FKrho(s1,s2,-1)*sqrt(2.);
      else
	F5 = 0.;
    }
    else if(ichan%8==0) F1 = -a1fact*BKstarF123(s1,ichan/8);
    else if(ichan%8==1) F2 = a1fact*BrhoF123(s2,(ichan-1)/8);
    else if(ichan%8>=2) F5 = -BrhoF5(q2,ichan/8)*FKrho(s1,s2,(ichan-2)%8)*sqrt(2.);
  }
  // calculate the K- pi0 k0
  else if(imode==2) {
    if(ichan<0) F2 =-a1fact*BrhoF123(s2,-1);
    else        F2 =-a1fact*BrhoF123(s2,ichan);
  }
  // the first three form-factors
  LorentzPolarizationVectorE vect =
    (F2-F1)*momenta[2] + F1*momenta[1] - F2*momenta[0];
  // multiply by the transverse projection operator
  Complex dot=(vect*q)/q2;
  // scalar and parity violating terms
  vect -=  dot*q;
  if(F5!=0.) 
    vect -= Complex(0.,1.)*F5/sqr(Constants::twopi)/sqr(_fpi)*
      Helicity::epsilon(momenta[0],momenta[1],momenta[2]);
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,q.mass()/_fpi*vect);
}

bool TwoKaonOnePionDefaultCurrent::accept(vector<int> id) {
  int npip(0),npim(0),nkp(0),nkm(0);
  int npi0(0),nk0(0),nk0bar(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
  }
  if      ( (nkp==1&&nkm==1&&npip==1) ||
	   (nkp==1&&nkm==1&&npim==1))       return true;
  else if( (nk0==1&&nk0bar==1&&npip==1) ||
	   (nk0==1&&nk0bar==1&&npim==1))    return true;
  else if( (nkp==1&&nk0bar==1&&npi0==1) ||
	   (nkm==1&&npi0==1&&nk0==1))       return true;
  return false;
}

unsigned int TwoKaonOnePionDefaultCurrent::decayMode(vector<int> id) {
  int npip(0),npim(0),nkp(0),nkm(0);
  int npi0(0),nk0(0),nk0bar(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
  }
  if      ( (nkp==1&&nkm==1&&npip==1) ||
	   (nkp==1&&nkm==1&&npim==1))       return 0;
  else if( (nk0==1&&nk0bar==1&&npip==1) ||
	   (nk0==1&&nk0bar==1&&npim==1))    return 1;
  else if( (nkp==1&&nk0bar==1&&npi0==1) ||
	   (nkm==1&&npi0==1&&nk0==1))       return 2;
  assert(false);
}

tPDVector TwoKaonOnePionDefaultCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart(3);
  if(imode==0) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::Kplus);
  }
  else if(imode==1) {
    extpart[0]=getParticleData(ParticleID::K0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::Kbar0);
  }
  else if(imode==2) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::K0);
  }
  // conjugate the particles if needed
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(extpart[ix]->CC()) extpart[ix]=extpart[ix]->CC();
    }
  }
  // return the answer
  return extpart;
}
#line 1 "./EtaPiPiDefaultCurrent.cc"
// -*- C++ -*-
//
// EtaPiPiDefaultCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiDefaultCurrent class.
//

#include "EtaPiPiDefaultCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/epsilon.h"
using namespace Herwig;
using namespace ThePEG;

DescribeClass<EtaPiPiDefaultCurrent,WeakCurrent>
describeHerwigEtaPiPiDefaultCurrent("Herwig::EtaPiPiDefaultCurrent",
				    "HwWeakCurrents.so");

EtaPiPiDefaultCurrent::EtaPiPiDefaultCurrent() {
  // set up for the modes in the base class
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
  // the pion decay constant
  _fpi=130.7*MeV/sqrt(2.);
  _mpi=ZERO;
  // set the initial weights for the resonances
  // the rho weights
  _rhoF123wgts = { 1.0,-0.145,0.};
  _rhoF5wgts   = {-26.,   6.5,1.};
  // local values of the rho parameters
  _rhoF123masses = {0.773*GeV,1.370*GeV,1.750*GeV};
  _rhoF123widths = {0.145*GeV,0.510*GeV,0.120*GeV};
  _rhoF5masses   = {0.773*GeV,1.500*GeV,1.750*GeV};
  _rhoF5widths   = {0.145*GeV,0.220*GeV,0.120*GeV};
}

void EtaPiPiDefaultCurrent::doinit() {
  WeakCurrent::doinit();
  // the particles we will use a lot
  tPDPtr piplus(getParticleData(ParticleID::piplus));
  // masses for the running widths
  _mpi=piplus->mass();
}

void EtaPiPiDefaultCurrent::persistentOutput(PersistentOStream & os) const {
  os << _rhoF123wgts << _rhoF5wgts << ounit(_fpi,GeV) << ounit(_mpi,GeV)
     << ounit(_rhoF123masses,GeV) << ounit(_rhoF5masses,GeV) 
     << ounit(_rhoF123widths,GeV) << ounit(_rhoF5widths,GeV);
}

void EtaPiPiDefaultCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _rhoF123wgts >> _rhoF5wgts >> iunit(_fpi,GeV) >> iunit(_mpi,GeV)
     >> iunit(_rhoF123masses,GeV) >> iunit(_rhoF5masses,GeV) 
     >> iunit(_rhoF123widths,GeV) >> iunit(_rhoF5widths,GeV);
}

void EtaPiPiDefaultCurrent::Init() {
        
  static ClassDocumentation<EtaPiPiDefaultCurrent> documentation
    ("The EtaPiPiDefaultCurrent class is designed to implement "
     "the three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, " 
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "pi- Kbar0 pi0, pi- pi0 eta. It uses the same currents as those in TAUOLA.",
     "The three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, "
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "and pi- Kbar0 pi0, pi- pi0 eta "
     "use the same currents as \\cite{Jadach:1993hs,Kuhn:1990ad,Decker:1992kj}.",
     "%\\cite{Jadach:1993hs}\n"
     "\\bibitem{Jadach:1993hs}\n"
     "  S.~Jadach, Z.~Was, R.~Decker and J.~H.~Kuhn,\n"
     "  %``The Tau Decay Library Tauola: Version 2.4,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 76}, 361 (1993).\n"
     "  %%CITATION = CPHCB,76,361;%%\n"
     "%\\cite{Kuhn:1990ad}\n"
     "\\bibitem{Kuhn:1990ad}\n"
     "  J.~H.~Kuhn and A.~Santamaria,\n"
     "  %``Tau decays to pions,''\n"
     "  Z.\\ Phys.\\  C {\\bf 48}, 445 (1990).\n"
     "  %%CITATION = ZEPYA,C48,445;%%\n"
     "%\\cite{Decker:1992kj}\n"
     "\\bibitem{Decker:1992kj}\n"
     "  R.~Decker, E.~Mirkes, R.~Sauer and Z.~Was,\n"
     "  %``Tau decays into three pseudoscalar mesons,''\n"
     "  Z.\\ Phys.\\  C {\\bf 58}, 445 (1993).\n"
     "  %%CITATION = ZEPYA,C58,445;%%\n"
     );
  
  static ParVector<EtaPiPiDefaultCurrent,double> interfaceF123RhoWgt
    ("F123RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &EtaPiPiDefaultCurrent::_rhoF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<EtaPiPiDefaultCurrent,double> interfaceF5RhoWgt
    ("F5RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &EtaPiPiDefaultCurrent::_rhoF5wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<EtaPiPiDefaultCurrent,Energy> interfacerhoF123masses
    ("rhoF123masses",
     "The masses for the rho resonances if used local values",
     &EtaPiPiDefaultCurrent::_rhoF123masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<EtaPiPiDefaultCurrent,Energy> interfacerhoF123widths
    ("rhoF123widths",
     "The widths for the rho resonances if used local values",
     &EtaPiPiDefaultCurrent::_rhoF123widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<EtaPiPiDefaultCurrent,Energy> interfacerhoF5masses
    ("rhoF5masses",
     "The masses for the rho resonances if used local values",
     &EtaPiPiDefaultCurrent::_rhoF5masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<EtaPiPiDefaultCurrent,Energy> interfacerhoF5widths
    ("rhoF5widths",
     "The widths for the rho resonances if used local values",
     &EtaPiPiDefaultCurrent::_rhoF5widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<EtaPiPiDefaultCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &EtaPiPiDefaultCurrent::_fpi, MeV, 92.4*MeV, ZERO, 200.0*MeV,
     false, false, true);
}

// complete the construction of the decay mode for integration
bool EtaPiPiDefaultCurrent::createMode(int icharge, tcPDPtr resonance,
				       FlavourInfo flavour,
				       unsigned int imode,PhaseSpaceModePtr mode,
				       unsigned int iloc,int ires,
				       PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((imode==0 && abs(icharge)!=3) ||
     (imode>0  && icharge !=0)) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return false;
      break;
    case IsoSpin::I3One:
      if(imode==1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode==1 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero ) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       ) return false;
  // make sure that the decays are kinematically allowed
  int iq(0),ia(0);
  tPDVector part = particles(icharge,imode,iq,ia);
  tPDVector extpart(particles(1,imode,iq,ia));
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  // set up the resonances
  tPDPtr res[3];
  if(icharge==0) {
    res[0] =getParticleData(113);
    res[1] =getParticleData(100113);
    res[2] =getParticleData(30113);
  }
  else {
    res[0] =getParticleData(213);
    res[1] =getParticleData(100213);
    res[2] =getParticleData(30213);
    if(icharge==-3) {
      for(unsigned int ix=0;ix<3;++ix) {
  	if(res[ix]&&res[ix]->CC()) res[ix]=res[ix]->CC();
      }
    }
  }
  // channels for pi- pi0 eta
  for(unsigned int ix=0;ix<3;++ix) {
    if(resonance && resonance != res[ix]) continue;
    for(unsigned int iy=0;iy<3;++iy) {
      mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],ires+1,iloc+3,ires+1,res[iy],
			ires+2,iloc+1,ires+2,iloc+2));
    }
  }
  // reset the rho masses
  for(unsigned int ix=0;ix<_rhoF5masses.size();++ix)
    mode->resetIntermediate(res[ix],_rhoF5masses[ix],_rhoF5widths[ix]);
  return true;
}

void EtaPiPiDefaultCurrent::dataBaseOutput(ofstream & output,bool header,
					      bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::EtaPiPiDefaultCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  output << "newdef " << name() << ":FPi "     << _fpi/MeV     << "\n";
  for(unsigned int ix=0;ix<_rhoF123wgts.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F123RhoWeight " << ix << " " << _rhoF123wgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF5wgts.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F5RhoWeight " << ix << " " << _rhoF5wgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF123masses.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF123masses " << ix 
	   << " " << _rhoF123masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF123widths.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF123widths " << ix << " " 
	   << _rhoF123widths[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF5masses.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF5masses " << ix << " " 
	   << _rhoF5masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF5widths.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF5widths " << ix << " " 
	   << _rhoF5widths[ix]/GeV << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
EtaPiPiDefaultCurrent::current(tcPDPtr resonance,
			       FlavourInfo flavour,
			       const int imode, const int ichan, Energy & scale, 
			       const tPDVector & outgoing,
			       const vector<Lorentz5Momentum> & momenta,
			       DecayIntegrator::MEOption) const {
  useMe();
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IOne)
    return vector<LorentzPolarizationVectorE>();
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode==1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode==1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return vector<LorentzPolarizationVectorE>();
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q = momenta[0] + momenta[1] + momenta[2];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 s3 = (momenta[0]+momenta[1]).m2();
  // the form factor
  Complex F5(0.);
  int ires1(-1),ires2(-1);
  if(ichan>=0) {
    ires1 = ichan/3;
    ires2 = ichan%3;
  }
  else {
    if(resonance) {
      switch(resonance->id()/1000) {
      case 0:
	ires1 = 0;
	break;
      case 100:
	ires1 = 1;
	break;
      case 30 :
	ires1 = 2;
	break;
      default:
	assert(false);
      }
    }
  }
  F5 = BrhoF5(q2,ires1)*BrhoF123(s3,ires2)*sqrt(2./3.);
  // constant the current
  LorentzPolarizationVector vect = -Complex(0.,1.)*F5/sqr(Constants::twopi)/pow<3,1>(_fpi)*
    Helicity::epsilon(momenta[0],momenta[1],momenta[2]);
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,q.mass()*vect);
}

bool EtaPiPiDefaultCurrent::accept(vector<int> id) {
  // check there are only three particles
  if(id.size()!=3) return false;
  unsigned int npip(0),npim(0),npi0(0),neta(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::eta)     ++neta;
  }
  if( (npip==1&&npim==1&&neta==1) ||
      (npi0==1&&npim+npip==1&&neta==1))
    return true;
  else
    return false;
}

unsigned int EtaPiPiDefaultCurrent::decayMode(vector<int> id) {
  unsigned int npi(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::piplus) ++npi;
  }
  if(npi==2) return 1;
  else       return 0;
}

tPDVector EtaPiPiDefaultCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector output(3);
  output[0]=getParticleData(ParticleID::piplus);
  output[2]=getParticleData(ParticleID::eta);
  if(imode==0) {
    output[1]=getParticleData(ParticleID::pi0);
    if(icharge==-3) {
      for(unsigned int ix=0;ix<output.size();++ix) {
	if(output[ix]->CC()) output[ix]=output[ix]->CC();
      }
    }
  }
  else {
    output[1]=getParticleData(ParticleID::piminus);
  }
  return output;
}
#line 1 "./ThreePionCLEOCurrent.cc"
// -*- C++ -*-
//
// ThreePionCLEOCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreePionCLEOCurrent class.
//

#include "ThreePionCLEOCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/ResonanceHelpers.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeClass<ThreePionCLEOCurrent,WeakCurrent>
describeHerwigThreePionCLEOCurrent("Herwig::ThreePionCLEOCurrent",
				   "HwWeakCurrents.so");
HERWIG_INTERPOLATOR_CLASSDESC(ThreePionCLEOCurrent,Energy,Energy2)

ThreePionCLEOCurrent::ThreePionCLEOCurrent() {
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(2,-1);
  setInitialModes(6);
  // rho masses and widths
  _rhomass  = {0.7743*GeV,1.370*GeV};
  _rhowidth = {0.1491*GeV,0.386*GeV};
  // f_2 mass and width
  _f2mass  = 1.275*GeV;
  _f2width = 0.185*GeV;
  // f_0(1370) mass and width
  _f0mass  = 1.186*GeV;
  _f0width = 0.350*GeV;
  // sigma mass and width
  _sigmamass  = 0.860*GeV;
  _sigmawidth = 0.880*GeV;
  // a1 mass and width
  _a1mass  = 1.331*GeV;
  _a1width = 0.814*GeV;
  // parameters for the K K* contribution to the a_1 running width
  _mKstar = 894*MeV;
  _mK     = 496*MeV;
  _gammk  = 3.32;
  // pion decay constant
  _fpi = 130.7*MeV/sqrt(2.);
  // couplings and phases for the different channels
  // p-wave rho and rho prime
  using Constants::pi;
  _rhomagP   = {1.,0.12};
  _rhophaseP = {0.,0.99*pi};
  // d-wave rho and rho prime
  _rhomagD   = {0.37/GeV2,0.87/GeV2};
  _rhophaseD = {-0.15*pi, 0.53*pi};
  // f_2
  _f2mag   = 0.71/GeV2;
  _f2phase = 0.56*pi;
  _f2coup  = ZERO;
  // sigma
  _sigmamag   = 2.10;
  _sigmaphase = 0.23*pi;
  _sigmacoup  = 0.;
  // f_0
  _f0mag   = 0.77;
  _f0phase = -0.54*pi;
  _f0coup  = 0.;
  // initialize the a_1 width
  _initializea1=false;
  _a1opt=true;
  double  a1q2in[200]={0      ,15788.6,31577.3,47365.9,63154.6,78943.2,
		       94731.9,110521 ,126309 ,142098 ,157886 ,173675 ,
		       189464 ,205252 ,221041 ,236830 ,252618 ,268407 ,
		       284196 ,299984 ,315773 ,331562 ,347350 ,363139 ,
		       378927 ,394716 ,410505 ,426293 ,442082 ,457871 ,
		       473659 ,489448 ,505237 ,521025 ,536814 ,552603 ,
		       568391 ,584180 ,599969 ,615757 ,631546 ,647334 ,
		       663123 ,678912 ,694700 ,710489 ,726278 ,742066 ,
		       757855 ,773644 ,789432 ,805221 ,821010 ,836798 ,
		       852587 ,868375 ,884164 ,899953 ,915741 ,931530 ,
		       947319 ,963107 ,978896 ,994685 ,
		       1.01047e+06,1.02626e+06,1.04205e+06,1.05784e+06,
		       1.07363e+06,1.08942e+06,1.10521e+06,1.12099e+06,
		       1.13678e+06,1.15257e+06,1.16836e+06,1.18415e+06,
		       1.19994e+06,1.21573e+06,1.23151e+06,1.24730e+06,
		       1.26309e+06,1.27888e+06,1.29467e+06,1.31046e+06,
		       1.32625e+06,1.34203e+06,1.35782e+06,1.37361e+06,
		       1.38940e+06,1.40519e+06,1.42098e+06,1.43677e+06,
		       1.45256e+06,1.46834e+06,1.48413e+06,1.49992e+06,
		       1.51571e+06,1.53150e+06,1.54729e+06,1.56308e+06,
		       1.57886e+06,1.59465e+06,1.61044e+06,1.62623e+06,
		       1.64202e+06,1.65781e+06,1.67360e+06,1.68939e+06,
		       1.70517e+06,1.72096e+06,1.73675e+06,1.75254e+06,
		       1.76833e+06,1.78412e+06,1.79991e+06,1.81569e+06,
		       1.83148e+06,1.84727e+06,1.86306e+06,1.87885e+06,
		       1.89464e+06,1.91043e+06,1.92621e+06,1.94200e+06,
		       1.95779e+06,1.97358e+06,1.98937e+06,2.00516e+06,
		       2.02095e+06,2.03674e+06,2.05252e+06,2.06831e+06,
		       2.08410e+06,2.09989e+06,2.11568e+06,2.13147e+06,
		       2.14726e+06,2.16304e+06,2.17883e+06,2.19462e+06,
		       2.21041e+06,2.22620e+06,2.24199e+06,2.25778e+06,
		       2.27356e+06,2.28935e+06,2.30514e+06,2.32093e+06,
		       2.33672e+06,2.35251e+06,2.36830e+06,2.38409e+06,
		       2.39987e+06,2.41566e+06,2.43145e+06,2.44724e+06,
		       2.46303e+06,2.47882e+06,2.49461e+06,2.51039e+06,
		       2.52618e+06,2.54197e+06,2.55776e+06,2.57355e+06,
		       2.58934e+06,2.60513e+06,2.62092e+06,2.63670e+06,
		       2.65249e+06,2.66828e+06,2.68407e+06,2.69986e+06,
		       2.71565e+06,2.73144e+06,2.74722e+06,2.76301e+06,
		       2.77880e+06,2.79459e+06,2.81038e+06,2.82617e+06,
		       2.84196e+06,2.85774e+06,2.87353e+06,2.88932e+06,
		       2.90511e+06,2.92090e+06,2.93669e+06,2.95248e+06,
		       2.96827e+06,2.98405e+06,2.99984e+06,3.01563e+06,
		       3.03142e+06,3.04721e+06,3.06300e+06,3.07879e+06,
		       3.09457e+06,3.11036e+06,3.12615e+06,3.14194e+06};
  double a1widthin[200]={0,0,0,0,0,0,0,0,
			 0,0,0,0.00021256,0.0107225,0.0554708,0.150142,0.303848,
			 0.522655,0.81121,1.1736,1.61381,2.13606,2.74499,3.44583,4.24454,
			 5.14795,6.16391,7.3014,8.57079,9.98398,11.5547,13.2987,15.2344,
			 17.3827,19.7683,22.4195,25.3695,28.6568,32.3264,36.4311,41.0322,
			 46.201,52.0203,58.5847,66.0011,74.3871,83.8666,94.5615,106.578,
			 119.989,134.807,150.968,168.315,186.615,205.576,224.893,244.28,
			 263.499,282.364,300.748,318.569,335.781,352.367,368.327,383.677,
			 398.438,412.638,426.306,439.472,452.167,464.421,476.263,487.719,
			 498.815,509.576,520.024,530.179,540.063,549.693,559.621,568.26,
			 577.229,586.005,594.604,603.035,611.314,619.447,627.446,635.321,
			 643.082,650.736,658.288,665.75,673.127,680.427,687.659,694.82,
			 701.926,708.977,715.983,722.944,729.862,736.752,743.619,750.452,
			 757.271,764.076,770.874,777.658,784.444,791.233,798.027,804.838,
			 811.649,818.485,825.342,832.224,839.139,846.082,853.059,860.079,
			 867.143,874.248,881.409,919.527,945.28,965.514,983.228,999.471,
			 1014.69,1029.15,1043.05,1056.49,1069.57,1082.36,1094.88,1107.2,
			 1120.89,1131.4,1143.33,1155.15,1166.92,1178.61,1190.27,1201.92,
			 1213.55,1225.18,1236.81,1250.06,1260.16,1271.86,1283.64,1295.46,
			 1307.36,1319.3,1331.34,1343.45,1355.64,1367.93,1380.31,1392.77,
			 1405.35,1418.03,1430.83,1443.75,1457.17,1469.94,1483.22,1496.64,
			 1510.18,1523.86,1537.67,1551.64,1565.72,1579.99,1594.38,1608.92,
			 1623.63,1642.08,1653.51,1668.69,1684.03,1699.53,1715.21,1731.04,
			 1747.05,1763.23,1779.59,1796.12,1812.83,1829.72,1846.79,1864.04,
			 1881.49,1899.11,1916.93,1934.93,1953.13,1971.52,1990.12,2008.89};
  if(_a1runwidth.empty()) {
    vector<double> tmp1(a1widthin,a1widthin+200);
    std::transform(tmp1.begin(), tmp1.end(),
		   back_inserter(_a1runwidth),
		   [](double x){return x*GeV;});

    vector<double> tmp2(a1q2in,a1q2in+200);
    _a1runq2.clear();
    std::transform(tmp2.begin(), tmp2.end(),
		   back_inserter(_a1runq2),
		   [](double x){return x*GeV2;});
  }
  // zero parameters which will be calculated later to avoid problems
  _mpi0=ZERO;
  _mpic=ZERO;
  _fact=ZERO;
  _maxmass=ZERO;
  _maxcalc=ZERO;
}

void ThreePionCLEOCurrent::doinit() {
  WeakCurrent::doinit();
  // parameters for the breit-wigners
  _mpic = getParticleData(ParticleID::piplus)->mass();
  _mpi0 = getParticleData(ParticleID::pi0)   ->mass();
  // couplings for the different modes
  Complex ii(0.,1.);
  _rhocoupP.resize(_rhomagP.size());
  for(unsigned int ix=0;ix<_rhomagP.size();++ix)
    _rhocoupP[ix]=_rhomagP[ix]*(cos(_rhophaseP[ix])+ii*sin(_rhophaseP[ix]));
  _rhocoupD.resize(_rhomagD.size());
  for(unsigned int ix=0;ix<_rhomagD.size();++ix)
    _rhocoupD[ix]=_rhomagD[ix]*(cos(_rhophaseD[ix])+ii*sin(_rhophaseD[ix]));
  _f0coup=_f0mag*(cos(_f0phase)+ii*sin(_f0phase));
  _f2coup=_f2mag*(cos(_f2phase)+ii*sin(_f2phase));
  _sigmacoup=_sigmamag*(cos(_sigmaphase)+ii*sin(_sigmaphase));
  // overall coupling
  _fact = 2.*sqrt(2.)/_fpi/3.;
  // initialise the a_1 running width calculation
  inita1Width(-1);
}

void ThreePionCLEOCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(_rhomass,GeV) << ounit(_rhowidth,GeV)
     << ounit(_f2mass,GeV) << ounit(_f2width,GeV)
     << ounit(_f0mass,GeV) << ounit(_f0width,GeV)
     << ounit(_sigmamass,GeV) << ounit(_sigmawidth,GeV)
     << ounit(_mpi0,GeV) << ounit(_mpic,GeV)
     << ounit(_fpi,GeV) << ounit(_fact,1/GeV)
     << _rhomagP << _rhophaseP
     << _rhocoupP << ounit(_rhomagD,1/GeV2) << _rhophaseD
     << ounit(_rhocoupD,1/GeV2) <<ounit(_f2mag,1/GeV2) << _f2phase << ounit(_f2coup,1/GeV2)
     << _f0mag << _f0phase << _f0coup << _sigmamag << _sigmaphase << _sigmacoup
     << ounit(_a1mass,GeV) << ounit(_a1width,GeV) << ounit(_a1runwidth,GeV)
     << ounit(_a1runq2,GeV2) <<  _initializea1
     << ounit(_mKstar,GeV) << ounit(_mK,GeV) << _gammk << _a1opt
     << ounit(_maxmass,GeV) << ounit(_maxcalc,GeV) << _a1runinter;
}

void ThreePionCLEOCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_rhomass,GeV) >> iunit(_rhowidth,GeV)
     >> iunit(_f2mass,GeV) >> iunit(_f2width,GeV)
     >> iunit(_f0mass,GeV) >> iunit(_f0width,GeV)
     >> iunit(_sigmamass,GeV) >> iunit(_sigmawidth,GeV)
     >> iunit(_mpi0,GeV) >> iunit(_mpic,GeV)
     >> iunit(_fpi,GeV) >> iunit(_fact,1/GeV)
     >> _rhomagP >> _rhophaseP
     >> _rhocoupP >> iunit(_rhomagD,1/GeV2) >> _rhophaseD >> iunit(_rhocoupD,1/GeV2)
     >> iunit(_f2mag,1/GeV2) >> _f2phase >> iunit(_f2coup,1/GeV2)
     >> _f0mag >> _f0phase >> _f0coup >> _sigmamag >> _sigmaphase >> _sigmacoup
     >> iunit(_a1mass,GeV) >> iunit(_a1width,GeV) >>  iunit(_a1runwidth,GeV)
     >> iunit(_a1runq2,GeV2) >>  _initializea1
     >> iunit(_mKstar,GeV) >> iunit(_mK,GeV) >> _gammk >> _a1opt
     >> iunit(_maxmass,GeV) >> iunit(_maxcalc,GeV) >> _a1runinter;
}

void ThreePionCLEOCurrent::Init() {

  static ClassDocumentation<ThreePionCLEOCurrent> documentation
    ("The ThreePionCLEOCurrent class performs the decay of the"
     " tau to three pions using the currents from CLEO",
     "The decay of tau to three pions is modelled using the currents from "
     "\\cite{Asner:1999kj}.",
     "  %\\cite{Asner:1999kj}\n"
     "\\bibitem{Asner:1999kj}\n"
     "  D.~M.~Asner {\\it et al.}  [CLEO Collaboration],\n"
     "   ``Hadronic structure in the decay tau- --> nu/tau pi- pi0 pi0 and the  sign\n"
     "  %of the tau neutrino helicity,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 61}, 012002 (2000)\n"
     "  [arXiv:hep-ex/9902022].\n"
     "  %%CITATION = PHRVA,D61,012002;%%\n"
     );

  static ParVector<ThreePionCLEOCurrent,Energy> interfacerhomass
    ("RhoMasses",
     "The masses of the different rho resonnaces",
     &ThreePionCLEOCurrent::_rhomass,
     MeV, 0, ZERO, -10000*MeV, 10000*MeV, false, false, true);

  static ParVector<ThreePionCLEOCurrent,Energy> interfacerhowidth
    ("RhoWidths",
     "The widths of the different rho resonnaces",
     &ThreePionCLEOCurrent::_rhowidth,
     MeV, 0, ZERO, -10000*MeV, 10000*MeV, false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacef_2Mass
    ("f_2Mass",
     "The mass of the f_2 meson",
     &ThreePionCLEOCurrent::_f2mass, GeV, 1.275*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacef_2Width
    ("f_2Width",
     "The width of the f_2 meson",
     &ThreePionCLEOCurrent::_f2width, GeV, 0.185*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacef_0Mass
    ("f_0Mass",
     "The mass of the f_0 meson",
     &ThreePionCLEOCurrent::_f0mass, GeV, 1.186*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacef_0Width
    ("f_0Width",
     "The width of the f_0 meson",
     &ThreePionCLEOCurrent::_f0width, GeV, 0.350*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacesigmaMass
    ("sigmaMass",
     "The mass of the sigma meson",
     &ThreePionCLEOCurrent::_sigmamass, GeV, 0.860*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacesigmaWidth
    ("sigmaWidth",
     "The width of the sigma meson",
     &ThreePionCLEOCurrent::_sigmawidth, GeV, 0.880*GeV, ZERO, 2.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacea1Mass
    ("a1Mass",
     "The mass of the a_1 meson",
     &ThreePionCLEOCurrent::_a1mass, GeV, 1.331*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacea1Width
    ("a1Width",
     "The width of the a_1 meson",
     &ThreePionCLEOCurrent::_a1width, GeV, 0.814*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfaceKaonMass
    ("KaonMass",
     "The mass of the kaon",
     &ThreePionCLEOCurrent::_mK, GeV, 0.496*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfaceKStarMass
    ("KStarMass",
     "The mass of the k* meson",
     &ThreePionCLEOCurrent::_mKstar, GeV, 0.894*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfaceKaonCoupling
    ("KaonCoupling",
     "The relative coupling for the kaon in the a_1 running width",
     &ThreePionCLEOCurrent::_gammk, 3.32, 0.0, 10.0,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant",
     &ThreePionCLEOCurrent::_fpi, MeV, 130.7*MeV/sqrt(2.), ZERO, 500.0*MeV,
     false, false, true);

  static ParVector<ThreePionCLEOCurrent,double> interfacerhomagP
    ("RhoPWaveMagnitude",
     "The magnitude of the couplings for the p-wave rho currents",
     &ThreePionCLEOCurrent::_rhomagP,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<ThreePionCLEOCurrent,double> interfacerhophaseP
    ("RhoPWavePhase",
     "The phase of the couplings for the p-wave rho currents",
     &ThreePionCLEOCurrent::_rhophaseP,
     0, 0, 0, -Constants::twopi, Constants::twopi, false, false, true);

  static ParVector<ThreePionCLEOCurrent,InvEnergy2> interfacerhomagD
    ("RhoDWaveMagnitude",
     "The magnitude of the couplings for the d-wave rho currents",
     &ThreePionCLEOCurrent::_rhomagD,
     1/MeV2, 0, ZERO, ZERO, 10000/MeV2, false, false, true);

  static ParVector<ThreePionCLEOCurrent,double> interfacerhophaseD
    ("RhoDWavePhase",
     "The phase of the couplings for the d-wave rho currents",
     &ThreePionCLEOCurrent::_rhophaseD,
     0, 0, 0, -Constants::twopi, Constants::twopi, false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacef0Phase
    ("f0Phase",
     "The phase of the f_0 scalar current",
     &ThreePionCLEOCurrent::_f0phase, 0.54*Constants::pi, -Constants::twopi, Constants::twopi,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacef2Phase
    ("f2Phase",
     "The phase of the f_2 tensor current",
     &ThreePionCLEOCurrent::_f2phase, 0.56*Constants::pi,-Constants::twopi, Constants::twopi,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacesigmaPhase
    ("sigmaPhase",
     "The phase of the sigma scalar current",
     &ThreePionCLEOCurrent::_sigmaphase, 0.23*Constants::pi, -Constants::twopi, Constants::twopi,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacef0Magnitude
    ("f0Magnitude",
     "The magnitude of the f_0 scalar current",
     &ThreePionCLEOCurrent::_f0mag, 0.77, 0.0, 10,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,InvEnergy2> interfacef2Magnitude
    ("f2Magnitude",
     "The magnitude of the f_2 tensor current",
     &ThreePionCLEOCurrent::_f2mag, 1./GeV2, 0.71/GeV2, ZERO, 10./GeV2,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacesigmaMagnitude
    ("sigmaMagnitude",
     "The magnitude of the sigma scalar current",
     &ThreePionCLEOCurrent::_sigmamag, 2.1, 0.0, 10,
     false, false, true);

  static ParVector<ThreePionCLEOCurrent,Energy> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &ThreePionCLEOCurrent::_a1runwidth,
     MeV, 0, ZERO, ZERO, 10000000*MeV, false, false, true);

  static ParVector<ThreePionCLEOCurrent,Energy2> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &ThreePionCLEOCurrent::_a1runq2,
     MeV2, 0, ZERO, ZERO, 10000000*MeV2, false, false, true);

  static Switch<ThreePionCLEOCurrent,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &ThreePionCLEOCurrent::_initializea1, false, false, false);
  static SwitchOption interfaceInitializea1Initialization
    (interfaceInitializea1,
     "Yes",
     "Initialize the calculation",
     true);
  static SwitchOption interfaceInitializea1NoInitialization
    (interfaceInitializea1,
     "No",
     "Use the default values",
     false);

  static Switch<ThreePionCLEOCurrent,bool> interfacea1WidthOption
    ("a1WidthOption",
     "Option for the treatment of the a1 width",
     &ThreePionCLEOCurrent::_a1opt, true, false, false);
  static SwitchOption interfacea1WidthOptionLocal
    (interfacea1WidthOption,
     "Local",
     "Use a calculation of the running width based on the parameters as"
     " interpolation table.",
     true);
  static SwitchOption interfacea1WidthOptionParam
    (interfacea1WidthOption,
     "Kuhn",
     "Use the parameterization of Kuhn and Santamaria for default parameters."
     " This should only be used for testing vs TAUOLA",
     false);
}

// initialisation of the a_1 width
// (iopt=-1 initialises, iopt=0 starts the interpolation)
void ThreePionCLEOCurrent::inita1Width(int iopt) {
  if(iopt==-1) {
    _maxcalc=_maxmass;
    if(!_initializea1||_maxmass==ZERO) return;
    // parameters for the table of values
    Energy2 step=sqr(_maxmass)/200.;
    // function to be integrated to give the matrix element
    // integrator to perform the integral
    vector<double> inweights;inweights.push_back(0.5);inweights.push_back(0.5);
    vector<int> intype;intype.push_back(2);intype.push_back(3);
    Energy mrho=getParticleData(ParticleID::rhoplus)->mass();
    Energy wrho=getParticleData(ParticleID::rhoplus)->width();
    vector<Energy> inmass;inmass.push_back(mrho);inmass.push_back(mrho);
    vector<Energy> inwidth;inwidth.push_back(wrho);inwidth.push_back(wrho);
    vector<double> inpow(2,0.0);
    ThreeBodyAllOnCalculator<ThreePionCLEOCurrent>
      widthgenN(inweights,intype,inmass,inwidth,inpow,*this,0,_mpi0,_mpi0,_mpic);
    ThreeBodyAllOnCalculator<ThreePionCLEOCurrent>
      widthgenC(inweights,intype,inmass,inwidth,inpow,*this,1,_mpic,_mpic,_mpic);
    // normalisation constant to give physical width if on shell
    double a1const = _a1width/(widthgenN.partialWidth(sqr(_a1mass))+
			       widthgenC.partialWidth(sqr(_a1mass)));
    // loop to give the values
    _a1runq2.clear();_a1runwidth.clear();
    for(Energy2 moff2=ZERO; moff2<=sqr(_maxmass); moff2+=step) {
      Energy moff=sqrt(moff2);
      _a1runq2.push_back(moff2);
      Energy charged=a1const*widthgenC.partialWidth(moff2);
      Energy neutral=a1const*widthgenN.partialWidth(moff2);
      Energy kaon = moff<=_mK+_mKstar ? ZERO : 2.870*_gammk*_gammk/8./Constants::pi*
	Kinematics::pstarTwoBodyDecay(moff,_mK,_mKstar)/moff2*GeV2;
      Energy total = charged + neutral + kaon;
      _a1runwidth.push_back(total);
    }
  }
  // set up the interpolator
  else if(iopt==0) {
    _a1runinter = make_InterpolatorPtr(_a1runwidth,_a1runq2,3);
  }
}

void ThreePionCLEOCurrent::CLEOFormFactor(int imode,int ichan,
					  Energy2 q2,Energy2 s1, Energy2 s2, Energy2 s3,
					  Complex & F1, Complex & F2,
					  Complex & F3) const {
  useMe();
  double fact=1.;
  if(imode<=1) {
    // identical particle factors
    fact = 1./sqrt(6.);
    // compute the breit wigners we need
    Complex sigbws1 = Resonance::BreitWignerSWave(s1,_sigmamass,_sigmawidth,_mpi0,_mpi0);
    Complex sigbws2 = Resonance::BreitWignerSWave(s2,_sigmamass,_sigmawidth,_mpi0,_mpi0);
    Complex sigbws3 = Resonance::BreitWignerSWave(s3,_sigmamass,_sigmawidth,_mpi0,_mpi0);
    Complex f0bws1  = Resonance::BreitWignerSWave(s1,   _f0mass,   _f0width,_mpi0,_mpi0);
    Complex f0bws2  = Resonance::BreitWignerSWave(s2,   _f0mass,   _f0width,_mpi0,_mpi0);
    Complex f0bws3  = Resonance::BreitWignerSWave(s3,   _f0mass,   _f0width,_mpi0,_mpi0);
    Complex f2bws1  = Resonance::BreitWignerDWave(s1,   _f2mass,   _f2width,_mpi0,_mpi0);
    Complex f2bws2  = Resonance::BreitWignerDWave(s2,   _f2mass,   _f2width,_mpi0,_mpi0);
    Complex f2bws3  = Resonance::BreitWignerDWave(s3,   _f2mass,   _f2width,_mpi0,_mpi0);
    if(ichan<0) {
      // the scalar terms
      F1=2./3.*(_sigmacoup*sigbws3+_f0coup*f0bws3)
	-2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
      F2=2./3.*(_sigmacoup*sigbws3+_f0coup*f0bws3)
	-2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1);
      F3=-2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1)
	+2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
      // the tensor terms
      complex<Energy2> Dfact1 = 1./18.*(4.*_mpi0*_mpi0-s1)*(q2+s1-_mpi0*_mpi0)/s1*f2bws1;
      complex<Energy2> Dfact2 = 1./18.*(4.*_mpi0*_mpi0-s2)*(q2+s2-_mpi0*_mpi0)/s2*f2bws2;
      complex<Energy2> Dfact3 = 1./18.*(4.*_mpi0*_mpi0-s3)*(q2-_mpi0*_mpi0+s3)/s3*f2bws3;
      F1+=Complex(_f2coup*( 0.5*(s3-s2)*f2bws1-Dfact2+Dfact3));
      F2+=Complex(_f2coup*( 0.5*(s3-s1)*f2bws2-Dfact1+Dfact3));
      F3+=Complex(_f2coup*(-0.5*(s1-s2)*f2bws3-Dfact1+Dfact2));
    }
    else if(ichan==0) {
      F2=-2./3.*_sigmacoup*sigbws1;
      F3=-2./3.*_sigmacoup*sigbws1;
    }
    else if(ichan==1) {
      F1=-2./3.*_sigmacoup*sigbws2;
      F3=+2./3.*_sigmacoup*sigbws2;
    }
    else if(ichan==2) {
      F1= 2./3.*_sigmacoup*sigbws3;
      F2= 2./3.*_sigmacoup*sigbws3;
    }
    else if(ichan==3) {
      complex<Energy2> Dfact1 = 1./18.*(4.*_mpi0*_mpi0-s1)*(q2+s1-_mpi0*_mpi0)/s1*f2bws1;
      F1+=Complex(_f2coup*0.5*(s3-s2)*f2bws1);
      F2-=Complex(_f2coup*Dfact1);
      F3-=Complex(_f2coup*Dfact1);
    }
    else if(ichan==4) {
      complex<Energy2> Dfact2 = 1./18.*(4.*_mpi0*_mpi0-s2)*(q2+s2-_mpi0*_mpi0)/s2*f2bws2;
      F2+=Complex(_f2coup*0.5*(s3-s1)*f2bws2);
      F1-=Complex(_f2coup*Dfact2);
      F3+=Complex(_f2coup*Dfact2);
    }
    else if(ichan==5) {
      complex<Energy2> Dfact3 = 1./18.*(4.*_mpi0*_mpi0-s3)*(q2-_mpi0*_mpi0+s3)/s3*f2bws3;
      F3+=Complex(-_f2coup*0.5*(s1-s2)*f2bws3);
      F1+=Complex(_f2coup*Dfact3);
      F2+=Complex(_f2coup*Dfact3);
    }
    else if(ichan==6) {
      F2=-2./3.*_f0coup*f0bws1;
      F3=-2./3.*_f0coup*f0bws1;
    }
    else if(ichan==7) {
      F1=-2./3.*_f0coup*f0bws2;
      F3=+2./3.*_f0coup*f0bws2;
    }
    else if(ichan==8) {
      F1= 2./3.*_f0coup*f0bws3;
      F2= 2./3.*_f0coup*f0bws3;
    }
  }
  // calculate the pi0 pi0 pi+ factor
  else if(imode==2) {
    // identical particle factors
    fact = 1./sqrt(2.);
    // compute the breit wigners we need
    Complex rhos1bw[3],rhos2bw[3];
    for(unsigned int ix=0,N=max(_rhocoupP.size(),_rhocoupD.size());ix<N;++ix) {
      rhos1bw[ix] = Resonance::BreitWignerPWave(s1,_rhomass[ix], _rhowidth[ix],_mpic,_mpi0);
      rhos2bw[ix] = Resonance::BreitWignerPWave(s2,_rhomass[ix], _rhowidth[ix],_mpic,_mpi0);
    }
    Complex f0bw  = Resonance::BreitWignerSWave(s3,   _f0mass,   _f0width,_mpi0,_mpi0);
    Complex sigbw = Resonance::BreitWignerSWave(s3,_sigmamass,_sigmawidth,_mpi0,_mpi0);
    Complex f2bw  = Resonance::BreitWignerDWave(s3,   _f2mass,   _f2width,_mpi0,_mpi0);
    if(ichan<0) {
      // the p-wave rho terms
      for(unsigned int ix=0;ix<_rhocoupP.size();++ix) {
	F1+=_rhocoupP[ix]*rhos1bw[ix];
	F2+=_rhocoupP[ix]*rhos2bw[ix];
      }
      // the D-wave rho terms
      Energy2 Dfact1=-1./3.*((s3-_mpic*_mpic)-(s1-_mpi0*_mpi0));
      Energy2 Dfact2=-1./3.*((s3-_mpic*_mpic)-(s2-_mpi0*_mpi0));
      for(unsigned int ix=0;ix<_rhocoupD.size();++ix) {
	F1+=Complex(Dfact1*_rhocoupD[ix]*rhos2bw[ix]);
	F2+=Complex(Dfact2*_rhocoupD[ix]*rhos1bw[ix]);
	F3+=Complex(_rhocoupD[ix]*(Dfact2*rhos1bw[ix]-Dfact1*rhos2bw[ix]));
      }
      // the scalar terms
      Complex scalar=2./3.*(_sigmacoup*sigbw+_f0coup*f0bw);
      F1+=scalar;F2+=scalar;
      // the tensor terms
      Complex Dfact3=1./18./s3*_f2coup*(q2-_mpic*_mpic+s3)*(4.*_mpi0*_mpi0-s3)*f2bw;
      F1+=Dfact3;F2+=Dfact3;
      F3-=Complex(0.5*_f2coup*(s1-s2)*f2bw);
    }
    else if(ichan%2==0&&ichan<=4) {
      unsigned int ires=ichan/2;
      if(ires<_rhocoupP.size()){F1+=_rhocoupP[ires]*rhos1bw[ires];}
      Energy2 Dfact2=-1./3.*((s3-_mpic*_mpic)-(s2-_mpi0*_mpi0));
      if(ires<_rhocoupD.size()) {
	F2+=Complex(Dfact2*_rhocoupD[ires]*rhos1bw[ires]);
	F3+=Complex(_rhocoupD[ires]*Dfact2*rhos1bw[ires]);
      }
    }
    else if(ichan%2==1&&ichan<=5) {
      unsigned int ires=(ichan-1)/2;
      if(ires<_rhocoupP.size()){F2+=_rhocoupP[ires]*rhos2bw[ires];}
      Energy2 Dfact1=-1./3.*((s3-_mpic*_mpic)-(s1-_mpi0*_mpi0));
      if(ires<_rhocoupD.size()) {
	F1+=Complex(Dfact1*_rhocoupD[ires]*rhos2bw[ires]);
	F3-=Complex(_rhocoupD[ires]*Dfact1*rhos2bw[ires]);
      }
    }
    else if(ichan==6) {
      F1+=2./3.*_sigmacoup*sigbw;
      F2+=2./3.*_sigmacoup*sigbw;
    }
    else if(ichan==7) {
      Complex Dfact3=1./18./s3*_f2coup*(q2-_mpic*_mpic+s3)*(4.*_mpi0*_mpi0-s3)*f2bw;
      F1+=Dfact3;
      F2+=Dfact3;
      F3-=Complex(0.5*_f2coup*(s1-s2)*f2bw);
    }
    else if(ichan==8) {
      F1+=2./3.*_f0coup*f0bw;
      F2+=2./3.*_f0coup*f0bw;
    }
  }
  // a_1^0 ->pi+pi-pi0
  else if(imode==3||imode==4) {
    // compute the breit wigners we need
    Complex rhos1bw[3],rhos2bw[3];
    for(unsigned int ix=0,N=max(_rhocoupP.size(),_rhocoupD.size());ix<N;++ix) {
      rhos1bw[ix] = Resonance::BreitWignerPWave(s1,_rhomass[ix], _rhowidth[ix],_mpic,_mpi0);
      rhos2bw[ix] = Resonance::BreitWignerPWave(s2,_rhomass[ix], _rhowidth[ix],_mpic,_mpi0);
    }
    Complex f0bw  = Resonance::BreitWignerSWave(s3,   _f0mass,   _f0width,_mpic,_mpic);
    Complex sigbw = Resonance::BreitWignerSWave(s3,_sigmamass,_sigmawidth,_mpic,_mpic);
    Complex f2bw  = Resonance::BreitWignerDWave(s3,   _f2mass,   _f2width,_mpic,_mpic);
    if(ichan<0) {
      // the p-wave rho terms
      for(unsigned int ix=0;ix<_rhocoupP.size();++ix) {
	F1+=_rhocoupP[ix]*rhos1bw[ix];
	F2+=_rhocoupP[ix]*rhos2bw[ix];
      }
      // the D-wave rho terms
      Energy2 Dfact1=-1./3.*(s3-_mpi0*_mpi0-s1+_mpic*_mpic);
      Energy2 Dfact2=-1./3.*(s3-_mpi0*_mpi0-s2+_mpic*_mpic);
      for(unsigned int ix=0;ix<_rhocoupD.size();++ix) {
	F1+=Complex(Dfact1*_rhocoupD[ix]*rhos2bw[ix]);
	F2+=Complex(Dfact2*_rhocoupD[ix]*rhos1bw[ix]);
	F3+=Complex(_rhocoupD[ix]*(Dfact2*rhos1bw[ix]-Dfact1*rhos2bw[ix]));
      }
      // the scalar terms
      Complex scalar=2./3.*(_sigmacoup*sigbw+_f0coup*f0bw);
      F1+=scalar;
      F2+=scalar;
      // the tensor terms
      Complex Dfact3=1./18./s3*_f2coup*(q2-_mpi0*_mpi0+s3)*(4.*_mpic*_mpic-s3)*f2bw;
      F1+=Dfact3;
      F2+=Dfact3;
      F3-=Complex(0.5*_f2coup*(s1-s2)*f2bw);
    }
    else if(ichan%2==0&&ichan<=4) {
      unsigned int ires=ichan/2;
      if(ires<_rhocoupP.size()) F1+=_rhocoupP[ires]*rhos1bw[ires];
      Energy2 Dfact2=-1./3.*(s3-_mpi0*_mpi0-s2+_mpic*_mpic);
      if(ires<_rhocoupD.size()) {
	F2+=Complex(Dfact2*_rhocoupD[ires]*rhos1bw[ires]);
	F3+=Complex(_rhocoupD[ires]*Dfact2*rhos1bw[ires]);
      }
    }
    else if(ichan%2==1&&ichan<=5) {
      unsigned int ires=(ichan-1)/2;
      if(ires<_rhocoupP.size()) F2+=_rhocoupP[ires]*rhos2bw[ires];
      Energy2 Dfact1=-1./3.*(s3-_mpi0*_mpi0-s1+_mpic*_mpic);
      if(ires<_rhocoupD.size()) {
	F1+=Complex(Dfact1*_rhocoupD[ires]*rhos2bw[ires]);
	F3-=Complex(_rhocoupD[ires]*-Dfact1*rhos2bw[ires]);
      }
    }
    else if(ichan==6) {
      F1+=2./3.*_sigmacoup*sigbw;
      F2+=2./3.*_sigmacoup*sigbw;
    }
    else if(ichan==7) {
      Complex Dfact3=1./18./s3*_f2coup*(q2-_mpi0*_mpi0+s3)*(4.*_mpic*_mpic-s3)*f2bw;
      F1+=Dfact3;
      F2+=Dfact3;
      F3-=Complex(0.5*_f2coup*(s1-s2)*f2bw);
    }
    else if(ichan==8) {
      F1+=2./3.*_f0coup*f0bw;
      F2+=2./3.*_f0coup*f0bw;
    }
  }
  else if(imode==5) {
    // identical particle factors
    fact = 1./sqrt(2.);
    // compute the breit wigners we need
    Complex rhos1bw[3],rhos2bw[3];
    for(unsigned int ix=0,N=max(_rhocoupP.size(),_rhocoupD.size());ix<N;++ix) {
      rhos1bw[ix] = Resonance::BreitWignerPWave(s1,_rhomass[ix], _rhowidth[ix],_mpic,_mpic);
      rhos2bw[ix] = Resonance::BreitWignerPWave(s2,_rhomass[ix], _rhowidth[ix],_mpic,_mpic);
    }
    Complex f0bws1  = Resonance::BreitWignerSWave(s1,   _f0mass,   _f0width,_mpic,_mpic);
    Complex sigbws1 = Resonance::BreitWignerSWave(s1,_sigmamass,_sigmawidth,_mpic,_mpic);
    Complex f2bws1  = Resonance::BreitWignerDWave(s1,   _f2mass,   _f2width,_mpic,_mpic);
    Complex f0bws2  = Resonance::BreitWignerSWave(s2,   _f0mass,   _f0width,_mpic,_mpic);
    Complex sigbws2 = Resonance::BreitWignerSWave(s2,_sigmamass,_sigmawidth,_mpic,_mpic);
    Complex f2bws2  = Resonance::BreitWignerDWave(s2,   _f2mass,   _f2width,_mpic,_mpic);
    if(ichan<0) {
      // the p-wave rho terms
      for(unsigned int ix=0;ix<_rhocoupP.size();++ix) {
	F1-=_rhocoupP[ix]*rhos1bw[ix];
	F2-=_rhocoupP[ix]*rhos2bw[ix];
      }
      // the D-wave rho terms
      Energy2 Dfact1=1./3.*(s1-s3);
      Energy2 Dfact2=1./3.*(s2-s3);
      for(unsigned int ix=0;ix<_rhocoupD.size();++ix) {
	F1-=Complex(Dfact1*_rhocoupD[ix]*rhos2bw[ix]);
	F2-=Complex(Dfact2*_rhocoupD[ix]*rhos1bw[ix]);
	F3-=Complex(_rhocoupD[ix]*(Dfact2*rhos1bw[ix]-Dfact1*rhos2bw[ix]));
      }
      // the scalar terms
      F1-=2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
      F2-=2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1);
      F3+=-2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1)
	  +2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
      // the tensor terms
      complex<Energy2> sfact1 = 1./18.*(4.*_mpic*_mpic-s1)*(q2+s1-_mpic*_mpic)/s1*f2bws1;
      complex<Energy2> sfact2 = 1./18.*(4.*_mpic*_mpic-s2)*(q2+s2-_mpic*_mpic)/s2*f2bws2;
      F1+=Complex(_f2coup*(0.5*(s3-s2)*f2bws1-sfact2));
      F2+=Complex(_f2coup*(0.5*(s3-s1)*f2bws2-sfact1));
      F3+=Complex(_f2coup*(-sfact1+sfact2));
    }
    else if(ichan%2==0&&ichan<=4) {
      unsigned int ires=ichan/2;
      Energy2 Dfact2=1./3.*(s2-s3);
      if(ires<_rhocoupP.size()) F1-=_rhocoupP[ires]*rhos1bw[ires];
      if(ires<_rhocoupD.size()) {
	F2-=Complex(Dfact2*_rhocoupD[ires]*rhos1bw[ires]);
	F3-=Complex(_rhocoupD[ires]*Dfact2*rhos1bw[ires]);
      }
    }
    else if(ichan%2==1&&ichan<=5) {
      unsigned int ires=(ichan-1)/2;
      Energy2 Dfact1=1./3.*(s1-s3);
      if(ires<_rhocoupP.size()) {
	F2-=_rhocoupP[ires]*rhos2bw[ires];
      }
      if(ires<_rhocoupD.size()) {
	F1-=Complex(Dfact1*_rhocoupD[ires]*rhos2bw[ires]);
	F3+=Complex(_rhocoupD[ires]*Dfact1*rhos2bw[ires]);
      }
    }
    else if(ichan==6) {
      F2-=2./3.*_sigmacoup*sigbws1;
      F3-=2./3.*_sigmacoup*sigbws1;
    }
    else if(ichan==7) {
      F1-=2./3.*_sigmacoup*sigbws2;
      F3+=2./3.*_sigmacoup*sigbws2;
    }
    else if(ichan==8) {
      complex<Energy2> sfact1 = 1./18.*(4.*_mpic*_mpic-s1)*(q2+s1-_mpic*_mpic)/s1*f2bws1;
      F1+=Complex(_f2coup*0.5*(s3-s2)*f2bws1);
      F2-=Complex(_f2coup*sfact1);
      F3-=Complex(_f2coup*sfact1);
    }
    else if(ichan==9) {
      complex<Energy2> sfact2 = 1./18.*(4.*_mpic*_mpic-s2)*(q2+s2-_mpic*_mpic)/s2*f2bws2;
      F1-=Complex(_f2coup*sfact2);
      F2+=Complex(_f2coup*0.5*(s3-s1)*f2bws2);
      F3+=Complex(_f2coup*sfact2);
    }
    else if(ichan==10) {
      F2-=2./3.*_f0coup*f0bws1;
      F3-=2./3.*_f0coup*f0bws1;
    }
    else if(ichan==11) {
      F1-=2./3.*_f0coup*f0bws2;
      F3+=2./3.*_f0coup*f0bws2;
    }
  }
  else {
    throw Exception() << "ThreePionCLEOCurrent Unknown Decay" << imode
				 << Exception::abortnow;
  }
  F1 *= fact;
  F2 *= fact;
  F3 *= fact;
}

// complete the construction of the decay mode for integration
bool ThreePionCLEOCurrent::createMode(int icharge, tcPDPtr resonance,
				      FlavourInfo flavour,
				      unsigned int imode,PhaseSpaceModePtr mode,
				      unsigned int iloc,int ires,
				      PhaseSpaceChannel phase, Energy upp ) {
  // check the charge and resonance
  if(imode<=1||imode==3||imode==4) {
    if(icharge!=0) return false;
    if(resonance && resonance->id()!=ParticleID::a_10) return false;
  }
  else if(imode==2||imode==5) {
    if(abs(icharge)!=3) return false;
    if(resonance && abs(resonance->id())!=ParticleID::a_1plus) return false;
  }
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==2||imode==5) return false;
      break;
    case IsoSpin::I3One:
      if((imode!=2&&imode!=5) || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if((imode!=2&&imode!=5) || icharge == 3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
  // get the particles and check the masses
  int iq(0),ia(0);
  tPDVector extpart=particles(1,imode,iq,ia);
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  _maxmass=max(_maxmass,upp);
  // pointers to the particles we need
  tPDPtr a1m = getParticleData(ParticleID::a_1minus);
  tPDPtr a10 = getParticleData(ParticleID::a_10);
  // the different rho resonances
  tPDPtr rhom[3] = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) rhom[ix]=rhom[ix]->CC();
    a1m = a1m->CC();
  }
  tPDPtr rho0[3] = {getParticleData(113),getParticleData(100113),getParticleData(30113)};
  // the sigma
  tPDPtr sigma = getParticleData(9000221);
  // the f_2
  tPDPtr f2=getParticleData(225);
  // the f_0
  tPDPtr f0=getParticleData(10221);
  assert(f2 && f0 && sigma);
  // a0 -> pi0 pi0 pi0
  if(imode<=1) {
    for(unsigned int ix=0;ix<3;++ix) {
      tPDPtr temp;
      if(ix==0)      temp = sigma;
      else if(ix==1) temp = f2;
      else if(ix==2) temp = f0;
      mode->addChannel((PhaseSpaceChannel(phase),ires,a10,ires+1,temp,ires+1,iloc+1,
			ires+2,iloc+2,ires+2,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a10,ires+1,temp,ires+1,iloc+2,
			ires+2,iloc+1,ires+2,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a10,ires+1,temp,ires+1,iloc+3,
			ires+2,iloc+1,ires+2,iloc+2));
    }
  }
  // decay mode a_1- -> pi0 pi0 pi-
  else if(imode==2) {
    for(unsigned int ix=0;ix<3;++ix) {
      // first rho+ channel
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,ires+1,rhom[ix],ires+1,iloc+1,
			ires+2,iloc+2,ires+2,iloc+3));
      // second rho+ channel
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,ires+1,rhom[ix],ires+1,iloc+2,
			ires+2,iloc+1,ires+2,iloc+3));
    }
    // I=0 channels
    for(unsigned int iy=0;iy<3;++iy) {
      tPDPtr temp;
      if(iy==0)      temp = sigma;
      else if(iy==1) temp = f2;
      else if(iy==2) temp = f0;
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,ires+1,temp,ires+1,iloc+3,
			ires+2,iloc+1,ires+2,iloc+2));
    }
  }
  // decay mode a_10 -> pi+ pi- pi0
  else if(imode==3||imode==4) {
    // rho modes
    for(unsigned int ix=0;ix<3;++ix) {
      // first rho channel
      mode->addChannel((PhaseSpaceChannel(phase),ires,a10,ires+1,rhom[ix],ires+1,iloc+1,
			ires+2,iloc+2,ires+2,iloc+3));
      // second channel
      mode->addChannel((PhaseSpaceChannel(phase),ires,a10,ires+1,rhom[ix],ires+1,iloc+2,
			ires+2,iloc+1,ires+2,iloc+3));
    }
    // I=0 channels
    for(unsigned int iy=0;iy<3;++iy) {
      tPDPtr temp;
      if(iy==0)      temp = sigma;
      else if(iy==1) temp = f2;
      else if(iy==2) temp = f0;
      mode->addChannel((PhaseSpaceChannel(phase),ires,a10,ires+1,temp,ires+1,iloc+3,
			ires+2,iloc+1,ires+2,iloc+2));
    }
  }
  else if(imode==5) {
    for(unsigned int ix=0;ix<3;++ix) {
      // the neutral rho channels
      // first channel
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,ires+1,rho0[ix],ires+1,iloc+1,
			ires+2,iloc+2,ires+2,iloc+3));
      // interchanged channel
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,ires+1,rho0[ix],ires+1,iloc+2,
			ires+2,iloc+1,ires+2,iloc+3));
    }
    for(unsigned int iy=0;iy<3;++iy) {
      tPDPtr temp;
      if(iy==0)      temp = sigma;
      else if(iy==1) temp = f2;
      else if(iy==2) temp = f0;
      // first channel
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,ires+1,temp,ires+1,iloc+1,
			ires+2,iloc+2,ires+2,iloc+3));
      // interchanged channel
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,ires+1,temp,ires+1,iloc+2,
			ires+2,iloc+1,ires+2,iloc+3));
    }
  }
  // reset the integration parameters
  for(unsigned int iy=0;iy<_rhomass.size();++iy) {
    mode->resetIntermediate(rho0[iy],_rhomass[iy],_rhowidth[iy]);
    mode->resetIntermediate(rhom[iy],_rhomass[iy],_rhowidth[iy]);
  }
  mode->resetIntermediate(sigma,_sigmamass,_sigmawidth);
  mode->resetIntermediate(f2,_f2mass,_f2width);
  mode->resetIntermediate(f0,_f0mass,_f0width);
  mode->resetIntermediate(a10,_a1mass,_a1width);
  mode->resetIntermediate(a10,_a1mass,_a1width);
  return true;
}

void ThreePionCLEOCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const {
  if(header){output << "update decayers set parameters=\"";}
  if(create) {
    output << "create Herwig::ThreePionCLEOCurrent " << name()
	   << " HwWeakCurrents.so\n";
  }
  for(unsigned int ix=0;ix<_rhomass.size();++ix) {
    if(ix<2) {
      output << "newdef    " << name() << ":RhoMasses " << ix
	     << " " << _rhomass[ix]/MeV << "\n";
    }
    else {
      output << "insert " << name() << ":RhoMasses " << ix
	     << " " << _rhomass[ix]/MeV << "\n";
    }
  }
  for(unsigned int ix=0;ix<_rhowidth.size();++ix) {
    if(ix<2) {
      output << "newdef    " << name() << ":RhoWidths " << ix
	     << " " << _rhowidth[ix]/MeV << "\n";
    }
    else {
      output << "insert " << name() << ":RhoWidths " << ix
	     << " " << _rhowidth[ix]/MeV << "\n";
    }
  }
  output << "newdef " << name() << ":f_2Mass " << _f2mass/GeV << "\n";
  output << "newdef " << name() << ":f_2Width " << _f2width/GeV << "\n";
  output << "newdef " << name() << ":f_0Mass " << _f0mass/GeV << "\n";
  output << "newdef " << name() << ":f_0Width " << _f0width/GeV << "\n";
  output << "newdef " << name() << ":sigmaMass " << _sigmamass/GeV << "\n";
  output << "newdef " << name() << ":sigmaWidth " << _sigmawidth/GeV << "\n";
  output << "newdef " << name() << ":a1Mass " << _a1mass/GeV << "\n";
  output << "newdef " << name() << ":a1Width " <<_a1width /GeV << "\n";
  output << "newdef " << name() << ":KaonMass " << _mK/GeV << "\n";
  output << "newdef " << name() << ":KStarMass " << _mKstar/GeV << "\n";
  output << "newdef " << name() << ":KaonCoupling " << _gammk << "\n";
  output << "newdef " << name() << ":Fpi " << _fpi/MeV << "\n";
  output << "newdef " << name() << ":a1WidthOption " << _a1opt << "\n";
  for(unsigned int ix=0;ix<_rhomagP.size();++ix) {
      if(ix<2) {
	output << "newdef    " << name() << ":RhoPWaveMagnitude " << ix
	       << " " << _rhomagP[ix] << "\n";
      }
      else {
	output << "insert " << name() << ":RhoPWaveMagnitude " << ix
	       << " " << _rhomagP[ix] << "\n";
      }
  }
  for(unsigned int ix=0;ix<_rhophaseP.size();++ix) {
    if(ix<2) {
      output << "newdef    " << name() << ":RhoPWavePhase " << ix
	     << " " << _rhophaseP[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":RhoPWavePhase " << ix
	     << " " << _rhophaseP[ix] << "\n";
    }
  }
  for(unsigned int ix=0;ix<_rhomagD.size();++ix) {
    if(ix<2) {
      output << "newdef    " << name() << ":RhoDWaveMagnitude " << ix
	     << " " << _rhomagD[ix]*MeV2 << "\n";
    }
    else {
      output << "insert " << name() << ":RhoDWaveMagnitude " << ix
	     << " " << _rhomagD[ix]*MeV2 << "\n";
    }
  }
  for(unsigned int ix=0;ix<_rhophaseD.size();++ix) {
    if(ix<2) {
      output << "newdef    " << name() << ":RhoDWavePhase " << ix
	     << " " << _rhophaseD[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":RhoDWavePhase " << ix
	     << " " << _rhophaseD[ix] << "\n";
    }
  }
  output << "newdef " << name() << ":f0Phase " << _f0phase << "\n";
  output << "newdef " << name() << ":f2Phase " <<_f2phase  << "\n";
  output << "newdef " << name() << ":sigmaPhase " <<_sigmaphase  << "\n";
  output << "newdef " << name() << ":f0Magnitude " << _f0mag << "\n";
  output << "newdef " << name() << ":f2Magnitude " << _f2mag*GeV2 << "\n";
  output << "newdef " << name() << ":sigmaMagnitude " <<_sigmamag  << "\n";
  output << "newdef " << name() << ":Initializea1 " <<_initializea1  << "\n";
  for(unsigned int ix=0;ix<_a1runwidth.size();++ix) {
    if(ix<200) {
      output << "newdef    " << name() << ":a1RunningWidth " << ix
	     << " " << _a1runwidth[ix]/MeV << "\n";
    }
    else {
      output << "insert " << name() << ":a1RunningWidth " << ix
	     << " " << _a1runwidth[ix]/MeV << "\n";
    }
  }
  for(unsigned int ix=0;ix<_a1runq2.size();++ix) {
    if(ix<200) {
      output << "newdef    " << name() << ":a1RunningQ2 " << ix
	     << " " << _a1runq2[ix]/MeV2 << "\n";
    }
    else {
      output << "insert " << name() << ":a1RunningQ2 " << ix
	     << " " << _a1runq2[ix]/MeV2 << "\n";
    }
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\""
		    << fullName() << "\";" << endl;
}

void ThreePionCLEOCurrent::doinitrun() {
  // set up the running a_1 width
  inita1Width(0);
  WeakCurrent::doinitrun();
}

void ThreePionCLEOCurrent::doupdate() {
  WeakCurrent::doupdate();
  // update running width if needed
  if ( !touched() ) return;
  if(_maxmass!=_maxcalc) inita1Width(-1);
}

Energy ThreePionCLEOCurrent::a1width(Energy2 q2) const {
  Energy output;
  if(_a1opt) output=(*_a1runinter)(q2);
  else {
    double gam(0.);
    if(q2<0.1753*GeV2) {
      gam =0.;
    }
    else if(q2<0.823*GeV2) {
      double p=q2/GeV2-0.1753;
      gam = 5.80900*p*sqr(p)*(1.-3.00980*p+4.57920*sqr(p));
    }
    else {
      double p=q2/GeV2;
      gam = -13.91400+27.67900*p-13.39300*sqr(p)
	+3.19240*sqr(p)*p-0.10487*sqr(sqr(p));
    }
    if(q2<0.1676*GeV2) {
      gam+=0.;
    }
    else if(q2<0.823*GeV2) {
      double p=q2/GeV2-0.1676;
      gam+= 6.28450*p*sqr(p)*(1.-2.95950*p+4.33550*sqr(p));
    }
    else {
      double p=q2/GeV2;
      gam+= -15.41100+32.08800*p-17.66600*sqr(p)
	+4.93550*sqr(p)*p-0.37498*sqr(sqr(p));
    }
    Energy mkst=0.894*GeV,mk=0.496*GeV;
    Energy2 mk1sq=sqr(mkst+mk), mk2sq=sqr(mkst-mk);
    double c3pi=sqr(0.2384),ckst=sqr(4.7621)*c3pi;
    gam*=c3pi;
    if(q2>mk1sq) gam+=0.5*ckst*sqrt((q2-mk1sq)*(q2-mk2sq))/q2;
    gam = gam*_a1width*_a1mass/GeV2/1.331/0.814/1.0252088;
    output = gam*GeV2/sqrt(q2);
  }
  return output;
}

double
ThreePionCLEOCurrent::threeBodyMatrixElement(const int iopt, const Energy2 q2,
					     const Energy2 s3, const Energy2 s2,
					     const Energy2 s1, const Energy,
					     const Energy, const Energy) const {
  Energy p1[5],p2[5],p3[5];
  Energy2 p1sq, p2sq, p3sq;
  Energy q=sqrt(q2);
  Energy2 mpi2c=_mpic*_mpic;
  Energy2 mpi20=_mpi0*_mpi0;
  // construct the momenta for the 2 neutral 1 charged mode
  Complex F1,F2,F3;
  if(iopt==0) {
    // construct the momenta of the decay products
    p1[0] = 0.5*(q2+mpi20-s1)/q; p1sq=p1[0]*p1[0]; p1[4]=sqrt(p1sq-mpi20);
    p2[0] = 0.5*(q2+mpi20-s2)/q; p2sq=p2[0]*p2[0]; p2[4]=sqrt(p2sq-mpi20);
    p3[0] = 0.5*(q2+mpi2c-s3)/q; p3sq=p3[0]*p3[0]; p3[4]=sqrt(p3sq-mpi2c);
    // take momentum of 1 parallel to z axis
    p1[1]=ZERO;p1[2]=ZERO;p1[3]=p1[4];
    // construct 2
    double cos2 = 0.5*(p1sq+p2sq-p3sq-2.*mpi20+mpi2c)/p1[4]/p2[4];
    p2[1] = p2[4]*sqrt(1.-cos2*cos2); p2[2]=ZERO; p2[3]=-p2[4]*cos2;
    // construct 3
    double cos3 = 0.5*(p1sq-p2sq+p3sq-mpi2c)/p1[4]/p3[4];
    p3[1] =-p3[4]*sqrt(1.-cos3*cos3); p3[2]=ZERO; p3[3]=-p3[4]*cos3;
    // calculate the form factors
    CLEOFormFactor(1,-1,q2,s1,s2,s3,F1,F2,F3);
  }
  // construct the momenta for the 3 charged mode
  else {
    // construct the momenta of the decay products
    p1[0] = 0.5*(q2+mpi2c-s1)/q; p1sq=p1[0]*p1[0]; p1[4]=sqrt(p1sq-mpi2c);
    p2[0] = 0.5*(q2+mpi2c-s2)/q; p2sq=p2[0]*p2[0]; p2[4]=sqrt(p2sq-mpi2c);
    p3[0] = 0.5*(q2+mpi2c-s3)/q; p3sq=p3[0]*p3[0]; p3[4]=sqrt(p3sq-mpi2c);
    // take momentum of 1 parallel to z axis
    p1[1]=ZERO;p1[2]=ZERO;p1[3]=p1[4];
    // construct 2
    double cos2 = 0.5*(p1sq+p2sq-p3sq-mpi2c)/p1[4]/p2[4];
    p2[1] = p2[4]*sqrt(1.-cos2*cos2); p2[2]=ZERO; p2[3]=-p2[4]*cos2;
    // construct 3
    double cos3 = 0.5*(p1sq-p2sq+p3sq-mpi2c)/p1[4]/p3[4];
    p3[1] =-p3[4]*sqrt(1.-cos3*cos3); p3[2]=ZERO; p3[3]=-p3[4]*cos3;
    // calculate the form factors
    CLEOFormFactor(0,-1,q2,s1,s2,s3,F1,F2,F3);
  }
  // construct a vector with the current
  complex<Energy> current[4];
  for(unsigned int ix=0;ix<4;++ix)
    current[ix] = F1*(p2[ix]-p3[ix])-F2*(p3[ix]-p1[ix])+F3*(p1[ix]-p2[ix]);
  complex<Energy2> dot1=current[0]*conj(current[0]);
  for(unsigned int ix=1;ix<4;++ix) dot1-=current[ix]*conj(current[ix]);
  complex<Energy2> dot2=current[0]*q;
  return(-dot1+dot2*conj(dot2)/q2).real() / sqr(_rhomass[0]);
}

// the hadronic currents
vector<LorentzPolarizationVectorE>
ThreePionCLEOCurrent::current(tcPDPtr resonance,
			      FlavourInfo flavour,
			      const int imode, const int ichan, Energy & scale,
			      const tPDVector & ,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  useMe();
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IOne)
    return vector<LorentzPolarizationVectorE>();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==2||imode==5) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode!=2&&imode!=5) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode!=2&&imode!=5) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return vector<LorentzPolarizationVectorE>();
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q;
  for(unsigned int ix=0;ix<momenta.size();++ix)
    q+=momenta[ix];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 s1 = (momenta[1]+momenta[2]).m2();
  Energy2 s2 = (momenta[0]+momenta[2]).m2();
  Energy2 s3 = (momenta[0]+momenta[1]).m2();
  // form factors
  Complex F1(0.), F2(0.), F3(0.);
  CLEOFormFactor(imode,ichan,q2,s1,s2,s3,F1,F2,F3);
  // change sign of the F2 term
  F2 =- F2;
  // prefactor
  complex<InvEnergy> a1fact = _fact;
  if(!resonance) a1fact = a1fact * a1BreitWigner(q2);
  // current
  LorentzPolarizationVectorE vect = q.mass()*a1fact*
    ((F2-F1)*momenta[2] + (F1-F3)*momenta[1] + (F3-F2)*momenta[0]);
  // scalar piece
  Complex dot=(vect*q)/q2;
  vect -= dot*q;
  // return the answer
  return vector<LorentzPolarizationVectorE>(1,vect);
}

bool ThreePionCLEOCurrent::accept(vector<int> id) {
  if(id.size()!=3) return false;
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus) continue;
    else if(id[ix]==ParticleID::piminus) continue;
    else if(id[ix]==ParticleID::pi0)     continue;
    return false;
  }
  return true;
}

unsigned int ThreePionCLEOCurrent::decayMode(vector<int> id) {
  if(id.size()!=3) return -1;
  int npip(0),npim(0),npi0(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
  }
  if       (npi0==3)                                  return 0;
  else if( (npip==1&&npi0==2) || (npim==1&&npi0==2) ) return 2;
  else if( npi0==1 && npip==1 && npim==1 )            return 3;
  else if( (npip==2&&npim==1) || (npim==2&&npip==1) ) return 5;
  else return -1;
}

tPDVector ThreePionCLEOCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart(3);
  if(imode==0||imode==1) {
    extpart[0]=getParticleData(ParticleID::pi0);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::pi0);
  }
  else if(imode==2) {
    extpart[0]=getParticleData(ParticleID::pi0);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::piminus);
  }
  else if(imode==3||imode==4) {
    extpart[0]=getParticleData(ParticleID::piplus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::pi0);
  }
  else if(imode==5) {
    extpart[0]=getParticleData(ParticleID::piminus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::piplus);
  }
  else
    assert(false);
  // conjugate the particles if needed
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(extpart[ix]->CC()) extpart[ix]=extpart[ix]->CC();
    }
  }
  // return the answer
  return extpart;
}
#line 1 "./TwoPionRhoCurrent.cc"
// -*- C++ -*-
//
// TwoPionRhoCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoPionRhoCurrent class.
//
//  Author: Peter Richardson
//

#include "TwoPionRhoCurrent.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG::Helicity;

TwoPionRhoCurrent::TwoPionRhoCurrent() {
  // set up for the modes in the base class
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(4);
  // the weights of the different resonances in the matrix elements
  _pimag   = {1.0,0.167,0.05};
  _piphase = {0.0,180  ,0.0};
  // models to use
  _pimodel = 0;
  // parameter for the masses (use the parameters freom the CLEO fit 
  // rather than the PDG masses etc)
  _rhoparameters=true;
  _rhomasses = {774.6*MeV,1408*MeV,1700*MeV};
  _rhowidths = {149*MeV,502*MeV,235*MeV};
}

void TwoPionRhoCurrent::doinit() {
  WeakCurrent::doinit();
  // check consistency of parametrers
  if(_rhomasses.size()!=_rhowidths.size()) {
    throw InitException() << "Inconsistent parameters in TwoPionRhoCurrent"
			  << "::doinit()" << Exception::abortnow;
  }
  // the resonances
  tPDPtr res[3]={getParticleData(-213   ),
		 getParticleData(-100213),
		 getParticleData(-30213 )};
  // reset the masses in the form-factors if needed
  if(_rhoparameters&&_rhomasses.size()<3) {
    for(unsigned int ix=_rhomasses.size();ix<3;++ix) {
      if(res[ix]) _rhomasses.push_back(res[ix]->mass() );
      if(res[ix]) _rhowidths.push_back(res[ix]->width());
    }
  }
  else if(!_rhoparameters) {
    _rhomasses.clear();_rhowidths.clear();
    for(unsigned int ix=0;ix<3;++ix) {
      if(res[ix]) _rhomasses.push_back(res[ix]->mass() );
      if(res[ix]) _rhowidths.push_back(res[ix]->width());
    }
  }
  // set up for the Breit Wigners
  Energy mpi0(   getParticleData(ParticleID::pi0   )->mass());
  Energy mpiplus(getParticleData(ParticleID::piplus)->mass());
  // rho resonances
  for(unsigned int ix=0;ix<3;++ix) {
    _mass.push_back(_rhomasses[ix]);
    _width.push_back(_rhowidths[ix]);
    _massa.push_back(mpi0);
    _massb.push_back(mpiplus);
    _hres.push_back(Resonance::Hhat(sqr(_mass.back()),_mass.back(),_width.back(),_massa.back(),_massb.back()));
    _dh.push_back(Resonance::dHhatds(_mass.back(),_width.back(),_massa.back(),_massb.back()));
    _h0.push_back(Resonance::H(ZERO,_mass.back(),_width.back(),_massa.back(),_massb.back(),_dh.back(),_hres.back()));
  }
  // weights for the rho channels
  if(_pimag.size()!=_piphase.size()) 
    throw InitException() << "The vectors containing the weights and phase for the"
			  << " rho channel must be the same size in "
			  << "TwoPionRhoCurrent::doinit()" << Exception::runerror;
  _piwgt.resize(_pimag.size());
  for(unsigned int ix=0;ix<_pimag.size();++ix) {
    double angle = _piphase[ix]/180.*Constants::pi;
    _piwgt[ix] = _pimag[ix]*(cos(angle)+Complex(0.,1.)*sin(angle));
  }
}

void TwoPionRhoCurrent::persistentOutput(PersistentOStream & os) const {
  os << _pimodel << _piwgt << _pimag << _piphase
     << _rhoparameters << ounit(_rhomasses,GeV) << ounit(_rhowidths,GeV)
     << ounit(_mass,GeV) << ounit(_width,GeV)
     << ounit(_massa,GeV) <<ounit(_massb,GeV)
     << _dh << ounit(_hres,GeV2) << ounit(_h0,GeV2);
}

void TwoPionRhoCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _pimodel >> _piwgt >> _pimag >> _piphase
     >> _rhoparameters >> iunit(_rhomasses,GeV) >> iunit(_rhowidths,GeV) 
     >> iunit(_mass,GeV) >> iunit(_width,GeV)
     >> iunit(_massa,GeV) >> iunit(_massb,GeV)
     >> _dh >> iunit(_hres,GeV2) >> iunit(_h0,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TwoPionRhoCurrent,WeakCurrent>
describeHerwigTwoPionRhoCurrent("Herwig::TwoPionRhoCurrent", "HwWeakCurrents.so");

void TwoPionRhoCurrent::Init() {

  static ParVector<TwoPionRhoCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the pi pi channel",
     &TwoPionRhoCurrent::_rhomasses, MeV, -1, 775.8*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<TwoPionRhoCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the pi pi channel",
     &TwoPionRhoCurrent::_rhowidths, MeV, -1, 150.3*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static Switch<TwoPionRhoCurrent,bool> interfaceRhoParameters
    ("RhoParameters",
     "Use local values for the rho meson masses and widths",
     &TwoPionRhoCurrent::_rhoparameters, true, false, false);
  static SwitchOption interfaceRhoParameterstrue
    (interfaceRhoParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceRhoParametersParticleData
    (interfaceRhoParameters,
     "ParticleData",
     "Use the value from the particle data objects",
     false);
  
  static ParVector<TwoPionRhoCurrent,double> interfacePiMagnitude
    ("PiMagnitude",
     "Magnitude of the weight of the different resonances for the pi pi channel",
     &TwoPionRhoCurrent::_pimag, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoPionRhoCurrent,double> interfacePiPhase
    ("PiPhase",
     "Phase of the weight of the different resonances for the pi pi channel",
     &TwoPionRhoCurrent::_piphase, -1, 0., 0, 0,
     false, false, Interface::nolimits);
  
  static Switch<TwoPionRhoCurrent,int> interfacePiModel
    ("PiModel",
     "The model to use for the propagator for the pion modes.",
     &TwoPionRhoCurrent::_pimodel, 0, false, false);
  static SwitchOption interfacePiModelKuhn
    (interfacePiModel,
     "Kuhn",
     "The model of Kuhn and Santamaria",
     0);
  static SwitchOption interfacePiModelGounaris
    (interfacePiModel,
     "Gounaris",
     "The model of Gounaris and Sakurai.",
     1);

  static ClassDocumentation<TwoPionRhoCurrent> documentation
    ("The TwoPionRhoCurrent class is designed to implement weak"
     "decay to two scalar mesons using the models of either Kuhn and "
     "Santamaria (Z. Phys. C48, 445 (1990)) or Gounaris and Sakurai Phys. Rev. "
     "Lett. 21, 244 (1968).  The mixing parameters are taken from "
     "Phys. Rev. D61:112002,2000 (CLEO), although the PDG values for the "
     "masses and widths are used, for the decay pi+/- pi0.",
     "The weak "
     "decay current to two scalar mesons is implemented "
     "using the models of either Kuhn and "
     "Santamaria \\cite{Kuhn:1990ad} or Gounaris and Sakurai \\cite{Gounaris:1968mw}. "
     "The mixing parameters are taken from "
     "\\cite{Asner:1999kj}, although the PDG values for the "
     "masses and widths are used, for the decay pi+/- pi0.",
     "%\\cite{Kuhn:1990ad}\n"
     "\\bibitem{Kuhn:1990ad}\n"
     "  J.~H.~Kuhn and A.~Santamaria,\n"
     "  %``Tau decays to pions,''\n"
     "  Z.\\ Phys.\\  C {\\bf 48}, 445 (1990).\n"
     "  %%CITATION = ZEPYA,C48,445;%%\n"
     "%\\cite{Gounaris:1968mw}\n"
     "\\bibitem{Gounaris:1968mw}\n"
     "  G.~J.~Gounaris and J.~J.~Sakurai,\n"
     "   ``Finite width corrections to the vector meson dominance prediction for rho\n"
     "  %$\\to$ e+ e-,''\n"
     "  Phys.\\ Rev.\\ Lett.\\  {\\bf 21}, 244 (1968).\n"
     "  %%CITATION = PRLTA,21,244;%%\n"
     "%\\cite{Asner:1999kj}\n"
     "\\bibitem{Asner:1999kj}\n"
     "  D.~M.~Asner {\\it et al.}  [CLEO Collaboration],\n"
     "   ``Hadronic structure in the decay tau- --> nu/tau pi- pi0 pi0 and the  sign\n"
     "  %of the tau neutrino helicity,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 61}, 012002 (2000)\n"
     "  [arXiv:hep-ex/9902022].\n"
     "  %%CITATION = PHRVA,D61,012002;%%\n"
     );

}

// complete the construction of the decay mode for integration
bool TwoPionRhoCurrent::createMode(int icharge, tcPDPtr resonance,
				    FlavourInfo flavour,
				    unsigned int imode,PhaseSpaceModePtr mode,
				    unsigned int iloc,int ires,
				    PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((imode<1 && abs(icharge)!=3) ||
     (imode>1  && icharge !=0)) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode<=1) return false;
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return false;
  // make sure that the decays are kinematically allowed
  tPDPtr part[2];
  if(imode==0) {
    part[0]=getParticleData(ParticleID::piplus);
    part[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==1) {
    part[0]=getParticleData(ParticleID::Kplus);
    part[1]=getParticleData(ParticleID::Kbar0);
  }
  else if(imode==2 || imode==3 ) {
    part[0]=getParticleData(ParticleID::piplus);
    part[1]=getParticleData(ParticleID::piminus);
  }
  Energy min(part[0]->massMin()+part[1]->massMin());
  if(min>upp) return false;
  // set the resonances
  // two pion or  K+ K0 decay
  tPDPtr res[3]={getParticleData(213),getParticleData(100213),getParticleData(30213)};
  if(icharge==-3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(res[ix]&&res[ix]->CC()) res[ix]=res[ix]->CC();
    }
  }
  // create the channels
  for(unsigned int ix=0;ix<3;++ix) {
    if(!res[ix]) continue;
    if(resonance && resonance != res[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],ires+1,iloc+1,ires+1,iloc+2));
  }
  // reset the masses in the intergrators if needed
  // for the rho 
  if(_rhoparameters) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(ix<_rhomasses.size()&&res[ix]) {
  	mode->resetIntermediate(res[ix],_rhomasses[ix],_rhowidths[ix]);
      }
    }
  }
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector TwoPionRhoCurrent::particles(int icharge, unsigned int imode,
				       int,int) {
  tPDVector output(2);
  if(imode==0) {
    output[0]=getParticleData(ParticleID::piplus);
    output[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==1) {
    output[0]=getParticleData(ParticleID::Kplus);
    output[1]=getParticleData(ParticleID::Kbar0);
  }
  else {
    output[0]=getParticleData(ParticleID::piplus);
    output[1]=getParticleData(ParticleID::piminus);
  }
  if(icharge==-3) {
    for(unsigned int ix=0;ix<output.size();++ix) {
      if(output[ix]->CC()) output[ix]=output[ix]->CC();
    }
  }
  return output;
}

// hadronic current   
vector<LorentzPolarizationVectorE> 
TwoPionRhoCurrent::current(tcPDPtr resonance,
			   FlavourInfo flavour,
			   const int imode, const int ichan,Energy & scale, 
			   const tPDVector & outgoing,
			   const vector<Lorentz5Momentum> & momenta,
			   DecayIntegrator::MEOption) const {
  useMe();
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IOne)
    return vector<LorentzPolarizationVectorE>();
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode<=1) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return vector<LorentzPolarizationVectorE>();
  // momentum difference and sum of the mesons
  Lorentz5Momentum pdiff(momenta[0]-momenta[1]);
  Lorentz5Momentum psum (momenta[0]+momenta[1]);
  psum.rescaleMass();
  scale=psum.mass();
  // mass2 of vector intermediate state
  Energy2 q2(psum.m2());
  double dot(psum*pdiff/q2);
  psum *=dot;
  // calculate the current
  Complex FPI(0.);
  Complex denom = std::accumulate(_piwgt.begin(),_piwgt.end(),Complex(0.));
  unsigned int imin=0, imax = _piwgt.size();
  if(ichan>0) {
    imin = ichan;
    imax = ichan+1;
  }
  if(resonance) {
    switch(resonance->id()/1000) {
    case 0:
      imin = 0;
      break;
    case 100:
      imin = 1;
      break;
    case 30 :
      imin = 2;
      break;
    default:
      assert(false);
    }
    imax=imin+1;
  }
  // rho
  for(unsigned int ix=imin;ix<imax;++ix) {
    FPI+=_piwgt[ix]*BreitWigner(q2,_pimodel,ix);
  }
  // additional prefactors
  FPI/=denom;
  // pion mode
  if(imode==0)           FPI *= sqrt(2.0);
    // two kaon modes
  else if(imode==1)      FPI *= 1.       ;
  // compute the current
  pdiff-=psum;
  return vector<LorentzPolarizationVectorE>(1,FPI*pdiff);
}
   
bool TwoPionRhoCurrent::accept(vector<int> id) {
  // check there are only two particles
  if(id.size()!=2) return false;
  // pion modes
  if((abs(id[0])==ParticleID::piplus  &&     id[1] ==ParticleID::pi0   ) ||
     (    id[0] ==ParticleID::pi0     && abs(id[1])==ParticleID::piplus))
    return true;
  // two kaons
  else if((id[0]==ParticleID::Kminus && id[1]==ParticleID::K0)     ||
	  (id[0]==ParticleID::K0     && id[1]==ParticleID::Kminus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::Kbar0)  ||
	  (id[0]==ParticleID::Kbar0  && id[1]==ParticleID::Kplus))
    return true;
  else if((id[0]==ParticleID::piminus && id[1]==ParticleID::piplus) ||
	  (id[0]==ParticleID::piplus  && id[1]==ParticleID::piminus))
    return true;
  else
    return false;
}

// the decay mode
unsigned int TwoPionRhoCurrent::decayMode(vector<int> idout) {
  unsigned int nkaon(0),npi(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::K0) {
      ++nkaon;
    }
    else if (abs(idout[ix])==ParticleID::Kplus) {
      ++nkaon;
    }
    else if(abs(idout[ix])==ParticleID::piplus) {
      ++npi;
    }
  }
  if(nkaon==2)    return 1;
  else if(npi==2) return 2;
  else            return 0;
}

// output the information for the database
void TwoPionRhoCurrent::dataBaseOutput(ofstream & output,bool header,
					     bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::TwoPionRhoCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  unsigned int ix;
  for(ix=0;ix<_rhomasses.size();++ix) {
    if(ix<3)  output << "newdef ";
    else      output << "insert ";
    output << name() << ":RhoMasses " << ix << " " << _rhomasses[ix]/MeV << "\n";
  }
  for(ix=0;ix<_rhowidths.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidths " << ix << " " << _rhowidths[ix]/MeV << "\n";
  }
  output << "newdef " << name() << ":RhoParameters " << _rhoparameters << "\n";
  for(ix=0;ix<_piwgt.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":PiMagnitude " << ix << " " << _pimag[ix]   << "\n";
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":PiPhase "     << ix << " " << _piphase[ix] << "\n";
  }
  output << "newdef " << name() << ":PiModel " << _pimodel << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./KPiKStarCurrent.cc"
// -*- C++ -*-
//
// KPiKStarCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KPiKStarCurrent class.
//
//  Author: Peter Richardson
//

#include "KPiKStarCurrent.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

KPiKStarCurrent::KPiKStarCurrent() {
  // set up for the modes in the base class
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  setInitialModes(3);
  // the weights of the different resonances in the matrix elements
  _kmag    = {1.0,0.038,0.0};
  _kphase  = {0.0,180  ,0.0};
  // model to use
  _kmodel  = 0;
  // parameter for the masses (use the parameters freom the CLEO fit 
  // rather than the PDG masses etc)
  _kstarparameters=true;
  _kstarmasses = {0.8921*GeV,1.700*GeV};
  _kstarwidths = {0.0513*GeV,0.235*GeV};
}

void KPiKStarCurrent::doinit() {
  WeakCurrent::doinit();
  // check consistency of parametrers
  if(_kstarmasses.size()!=_kstarwidths.size()) {
    throw InitException() << "Inconsistent parameters in KPiKStarCurrent"
			  << "::doinit()" << Exception::abortnow;
  }
  // the resonances
  tPDPtr res[3]={getParticleData(-323   ),
		 getParticleData(-100323),
		 getParticleData(-30323 )};
  // reset the masses in the form-factors if needed
  if(_kstarparameters&&_kstarmasses.size()<3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(res[ix]) _kstarmasses.push_back(res[ix]->mass());
      if(res[ix]) _kstarwidths.push_back(res[ix]->width());
    }
  }
  else if(!_kstarparameters) {
    _kstarmasses.clear();_kstarwidths.clear();
    for(unsigned int ix=0;ix<3;++ix) {
      if(res[ix]) _kstarmasses.push_back(res[ix]->mass());
      if(res[ix]) _kstarwidths.push_back(res[ix]->width());
    }
  }
  // set up for the Breit Wigners
  Energy mpiplus(getParticleData(ParticleID::piplus)->mass());
  Energy mk0(    getParticleData(ParticleID::K0    )->mass());
  // Kstar resonances
  for(unsigned int ix=0;ix<3;++ix) {
    _mass.push_back(_kstarmasses[ix]);
    _width.push_back(_kstarwidths[ix]);
    _massa.push_back(mk0);
    _massb.push_back(mpiplus);
    _hres.push_back(Resonance::Hhat(sqr(_mass.back()),_mass.back(),_width.back(),_massa.back(),_massb.back()));
    _dh.push_back(Resonance::dHhatds(_mass.back(),_width.back(),_massa.back(),_massb.back()));
    _h0.push_back(Resonance::H(ZERO,_mass.back(),_width.back(),_massa.back(),_massb.back(),_dh.back(),_hres.back()));
  }
  // weights for the K* channels
  if(_kmag.size()!=_kphase.size()) 
    throw InitException() << "The vectors containing the weights and phase for the"
			  << " K* channel must be the same size in "
			  << "KPiKStarCurrent::doinit()" << Exception::runerror;
  _kwgt.resize(_kmag.size());
  for(unsigned int ix=0;ix<_kmag.size();++ix) {
    double angle = _kphase[ix]/180.*Constants::pi;
    _kwgt[ix] = _kmag[ix]*(cos(angle)+Complex(0.,1.)*sin(angle));
  }
}

void KPiKStarCurrent::persistentOutput(PersistentOStream & os) const {
  os << _kmodel << _kwgt << _kmag 
     << _kphase << _kstarparameters
     << ounit(_kstarmasses,GeV) << ounit(_kstarwidths,GeV) 
     << ounit(_mass,GeV) << ounit(_width,GeV)
     << ounit(_massa,GeV) <<ounit(_massb,GeV)
     << _dh << ounit(_hres,GeV2) << ounit(_h0,GeV2);
}

void KPiKStarCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _kmodel >> _kwgt >> _kmag 
     >> _kphase >> _kstarparameters
     >> iunit(_kstarmasses,GeV) >> iunit(_kstarwidths,GeV) 
     >> iunit(_mass,GeV) >> iunit(_width,GeV)
     >> iunit(_massa,GeV) >> iunit(_massb,GeV)
     >> _dh >> iunit(_hres,GeV2) >> iunit(_h0,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KPiKStarCurrent,WeakCurrent>
describeHerwigKPiKStarCurrent("Herwig::KPiKStarCurrent", "HwWeakCurrents.so");

void KPiKStarCurrent::Init() {
  
  static ParVector<KPiKStarCurrent,Energy> interfaceKstarMasses
    ("KstarMasses",
     "The masses of the different K* resonances for the pi pi channel",
     &KPiKStarCurrent::_kstarmasses, MeV, -1, 891.66*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<KPiKStarCurrent,Energy> interfaceKstarWidths
    ("KstarWidths",
     "The widths of the different K* resonances for the pi pi channel",
     &KPiKStarCurrent::_kstarwidths, MeV, -1, 50.8*MeV, ZERO, 1000.*MeV,
     false, false, true);

  static Switch<KPiKStarCurrent,bool> interfaceKstarParameters
    ("KstarParameters",
     "Use local values for the Kstar meson masses and widths",
     &KPiKStarCurrent::_kstarparameters, true, false, false);
  static SwitchOption interfaceKstarParameterstrue
    (interfaceKstarParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceKstarParametersParticleData
    (interfaceKstarParameters,
     "ParticleData",
     "Use the value from the particle data objects",
     false);

  static ParVector<KPiKStarCurrent,double> interfaceKMagnitude
    ("KMagnitude",
     "Magnitude of the weight of the different resonances for the K pi channel",
     &KPiKStarCurrent::_kmag, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<KPiKStarCurrent,double> interfaceKPhase
    ("KPhase",
     "Phase of the weight of the different resonances for the K pi channel",
     &KPiKStarCurrent::_kphase, -1, 0., 0, 0,
     false, false, Interface::nolimits);
  
  static Switch<KPiKStarCurrent,int> interfaceKModel
    ("KModel",
     "The model to use for the propagator for the kaon modes.",
     &KPiKStarCurrent::_kmodel, 0, false, false);
  static SwitchOption interfaceKModelKuhn
    (interfaceKModel,
     "Kuhn",
     "The model of Kuhn and Santamaria",
     0);
  static SwitchOption interfaceKModelGounaris
    (interfaceKModel,
     "Gounaris",
     "The model of Gounaris and Sakurai.",
     1);

  static ClassDocumentation<KPiKStarCurrent> documentation
    ("The KPiKStarCurrent class is designed to implement weak"
     "decay to two scalar mesons using the models of either Kuhn and "
     "Santamaria (Z. Phys. C48, 445 (1990)) or Gounaris and Sakurai Phys. Rev. "
     "Lett. 21, 244 (1968).  The mixing parameters are taken from "
     "Phys. Rev. D61:112002,2000 (CLEO), although the PDG values for the "
     "masses and widths are used, for the decay pi+/- pi0."
     " The decay K pi is assumed to  be dominated by the lowest lying K* resonance.",
     "The weak "
     "decay current to two scalar mesons is implemented "
     "using the models of either Kuhn and "
     "Santamaria \\cite{Kuhn:1990ad} or Gounaris and Sakurai \\cite{Gounaris:1968mw}. "
     "The mixing parameters are taken from "
     "\\cite{Asner:1999kj}, although the PDG values for the "
     "masses and widths are used, for the decay pi+/- pi0."
     " The decay K pi is assumed to  be dominated by the lowest lying K* resonance.",
     "%\\cite{Kuhn:1990ad}\n"
     "\\bibitem{Kuhn:1990ad}\n"
     "  J.~H.~Kuhn and A.~Santamaria,\n"
     "  %``Tau decays to pions,''\n"
     "  Z.\\ Phys.\\  C {\\bf 48}, 445 (1990).\n"
     "  %%CITATION = ZEPYA,C48,445;%%\n"
     "%\\cite{Gounaris:1968mw}\n"
     "\\bibitem{Gounaris:1968mw}\n"
     "  G.~J.~Gounaris and J.~J.~Sakurai,\n"
     "   ``Finite width corrections to the vector meson dominance prediction for rho\n"
     "  %$\\to$ e+ e-,''\n"
     "  Phys.\\ Rev.\\ Lett.\\  {\\bf 21}, 244 (1968).\n"
     "  %%CITATION = PRLTA,21,244;%%\n"
     "%\\cite{Asner:1999kj}\n"
     "\\bibitem{Asner:1999kj}\n"
     "  D.~M.~Asner {\\it et al.}  [CLEO Collaboration],\n"
     "   ``Hadronic structure in the decay tau- --> nu/tau pi- pi0 pi0 and the  sign\n"
     "  %of the tau neutrino helicity,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 61}, 012002 (2000)\n"
     "  [arXiv:hep-ex/9902022].\n"
     "  %%CITATION = PHRVA,D61,012002;%%\n"
     );

}

// complete the construction of the decay mode for integration
bool KPiKStarCurrent::createMode(int icharge, tcPDPtr resonance,
				 FlavourInfo flavour,
				 unsigned int imode,PhaseSpaceModePtr mode,
				 unsigned int iloc,int ires,
				 PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(abs(icharge)!=3) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IHalf) return false;
  } 
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Half:
      if(icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusHalf:
      if(icharge == 3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown) {
    if(icharge== 3 and flavour.strange != Strangeness::PlusOne ) return false;
    if(icharge==-3 and flavour.strange != Strangeness::MinusOne) return false;
  }
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
  // make sure that the decays are kinematically allowed
  tPDPtr part[2];
  if(imode==0) {
    part[0]=getParticleData(ParticleID::Kplus);
    part[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==1) {
    part[0]=getParticleData(ParticleID::K0);
    part[1]=getParticleData(ParticleID::piplus);
  }
  else if(imode==2) {
    part[0]=getParticleData(ParticleID::eta);
    part[1]=getParticleData(ParticleID::Kplus);
  }
  Energy min(part[0]->massMin()+part[1]->massMin());
  if(min>upp) return false;
  // set the resonances
  // K+ pi0 or K0 pi+ or K eta decay
  tPDPtr res[3]={getParticleData(323),getParticleData(100323),getParticleData(30323)};
  if(icharge==-3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(res[ix]&&res[ix]->CC()) res[ix]=res[ix]->CC();
    }
  }
  // create the channels
  for(unsigned int ix=0;ix<3;++ix) {
    if(!res[ix]) continue;
    if(resonance && resonance != res[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],
		      ires+1,iloc+1,ires+1,iloc+2));
  }
  // reset the masses in the intergrators if needed
  if(_kstarparameters) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(ix<_kstarmasses.size()&&res[ix]) {
	mode->resetIntermediate(res[ix],_kstarmasses[ix],_kstarwidths[ix]);
      }
    }
  }
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector KPiKStarCurrent::particles(int icharge, unsigned int imode,
				     int,int) {
  tPDVector output(2);
  if(imode==0) {
    output[0]=getParticleData(ParticleID::Kplus);
    output[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==1) {
    output[0]=getParticleData(ParticleID::K0);
    output[1]=getParticleData(ParticleID::piplus);
  }
  else if(imode==2) {
    output[0]=getParticleData(ParticleID::eta);
    output[1]=getParticleData(ParticleID::Kplus);
  }
  if(icharge==-3) {
    for(unsigned int ix=0;ix<output.size();++ix) {
      if(output[ix]->CC()) output[ix]=output[ix]->CC();
    }
  }
  return output;
}

// hadronic current   
vector<LorentzPolarizationVectorE> 
KPiKStarCurrent::current(tcPDPtr resonance,
		    FlavourInfo flavour,
		    const int imode, const int ichan,Energy & scale, 
		    const tPDVector & outgoing,
		    const vector<Lorentz5Momentum> & momenta,
		    DecayIntegrator::MEOption) const {
  useMe();
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IHalf)
    return vector<LorentzPolarizationVectorE>();
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Half:
      if(icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusHalf:
      if(icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown) {
    if(icharge== 3 and flavour.strange != Strangeness::PlusOne ) return vector<LorentzPolarizationVectorE>();
    if(icharge==-3 and flavour.strange != Strangeness::MinusOne) return vector<LorentzPolarizationVectorE>();
  }
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return vector<LorentzPolarizationVectorE>();
  // momentum difference and sum of the mesons
  Lorentz5Momentum pdiff(momenta[0]-momenta[1]);
  Lorentz5Momentum psum (momenta[0]+momenta[1]);
  psum.rescaleMass();
  scale=psum.mass();
  // mass2 of vector intermediate state
  Energy2 q2(psum.m2());
  double dot(psum*pdiff/q2);
  psum *=dot;
  LorentzPolarizationVector vect;
  // calculate the current
  unsigned int imin=0, imax=_kwgt.size();
  if(ichan>0) {
    imin = ichan;
    imax = ichan+1;
  }
  if(resonance) {
    switch(abs(resonance->id())/1000) {
    case 0:
      imin=0; break;
    case 100:
      imin=1; break;
    case  30:
      imin=2; break;
    default:
      assert(false);
    }
    imax = imin+1;
  }
  Complex denom=std::accumulate(_kwgt.begin(),_kwgt.end(),Complex(0.));
  Complex FK(0.);
  for(unsigned int ix=imin;ix<imax;++ix) {
    FK+=_kwgt[ix]*BreitWigner(q2,_kmodel,ix);
  }
  // additional prefactors
  FK/=denom;
  // single kaon/pion modes
  if     (imode==0)      FK *= sqrt(0.5);
  else if(imode==1)      FK *= 1.       ;
  // the kaon eta mode
  else if(imode==2)      FK *=sqrt(1.5);
  // compute the current
  pdiff-=psum;
  return vector<LorentzPolarizationVectorE>(1,FK*pdiff);
}
   
bool KPiKStarCurrent::accept(vector<int> id) {
  // check there are only two particles
  if(id.size()!=2) return false;
  // single charged kaon
  if((abs(id[0])==ParticleID::Kplus  &&     id[1] ==ParticleID::pi0  ) ||
     (    id[0] ==ParticleID::pi0    && abs(id[1])==ParticleID::Kplus))
    return true;
  // single neutral kaon
  else if((id[0]==ParticleID::piminus && id[1]==ParticleID::Kbar0)   ||
	  (id[0]==ParticleID::Kbar0   && id[1]==ParticleID::piminus) ||
	  (id[0]==ParticleID::piplus  && id[1]==ParticleID::K0)      ||
	  (id[0]==ParticleID::K0      && id[1]==ParticleID::piplus))
    return true;
  // charged kaon and eta
  else if((id[0]==ParticleID::Kminus && id[1]==ParticleID::eta)    ||
	  (id[0]==ParticleID::eta    && id[1]==ParticleID::Kminus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::eta)    ||
	  (id[0]==ParticleID::eta    && id[1]==ParticleID::Kplus))
    return true;
  else
    return false;
}

// the decay mode
unsigned int KPiKStarCurrent::decayMode(vector<int> idout) {
  unsigned int imode(0),nkaon(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::K0) {
      imode=1;
      ++nkaon;
    }
    else if (abs(idout[ix])==ParticleID::Kplus) {
      imode=0;
      ++nkaon;
    }
    else if (idout[ix]==ParticleID::eta) {
      imode=2;
      break;
    }
  }
  return imode;
}

// output the information for the database
void KPiKStarCurrent::dataBaseOutput(ofstream & output,bool header,
					     bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::KPiKStarCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  unsigned int ix;
  for(ix=0;ix<_kstarmasses.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarMasses " << ix << " " << _kstarmasses[ix]/MeV << "\n";
  }
  for(ix=0;ix<_kstarwidths.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarWidths " << ix << " " << _kstarwidths[ix]/MeV << "\n";
  }
  output << "newdef " << name() << ":KstarParameters " << _kstarparameters << "\n";
  for(ix=0;ix<_kwgt.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KMagnitude " << ix << " " << _kmag[ix]   << "\n";
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KPhase "     << ix << " " << _kphase[ix] << "\n";
  }
  output << "newdef " << name() << ":KModel  " << _kmodel  << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./EtaPhotonCurrent.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPhotonCurrent class.
//

#include "EtaPhotonCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/HelicityFunctions.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"

using namespace Herwig;

EtaPhotonCurrent::EtaPhotonCurrent() {
  // modes handled
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
  // Masses for the resonances
  resMasses_ = {0.77526*GeV,0.78284*GeV,1.01952*GeV,1.465*GeV,1.70*GeV};
  // widths for the resonances
  resWidths_ = {0.1491 *GeV,0.00868*GeV,0.00421*GeV,0.40*GeV,0.30*GeV};
  // amplitudes
  amp_   = {0.0861/GeV,0.00824/GeV,0.0158/GeV,0.0147/GeV,ZERO};
  // phases
  phase_ = {0.,11.3,170.,61.,0.};
}

IBPtr EtaPhotonCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr EtaPhotonCurrent::fullclone() const {
  return new_ptr(*this);
}

void EtaPhotonCurrent::doinit() {
  WeakCurrent::doinit();
  assert(phase_.size()==amp_.size());
  couplings_.clear();
  Complex ii(0.,1.);
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    double phi = phase_[ix]/180.*Constants::pi;
    couplings_.push_back(amp_[ix]*(cos(phi)+ii*sin(phi)));
  }
  mpi_ = getParticleData(ParticleID::piplus)->mass();
}

void EtaPhotonCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(resMasses_,GeV) << ounit(resWidths_,GeV)
     << ounit(amp_,1./GeV) << phase_ << ounit(couplings_,1./GeV)
     << ounit(mpi_,GeV);
}

void EtaPhotonCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(resMasses_,GeV) >> iunit(resWidths_,GeV)
     >> iunit(amp_,1./GeV) >> phase_ >> iunit(couplings_,1./GeV)
     >> iunit(mpi_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EtaPhotonCurrent,WeakCurrent>
describeHerwigEtaPhotonCurrent("Herwig::EtaPhotonCurrent",
				"HwWeakCurrents.so");

void EtaPhotonCurrent::Init() {

  static ClassDocumentation<EtaPhotonCurrent> documentation
    ("The EtaPhotonCurrent class implements a current based"
     " on the model of SND for pion+photon",
     "The current based on the model of \\cite{Achasov:2006dv}"
     " for eta and photon was used.",
     "\\bibitem{Achasov:2006dv}\n"
     "M.~N.~Achasov {\\it et al.},\n"
     "%``Study of the e+ e- ---> eta gamma process with SND detector at the VEPP-2M e+ e- collider,''\n"
     "Phys.\\ Rev.\\ D {\\bf 74} (2006) 014016\n"
     "doi:10.1103/PhysRevD.74.014016\n"
     "[hep-ex/0605109].\n"
     "%%CITATION = doi:10.1103/PhysRevD.74.014016;%%\n"
     "%25 citations counted in INSPIRE as of 23 Aug 2018\n");


  static ParVector<EtaPhotonCurrent,Energy> interfaceResonanceMasses
    ("ResonanceMasses",
     "The masses of the resonances for the form factor",
     &EtaPhotonCurrent::resMasses_, GeV, 5, 775.26*MeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<EtaPhotonCurrent,Energy> interfaceResonanceWidths
    ("ResonanceWidths",
     "The widths of the resonances for the form factor",
     &EtaPhotonCurrent::resWidths_, GeV, 5, 149.1*MeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<EtaPhotonCurrent,InvEnergy> interfaceAmplitude
    ("Amplitude",
     "The amplitudes of the couplings",
     &EtaPhotonCurrent::amp_, 1./GeV, 5, 1./GeV, 0.0/GeV, 100./GeV,
     false, false, Interface::limited);

  static ParVector<EtaPhotonCurrent,double> interfacePhase
    ("Phase",
     "The phases of the couplings in degrees",
     &EtaPhotonCurrent::phase_, 5, 0., 0.0, 360.0,
     false, false, Interface::limited);

}

// complete the construction of the decay mode for integration
bool EtaPhotonCurrent::createMode(int icharge, tcPDPtr resonance,
				   FlavourInfo flavour,
				   unsigned int, PhaseSpaceModePtr mode,
				   unsigned int iloc,int ires,
				   PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(icharge!=0) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IZero) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    if(flavour.I3!=IsoSpin::I3Zero) return false;
  }
  if(flavour.strange != Strangeness::Unknown)
     if(flavour.strange != Strangeness::Zero and flavour.strange != Strangeness::ssbar) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       ) return false;
  // check that the mode is are kinematical allowed
  Energy min = getParticleData(ParticleID::eta)->mass();
  if(min>upp) return false;
  // resonances for the intermediate channels
  tPDVector res = {getParticleData(113),
		   getParticleData(   223),
		   getParticleData(   333),
		   getParticleData(100113),
		   getParticleData( 100333)};
  // set up the integration channels;
  for(unsigned int ix=0;ix<res.size();++ix) {
    if(resonance && resonance!=res[ix]) continue;
    if(flavour.strange != Strangeness::Unknown) {
      if     (flavour.strange == Strangeness::Zero && (ix==2||ix==4)) {
	continue;
      }
      else if(flavour.strange == Strangeness::ssbar && (ix<2||ix==3)) {
	continue;
      }
    }
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],
 		      ires+1,iloc+1,ires+1,iloc+2));
  }
  // reset the masses and widths of the resonances if needed
  for(unsigned int ix=0;ix<res.size();++ix) {
    mode->resetIntermediate(res[ix],resMasses_[ix],resWidths_[ix]);
  }
  return true;
}

// the particles produced by the current
tPDVector EtaPhotonCurrent::particles(int icharge, unsigned int imode,int,int) {
  assert(icharge==0 && imode<=1);
  return {getParticleData(ParticleID::eta),getParticleData(ParticleID::gamma)};
}

void EtaPhotonCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) ++ix;
    temp[ix] = HelicityFunctions::polarizationVector(-decay[1]->momentum(),
						     ix,Helicity::outgoing);
  }
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  VectorWaveFunction::constructSpinInfo(temp,decay[1],
					outgoing,true,true);
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
EtaPhotonCurrent::current(tcPDPtr resonance,
			   FlavourInfo flavour,
			   const int, const int ichan, Energy & scale, 
			   const tPDVector & ,
			   const vector<Lorentz5Momentum> & momenta,
			   DecayIntegrator::MEOption) const {
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IZero) return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    if(flavour.I3!=IsoSpin::I3Zero) return vector<LorentzPolarizationVectorE>();
  }
  if(flavour.strange != Strangeness::Unknown)
     if(flavour.strange != Strangeness::Zero and flavour.strange != Strangeness::ssbar)
       return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return vector<LorentzPolarizationVectorE>();
  useMe();
  // polarization vectors of the photon
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) ++ix;
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // total momentum of the system
  Lorentz5Momentum q(momenta[0]+momenta[1]);
  // overall hadronic mass
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2(q.m2());
  unsigned int imin = 0;
  unsigned int imax = couplings_.size();
  if(flavour.strange != Strangeness::Unknown) {
    if     (flavour.strange == Strangeness::Zero) {
      if(imin==2) imin=4;
    }
    else if(flavour.strange == Strangeness::ssbar) {
      imin=2*imin+2;
    } 
  }
  if(ichan>0) {
    imin = ichan;
    imax = imin+1;
  }
  if(resonance) {
    switch(abs(resonance->id())) {
    case 113: case 213 :
      imin=0;
      break;
    case 223:
      imin=1;
      break;
    case 333:
      imin=2;
      break;
    case 100213:
      imin = 3;
      break;
    case 100333 :
      imin = 4;
      break;
    default:
      assert(false);
    }
    imax=imin+1;
  }
  // compute the form factor
  complex<InvEnergy> formFactor(ZERO);
  // loop over the resonances
  for(unsigned int ix=imin;ix<imax;++ix) {
    if(flavour.strange != Strangeness::Unknown) {
      if     (flavour.strange == Strangeness::Zero && (ix==2||ix==4)) continue;
      else if(flavour.strange == Strangeness::ssbar && (ix<2||ix==3)) continue; 
    }
    Energy2 mR2(sqr(resMasses_[ix]));
    // compute the width
    Energy width(ZERO);
    // rho
    if(ix==0) {
      width = resWidths_[0]*mR2/q2*pow(max(double((q2-4.*sqr(mpi_))/(mR2-4.*sqr(mpi_))),0.),1.5);
    }
    else {
      width = resWidths_[ix];
    }
    formFactor += couplings_[ix]*mR2/(mR2-q2-Complex(0.,1.)*q.mass()*width);
  }
  // calculate the current
  vector<LorentzPolarizationVectorE> ret(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) continue;
    ret[ix] += formFactor*Helicity::epsilon(q,temp[ix],momenta[1]);
  }
  return ret;
}

bool EtaPhotonCurrent::accept(vector<int> id) {
  if(id.size()!=2) return false;
  unsigned int neta(0),ngamma(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::eta)   ++neta;
    else if(id[ix]==ParticleID::gamma) ++ngamma;
  }
  return ngamma == 1 && neta==1;
}

unsigned int EtaPhotonCurrent::decayMode(vector<int>) {
  return 0;
}

// output the information for the database
void EtaPhotonCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::EtaPhotonCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<resMasses_.size();++ix) {
    if(ix<5) output << "newdef " << name() << ":ResonanceMasses " << ix 
		    << " " << resMasses_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":ResonanceMasses " << ix 
		    << " " << resMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<resWidths_.size();++ix) {
    if(ix<5) output << "newdef " << name() << ":ResonanceWidths " << ix 
		    << " " << resWidths_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":ResonanceWidths " << ix 
		    << " " << resWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    if(ix<5) output << "newdef " << name() << ":Amplitude " << ix 
		    << " " << amp_[ix]*GeV << "\n";
    else     output << "insert " << name() << ":Amplitude " << ix 
		    << " " << amp_[ix]*GeV << "\n";
  }
  for(unsigned int ix=0;ix<phase_.size();++ix) {
    if(ix<5) output << "newdef " << name() << ":Phase " << ix 
		    << " " << phase_[ix] << "\n";
    else     output << "insert " << name() << ":Phase " << ix 
		    << " " << phase_[ix] << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./PionPhotonCurrent.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PionPhotonCurrent class.
//

#include "PionPhotonCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/HelicityFunctions.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"

using namespace Herwig;

PionPhotonCurrent::PionPhotonCurrent() {
  // modes handled
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
  // Masses for the resonances
  resMasses_ = {0.77526*GeV,0.78284*GeV,1.45*GeV,1.70*GeV,1.01952*GeV};
  // widths for the resonances
  resWidths_ = {0.1491 *GeV,0.00868*GeV,0.40*GeV,0.30*GeV,0.00421*GeV};
  // amplitudes
  amp_   = {0.0426/GeV,0.0434/GeV,0.00523/GeV,ZERO,0.00303/GeV};
  // phases
  phase_ = {-12.7,0.,180.,0.,158.};
}

IBPtr PionPhotonCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr PionPhotonCurrent::fullclone() const {
  return new_ptr(*this);
}

void PionPhotonCurrent::doinit() {
  WeakCurrent::doinit();
  assert(phase_.size()==amp_.size());
  couplings_.clear();
  Complex ii(0.,1.);
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    double phi = phase_[ix]/180.*Constants::pi;
    couplings_.push_back(amp_[ix]*(cos(phi)+ii*sin(phi)));
  }
  mpi_ = getParticleData(ParticleID::piplus)->mass();
}

void PionPhotonCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(resMasses_,GeV) << ounit(resWidths_,GeV)
     << ounit(amp_,1./GeV) << phase_ << ounit(couplings_,1./GeV)
     << ounit(mpi_,GeV);
}

void PionPhotonCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(resMasses_,GeV) >> iunit(resWidths_,GeV)
     >> iunit(amp_,1./GeV) >> phase_ >> iunit(couplings_,1./GeV)
     >> iunit(mpi_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PionPhotonCurrent,WeakCurrent>
describeHerwigPionPhotonCurrent("Herwig::PionPhotonCurrent",
				"HwWeakCurrents.so");

void PionPhotonCurrent::Init() {

  static ClassDocumentation<PionPhotonCurrent> documentation
    ("The PionPhotonCurrent class implements a current based"
     " on the model of SND for pion+photon",
     "The current based on the model of \\cite{Achasov:2016bfr}"
     " for pion and photon was used.",
     "\\bibitem{Achasov:2016bfr}\n"
     "M.~N.~Achasov {\\it et al.} [SND Collaboration],\n"
     "%``Study of the reaction $e^+e^- \\to \\pi^0\\gamma$ with the SND detector at the VEPP-2M collider,''\n"
     "Phys.\\ Rev.\\ D {\\bf 93} (2016) no.9,  092001\n"
     "doi:10.1103/PhysRevD.93.092001\n"
     "[arXiv:1601.08061 [hep-ex]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.93.092001;%%\n"
     "%20 citations counted in INSPIRE as of 23 Aug 2018\n");

  static ParVector<PionPhotonCurrent,Energy> interfaceResonanceMasses
    ("ResonanceMasses",
     "The masses of the resonances for the form factor",
     &PionPhotonCurrent::resMasses_, GeV, 5, 775.26*MeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<PionPhotonCurrent,Energy> interfaceResonanceWidths
    ("ResonanceWidths",
     "The widths of the resonances for the form factor",
     &PionPhotonCurrent::resWidths_, GeV, 5, 149.1*MeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<PionPhotonCurrent,InvEnergy> interfaceAmplitude
    ("Amplitude",
     "The amplitudes of the couplings",
     &PionPhotonCurrent::amp_, 1./GeV, 5, 1./GeV, 0.0/GeV, 100./GeV,
     false, false, Interface::limited);

  static ParVector<PionPhotonCurrent,double> interfacePhase
    ("Phase",
     "The phases of the couplings in degrees",
     &PionPhotonCurrent::phase_, 5, 0., -360.0, 360.0,
     false, false, Interface::limited);

}

// complete the construction of the decay mode for integration
bool PionPhotonCurrent::createMode(int icharge, tcPDPtr resonance,
				   FlavourInfo flavour,
				   unsigned int imode,PhaseSpaceModePtr mode,
				   unsigned int iloc,int ires,
				   PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode >= 1))
    return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(imode==0) {
      if(flavour.I!=IsoSpin::IOne) return false;
    }
    else {
      if(flavour.I!=IsoSpin::IOne &&
	 flavour.I!=IsoSpin::IZero) return false;
    }
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode!=1) return false;
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown) {
    if(imode==0 and flavour.strange != Strangeness::Zero) return false;
    if(imode==1 and flavour.strange != Strangeness::Zero and flavour.strange != Strangeness::ssbar) return false;
  }
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       ) return false;
  // check that the mode is are kinematical allowed
  Energy min = imode==0 ?
    getParticleData(ParticleID::piplus)->mass() :
    getParticleData(ParticleID::pi0   )->mass();
  if(min>upp) return false;
  // resonances for the intermediate channels
  tPDVector res;
  if(imode==0) {
    if(icharge==-3) res.push_back(getParticleData(-213));
    else            res.push_back(getParticleData( 213));
  }
  else {
    if(flavour.I==IsoSpin::IUnknown||flavour.I==IsoSpin::IOne)
      res.push_back(getParticleData(113));
    if(flavour.I==IsoSpin::IUnknown||flavour.I==IsoSpin::IZero) {
      res.push_back(getParticleData(   223));
      res.push_back(getParticleData(100223));
      res.push_back(getParticleData( 30223));
      res.push_back(getParticleData(   333));
    }
  }
  // set up the integration channels;
  for(unsigned int ix=0;ix<res.size();++ix) {
    if(resonance && resonance!=res[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],
 		      ires+1,iloc+1,ires+1,iloc+2));
  }
  // reset the masses and widths of the resonances if needed
  for(unsigned int ix=0;ix<res.size();++ix) {
    int ires(0);
    if(res[ix]->id()==223)         ires=1;
    else if(res[ix]->id()==100223) ires=2;
    else if(res[ix]->id()== 30223) ires=3;
    else if(res[ix]->id()==   333) ires=4;
    mode->resetIntermediate(res[ix],resMasses_[ires],resWidths_[ires]);
  }
  return true;
}

// the particles produced by the current
tPDVector PionPhotonCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart = {tPDPtr(),
 		       getParticleData(ParticleID::gamma)};
  if(imode==0) {
    if(icharge==3)       extpart[0] = getParticleData(ParticleID::piplus );
    else if(icharge==-3) extpart[0] = getParticleData(ParticleID::piminus);
  }
  else {
    extpart[0] = getParticleData(ParticleID::pi0);
  }
  return extpart;
}

void PionPhotonCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) ++ix;
    temp[ix] = HelicityFunctions::polarizationVector(-decay[1]->momentum(),
						     ix,Helicity::outgoing);
  }
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  VectorWaveFunction::constructSpinInfo(temp,decay[1],
					outgoing,true,true);
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
PionPhotonCurrent::current(tcPDPtr resonance,
			   FlavourInfo flavour,
			   const int imode, const int ichan, Energy & scale, 
			   const tPDVector & outgoing,
			   const vector<Lorentz5Momentum> & momenta,
			   DecayIntegrator::MEOption) const {
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode == 1))
    return vector<LorentzPolarizationVectorE>();
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(imode==0) {
      if(flavour.I!=IsoSpin::IOne) return vector<LorentzPolarizationVectorE>();
    }
    else {
      if(flavour.I!=IsoSpin::IOne &&
	 flavour.I!=IsoSpin::IZero) return vector<LorentzPolarizationVectorE>();
    }
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode!=1) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  useMe();
  // polarization vectors of the photon
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) ++ix;
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // total momentum of the system
  Lorentz5Momentum q(momenta[0]+momenta[1]);
  // overall hadronic mass
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2(q.m2());
  unsigned int imin = 0;
  unsigned int imax = imode==0 ? 1 : 5;
  if(flavour.I==IsoSpin::IOne)
    imax = 1;
  else if(flavour.I==IsoSpin::IZero) {
    imin = 1;
  }
  if(ichan>0) {
    if(flavour.I==IsoSpin::IZero)
      imin = ichan+1;
    else
      imin = ichan;
    imax=imin+1;
  }
  if(resonance) {
    switch(abs(resonance->id())) {
    case 113: case 213 :
      imin=0;
      break;
    case 223:
      imin=1;
      break;
    case 333:
      imin=4;
      break;
    case 100223:
      imin = 2;
      break;
    case 30223 :
      imin = 3;
      break;
    default:
      assert(false);
    }
    imax=imin+1;
  }
  // compute the form factor
  complex<InvEnergy> formFactor(ZERO);
  // loop over the resonances
  for(unsigned int ix=imin;ix<imax;++ix) {
    if(ix==4 and flavour.strange==Strangeness::Zero) continue;
    if(ix<4  and flavour.strange==Strangeness::ssbar) continue;
    Energy2 mR2(sqr(resMasses_[ix]));
    // compute the width
    Energy width(ZERO);
    // rho
    if(ix==0) {
      width = resWidths_[0]*mR2/q2*pow(max(double((q2-4.*sqr(mpi_))/(mR2-4.*sqr(mpi_))),0.),1.5);
    }
    else {
      width = resWidths_[ix];
    }
    formFactor += couplings_[ix]*mR2/(mR2-q2-Complex(0.,1.)*q.mass()*width);
  }
  // calculate the current
  vector<LorentzPolarizationVectorE> ret(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) continue;
    ret[ix] += formFactor*Helicity::epsilon(q,temp[ix],momenta[1]);
  }
  return ret;
}

bool PionPhotonCurrent::accept(vector<int> id) {
  if(id.size()!=2) return false;
  unsigned int npiplus(0),npi0(0),ngamma(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::piplus) ++npiplus;
    else if(id[ix]==ParticleID::gamma)  ++ngamma;
    else if(id[ix]==ParticleID::pi0)    ++npi0;
  }
  return ngamma == 1 && (npiplus==1 || npi0==1);
}

unsigned int PionPhotonCurrent::decayMode(vector<int> id) {
  int npip(0),npim(0),npi0(0),ngamma(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)         ++npip;
    else if(id[ix]==ParticleID::piminus)   ++npim;
    else if(id[ix]==ParticleID::pi0)       ++npi0;
    else if(id[ix]==ParticleID::gamma)   ++ngamma;
  }
  if((npip==1 || npim == 1) && ngamma==1)
    return 0;
  else
    return 1;
}

// output the information for the database
void PionPhotonCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::PionPhotonCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<resMasses_.size();++ix) {
    if(ix<5) output << "newdef " << name() << ":ResonanceMasses " << ix 
		    << " " << resMasses_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":ResonanceMasses " << ix 
		    << " " << resMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<resWidths_.size();++ix) {
    if(ix<5) output << "newdef " << name() << ":ResonanceWidths " << ix 
		    << " " << resWidths_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":ResonanceWidths " << ix 
		    << " " << resWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    if(ix<5) output << "newdef " << name() << ":Amplitude " << ix 
		    << " " << amp_[ix]*GeV << "\n";
    else     output << "insert " << name() << ":Amplitude " << ix 
		    << " " << amp_[ix]*GeV << "\n";
  }
  for(unsigned int ix=0;ix<phase_.size();++ix) {
    if(ix<5) output << "newdef " << name() << ":Phase " << ix 
		    << " " << phase_[ix] << "\n";
    else     output << "insert " << name() << ":Phase " << ix 
		    << " " << phase_[ix] << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./TwoPionPhotonCurrent.cc"
// -*- C++ -*-
//
// TwoPionPhotonCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoPionPhotonCurrent class.
//
//  Author: Peter Richardson
//

#include "TwoPionPhotonCurrent.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

TwoPionPhotonCurrent::TwoPionPhotonCurrent() {
  // modes handled
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
  // weight of the resonances in the current
  _resweights = {1.0,-0.1};
  // parameters of the rho resonaces
  _rhomasses = {0.773*GeV,1.70*GeV};
  _rhowidths = {0.145*GeV,0.26*GeV};
  // parameters fo the omega resonance
  _omegamass=782*MeV;_omegawidth=8.5*MeV;
  // couplings
  _grho   = 0.11238947*GeV2;
  _grhoomegapi = 12.924/GeV;
  // parameters for the resonance used in the integration
  _intmass  = 1.2*GeV;
  _intwidth = 0.35*GeV;
}

void TwoPionPhotonCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(_grho,GeV2) << ounit(_grhoomegapi,1/GeV) << _resweights 
     << ounit(_rhomasses,GeV) << ounit(_rhowidths,GeV) << ounit(_omegamass,GeV) 
     << ounit(_omegawidth,GeV) << ounit(_intmass,GeV) 
     << ounit(_intwidth,GeV) ;
}

void TwoPionPhotonCurrent::persistentInput(PersistentIStream & is, int) { 
  is >> iunit(_grho,GeV2) >> iunit(_grhoomegapi,1/GeV) >> _resweights 
     >> iunit(_rhomasses,GeV) >> iunit(_rhowidths,GeV) >> iunit(_omegamass,GeV) 
     >> iunit(_omegawidth,GeV) >> iunit(_intmass,GeV)
     >> iunit(_intwidth,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TwoPionPhotonCurrent,WeakCurrent>
describeHerwigTwoPionPhotonCurrent("Herwig::TwoPionPhotonCurrent", "HwWeakCurrents.so");

void TwoPionPhotonCurrent::Init() {

  static ParVector<TwoPionPhotonCurrent,double> interfacereswgt
    ("Weights",
     "The weights of the different resonances for the decay tau -> nu pi pi gamma",
     &TwoPionPhotonCurrent::_resweights,
     0, 0, 0, -1000, 1000, false, false, true);

  static ParVector<TwoPionPhotonCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the decay tau ->  pi pi photon",
     &TwoPionPhotonCurrent::_rhomasses, MeV, -1, 773.*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<TwoPionPhotonCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the decay tau -> nu pi pi photon",
     &TwoPionPhotonCurrent::_rhowidths, MeV, -1, 145.*MeV, ZERO, 1000.*MeV,
     false, false, true);

  static Parameter<TwoPionPhotonCurrent,Energy> interfaceomegamass
    ("omegamass",
     "The mass of the omega",
     &TwoPionPhotonCurrent::_omegamass, GeV, 0.782*GeV, ZERO, 1.0*GeV,
     false, false, true);
  
  static Parameter<TwoPionPhotonCurrent,Energy> interfaceomegawidth
    ("omegawidth",
     "The width of the omega for the decay tau- -> pi pi photon",
     &TwoPionPhotonCurrent::_omegawidth, GeV, 0.0085*GeV, ZERO, 1.*GeV,
     false, false, false);
  
  static ClassDocumentation<TwoPionPhotonCurrent> documentation
    ("The TwoPionPhotonCurrent class implements the decay "
     "tau+/- -> pi+/- pi0 gamma via an omega.",
     "The decay $\\tau^\\pm \\to \\omega \\to \\pi^\\pm \\pi^0 \\gamma$ "
     "is modelled after \\cite{Jadach:1993hs}.",
     "  %\\cite{Jadach:1993hs}\n"
     "\\bibitem{Jadach:1993hs}\n"
     "  S.~Jadach, Z.~Was, R.~Decker and J.~H.~Kuhn,\n"
     "  %``The Tau Decay Library Tauola: Version 2.4,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 76}, 361 (1993).\n"
     "  %%CITATION = CPHCB,76,361;%%\n"
     );

  static Parameter<TwoPionPhotonCurrent,Energy2> interfacegrho
    ("grho",
     "The rho meson decay constant.",
     &TwoPionPhotonCurrent::_grho, GeV2, 0.11238947*GeV2, -1.*GeV2, 1.*GeV2,
     false, false, false);

  static Parameter<TwoPionPhotonCurrent,InvEnergy> interfacegrhoomegapi
    ("grhoomegapi",
     "The rho-omega-pi coupling",
     &TwoPionPhotonCurrent::_grhoomegapi, 1./GeV, 12.924/GeV,
     -100./GeV, 100./GeV,
     false, false, false);

  static Parameter<TwoPionPhotonCurrent,Energy> interfaceIntegrationMass
    ("IntegrationMass",
     "Mass of the pseudoresonance used to improve integration effciency",
     &TwoPionPhotonCurrent::_intmass, GeV, 1.4*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<TwoPionPhotonCurrent,Energy> interfaceIntegrationWidth
    ("IntegrationWidth",
     "Width of the pseudoresonance used to improve integration effciency",
     &TwoPionPhotonCurrent::_intwidth, GeV, 0.5*GeV, ZERO, 10.0*GeV,
     false, false, true);
}

// complete the construction of the decay mode for integration
bool TwoPionPhotonCurrent::createMode(int icharge, tcPDPtr resonance,
				      FlavourInfo flavour,
				      unsigned int imode,PhaseSpaceModePtr mode,
				      unsigned int iloc,int ires,
				      PhaseSpaceChannel phase, Energy upp ) {
  assert(!resonance);
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode >= 1))
    return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode!=1) return false;
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return false;
  // check that the mode is are kinematical allowed
  Energy min(getParticleData(ParticleID::piplus)->mass()+
  	     getParticleData(ParticleID::pi0   )->mass());
  if(min>upp) return false;
  // set up the integration channels;
  tPDPtr omega(getParticleData(ParticleID::omega));
  tPDPtr rho;
  if(icharge==-3)     rho = getParticleData(-213);
  else if(icharge==0) rho = getParticleData( 113);
  else if(icharge==3) rho = getParticleData( 213);
  mode->addChannel((PhaseSpaceChannel(phase),ires,rho,
		    ires+1,omega,ires+1,iloc+1,
		    ires+2,iloc+2,ires+2,iloc+3));
  // reset the masses and widths of the resonances if needed
  mode->resetIntermediate(rho,_intmass,_intwidth);
  // set up the omega masses and widths
  mode->resetIntermediate(omega,_omegamass,_omegawidth);
  return true;
}

// the particles produced by the current
tPDVector TwoPionPhotonCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart = {tPDPtr(),
		       getParticleData(ParticleID::pi0),
		       getParticleData(ParticleID::gamma)};
  if(imode==0) {
    if(icharge==3)       extpart[0] = getParticleData(ParticleID::piplus );
    else if(icharge==-3) extpart[0] = getParticleData(ParticleID::piminus);
  }
  else {
    extpart[0] = getParticleData(ParticleID::pi0);
  }
  return extpart;
}

void TwoPionPhotonCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) ++ix;
    temp[ix] = HelicityFunctions::polarizationVector(-decay[2]->momentum()
						     ,ix,Helicity::outgoing);
  }
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
  VectorWaveFunction::constructSpinInfo(temp,decay[2],
					outgoing,true,true);
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
TwoPionPhotonCurrent::current(tcPDPtr resonance,
			      FlavourInfo flavour,
			      const int imode, const int ,Energy & scale, 
			      const tPDVector & outgoing,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  assert(!resonance);
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge()+outgoing[2]->iCharge();
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode == 1))
    return vector<LorentzPolarizationVectorE>();
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode!=1) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return vector<LorentzPolarizationVectorE>();
  useMe();
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) ++ix;
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[2],ix,Helicity::outgoing);
  }
  // locate the particles
  Lorentz5Momentum pout(momenta[1]+momenta[2]+momenta[0]);
  // overall hadronic mass
  pout.rescaleMass();
  scale=pout.mass();
  Energy2 q2(pout.m2());
  // mass of the omega
  pout = momenta[1]+momenta[2];
  pout.rescaleMass();
  Energy2 s2(pout.m2());
  // compute the prefactor
  complex<InvEnergy3> prefactor(-FFunction(ZERO)*FFunction(q2)*scale*
				sqrt(Constants::twopi*generator()->standardModel()->alphaEM())*
				BreitWigner(s2,10));
  // dot products which don't depend on the polarization vector
  Energy2 dot12(momenta[2]*momenta[1]);
  Energy2 dot13(momenta[2]*momenta[0]);
  Energy2 dot23(momenta[1]*momenta[0]);
  Energy2 mpi2 = sqr(momenta[0].mass());
  vector<LorentzPolarizationVectorE> ret(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix!=1) {
      // obtain the dot products we need
      complex<Energy> dote2 = temp[ix]*momenta[1];
      complex<Energy> dote3 = temp[ix]*momenta[0];
      // now compute the coefficients
      complex<Energy4> coeffa = mpi2*dot13-dot12*(dot23-dot13);
      complex<Energy3> coeffb = dote2*dot13-dote3*dot12;
      complex<Energy3> coeffc = dote2*dot23-dote3*(mpi2+dot12);
      // finally compute the current
      ret[ix]= prefactor*(coeffa*temp[ix]
			  -coeffb*momenta[1]
			  +coeffc*momenta[2]);
    }
    else
      ret[ix]=LorentzPolarizationVectorE();
  }
  return ret;
}

bool TwoPionPhotonCurrent::accept(vector<int> id) {
  if(id.size()!=3){return false;}
  unsigned int npiplus(0),npi0(0),ngamma(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::piplus) ++npiplus;
    else if(id[ix]==ParticleID::gamma)  ++ngamma;
    else if(id[ix]==ParticleID::pi0)    ++npi0;
  }
  return (npiplus==1&&ngamma==1&&npi0==1) ||
    (npi0==2&&ngamma==1);
}

unsigned int TwoPionPhotonCurrent::decayMode(vector<int> id) {
  int npip(0),npim(0),npi0(0),ngamma(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)         ++npip;
    else if(id[ix]==ParticleID::piminus)   ++npim;
    else if(id[ix]==ParticleID::pi0)       ++npi0;
    else if(id[ix]==ParticleID::gamma)   ++ngamma;
  }
  if((npip==1 || npim == 1) && npi0==1 && ngamma==1)
    return 0;
  else
    return 1;
}

// output the information for the database
void TwoPionPhotonCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::TwoPionPhotonCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  output << "newdef " << name() << ":omegamass "    << _omegamass/GeV << "\n";
  output << "newdef " << name() << ":omegawidth "    << _omegawidth/GeV << "\n";
  output << "newdef " << name() << ":grho "    << _grho/GeV2 << "\n";
  output << "newdef " << name() << ":grhoomegapi "    << _grhoomegapi*GeV << "\n";
  output << "newdef " << name() << ":IntegrationMass "  << _intmass/GeV  << "\n";
  output << "newdef " << name() << ":IntegrationWidth " << _intwidth/GeV  << "\n";
  unsigned int ix;
  for(ix=0;ix<_resweights.size();++ix) {
    if(ix<2) output << "newdef " << name() << ":Weights " << ix 
		    << " " << _resweights[ix] << "\n";
    else     output << "insert " << name() << ":Weights " << ix 
		    << " " << _resweights[ix] << "\n";
  }
  for(ix=0;ix<_rhomasses.size();++ix) {
    if(ix<2) output << "newdef " << name() << ":RhoMasses " << ix 
		    << " " << _rhomasses[ix]/MeV << "\n";
    else     output << "insert " << name() << ":RhoMasses " << ix 
		    << " " << _rhomasses[ix]/MeV << "\n";
  }
  for(ix=0;ix<_rhowidths.size();++ix) {
    if(ix<2) output << "newdef " << name() << ":RhoWidths " << ix 
		    << " " << _rhowidths[ix]/MeV << "\n";
    else     output << "insert " << name() << ":RhoWidths " << ix 
		    << " " << _rhowidths[ix]/MeV << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./TwoPionPhotonSNDCurrent.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoPionPhotonSNDCurrent class.
//

#include "TwoPionPhotonSNDCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/epsilon.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"

using namespace Herwig;

TwoPionPhotonSNDCurrent::TwoPionPhotonSNDCurrent() {
  // modes handled
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
  // amplitudes for the weights in the current
  amp_   = {1.,0.175,0.014};
  phase_ = {0.,124.,-63.};
  // rho masses and widths
  rhoMasses_ = {0.77526*GeV,1.510*GeV,1.720*GeV};
  rhoWidths_ = {0.1491 *GeV,0.44 *GeV,0.25 *GeV};
  // coupling
  gRhoOmegaPi_   = 15.9/GeV;
  fRho_        = 4.9583;
  gGammaOmegaPi_ = 0.695821538653/GeV;
  fRho_        = 4.9583;
  // omega parameters
  omegaMass_  = 782.65*MeV;
  omegaWidth_ = 8.49 *MeV;
}

IBPtr TwoPionPhotonSNDCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr TwoPionPhotonSNDCurrent::fullclone() const {
  return new_ptr(*this);
}

void TwoPionPhotonSNDCurrent::doinit() {
  WeakCurrent::doinit();
  assert(phase_.size()==amp_.size());
  wgts_.clear();
  Complex ii(0.,1.);
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    double phi = phase_[ix]/180.*Constants::pi;
    wgts_.push_back(amp_[ix]*(cos(phi)+ii*sin(phi)));
  }
  mpi_ = getParticleData(ParticleID::piplus)->mass();
}

void TwoPionPhotonSNDCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV)
     << amp_ << phase_ << wgts_ << fRho_
     << ounit(gRhoOmegaPi_,1./GeV) << ounit(gGammaOmegaPi_,1./GeV)
     << ounit(omegaMass_,GeV) << ounit(omegaWidth_,GeV) << ounit(mpi_,GeV);
}

void TwoPionPhotonSNDCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV)
     >> amp_ >> phase_ >> wgts_ >> fRho_
     >> iunit(gRhoOmegaPi_,1./GeV) >> iunit(gGammaOmegaPi_,1./GeV)
     >> iunit(omegaMass_,GeV) >> iunit(omegaWidth_,GeV) >> iunit(mpi_,GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TwoPionPhotonSNDCurrent,WeakCurrent>
describeHerwigTwoPionPhotonSNDCurrent("Herwig::TwoPionPhotonSNDCurrent",
				      "HwWeakCurrents.so");

void TwoPionPhotonSNDCurrent::Init() {

  static ClassDocumentation<TwoPionPhotonSNDCurrent> documentation
    ("The TwoPionPhotonSNDCurrent class provides the weka current for"
     "pi pi gamma using the model of SND",
     "The current based on \\cite{Achasov:2016zvn} for $\\pi\\pi^0\\gamma$ was used.\n",
     "\\bibitem{Achasov:2016zvn}"
     "M.~N.~Achasov {\\it et al.},\n"
     "%``Updated measurement of the $e^+e^- \\to \\omega \\pi^0 \\to \\pi^0\\pi^0\\gamma$ cross section with the SND detector,''\n"
     "Phys.\\ Rev.\\ D {\\bf 94} (2016) no.11,  112001\n"
     "doi:10.1103/PhysRevD.94.112001\n"
     "[arXiv:1610.00235 [hep-ex]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.94.112001;%%\n"
     "%12 citations counted in INSPIRE as of 22 Aug 2018\n");

  static ParVector<TwoPionPhotonSNDCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the rho mesons",
     &TwoPionPhotonSNDCurrent::rhoMasses_, GeV, -1, 775.26*MeV,
     0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<TwoPionPhotonSNDCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the rho mesons",
     &TwoPionPhotonSNDCurrent::rhoWidths_, GeV, -1, 0.1491*GeV,
     0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<TwoPionPhotonSNDCurrent,double> interfaceAmplitudes
    ("Amplitudes",
     "THe amplitudes for the different rho resonances",
     &TwoPionPhotonSNDCurrent::amp_, -1, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static ParVector<TwoPionPhotonSNDCurrent,double> interfacePhase
    ("Phase",
     "The phases for the different rho resonances in degrees",
     &TwoPionPhotonSNDCurrent::phase_, -1, 0., 0.0, 360.,
     false, false, Interface::limited);

  static Parameter<TwoPionPhotonSNDCurrent,double> interfacefRho
    ("fRho",
     "The coupling of the photon and the rho meson",
     &TwoPionPhotonSNDCurrent::fRho_, 4.9583, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<TwoPionPhotonSNDCurrent,InvEnergy> interfacegRhoOmegaPi
    ("gRhoOmegaPi",
     "The coupling rho-omega-pi",
     &TwoPionPhotonSNDCurrent::gRhoOmegaPi_, 1./GeV,
     15.9/GeV, 0./GeV, 1000./GeV,
     false, false, Interface::limited);

  static Parameter<TwoPionPhotonSNDCurrent,InvEnergy> interfacegGammaOmegaPi
    ("gGammaOmegaPi",
     "The coupling gamma-omega-pi",
     &TwoPionPhotonSNDCurrent::gGammaOmegaPi_, 1./GeV,
     0.695821538653/GeV, 0./GeV, 1000./GeV,
     false, false, Interface::limited);

  static Parameter<TwoPionPhotonSNDCurrent,Energy> interfaceOmegaMass
    ("OmegaMass",
     "The mass of the omega meson",
     &TwoPionPhotonSNDCurrent::omegaMass_, GeV, 0.78265*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<TwoPionPhotonSNDCurrent,Energy> interfaceOmegaWidth
    ("OmegaWidth",
     "The width of the omega meson",
     &TwoPionPhotonSNDCurrent::omegaWidth_, GeV, 8.49*MeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
}

// complete the construction of the decay mode for integration
bool TwoPionPhotonSNDCurrent::createMode(int icharge, tcPDPtr resonance,
					 FlavourInfo flavour,
					 unsigned int imode,PhaseSpaceModePtr mode,
					 unsigned int iloc,int ires,
					 PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode >= 1))
    return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode!=1) return false;
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return false;
  // check that the mode is are kinematical allowed
  Energy min(getParticleData(ParticleID::piplus)->mass()+
  	     getParticleData(ParticleID::pi0   )->mass());
  if(min>upp) return false;
  // set up the integration channels;
  tPDPtr omega(getParticleData(ParticleID::omega));
  vector<tPDPtr> rho;
  if(icharge==-3)
    rho = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  else if(icharge==0)
    rho = {getParticleData( 113),getParticleData( 100113),getParticleData( 30113)};
  else if(icharge==3)
    rho = {getParticleData( 213),getParticleData( 100213),getParticleData( 30213)};
  for(unsigned int ix=0;ix<3;++ix) {
    if(resonance && resonance!=rho[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,rho[ix],
		      ires+1,omega,ires+1,iloc+1,
		      ires+2,iloc+2,ires+2,iloc+3));
    // channel with the pions exchanged
    if(icharge==0)
      mode->addChannel((PhaseSpaceChannel(phase),ires,rho[ix],
			ires+1,omega,ires+1,iloc+2,
			ires+2,iloc+1,ires+2,iloc+3));
  }
  // reset the masses and widths of the resonances if needed
  for(unsigned int ix=0;ix<3;++ix) {
    mode->resetIntermediate(rho[ix],rhoMasses_[ix],rhoWidths_[ix]);
  }
  // set up the omega masses and widths
  mode->resetIntermediate(omega,omegaMass_,omegaWidth_);
  return true;
}

// the particles produced by the current
tPDVector TwoPionPhotonSNDCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart = {tPDPtr(),
		       getParticleData(ParticleID::pi0),
		       getParticleData(ParticleID::gamma)};
  if(imode==0) {
    if(icharge==3)       extpart[0] = getParticleData(ParticleID::piplus );
    else if(icharge==-3) extpart[0] = getParticleData(ParticleID::piminus);
  }
  else {
    extpart[0] = getParticleData(ParticleID::pi0);
  }
  return extpart;
}

void TwoPionPhotonSNDCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) ++ix;
    temp[ix] = HelicityFunctions::polarizationVector(-decay[2]->momentum()
						     ,ix,Helicity::outgoing);
  }
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
  VectorWaveFunction::constructSpinInfo(temp,decay[2],
					outgoing,true,true);
}

// the hadronic currents
vector<LorentzPolarizationVectorE>
TwoPionPhotonSNDCurrent::current(tcPDPtr resonance,
			      FlavourInfo flavour,
			      const int imode, const int ichan, Energy & scale,
			      const tPDVector & outgoing,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge()+outgoing[2]->iCharge();
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode == 1))
    return vector<LorentzPolarizationVectorE>();
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode!=1) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return vector<LorentzPolarizationVectorE>();
  useMe();
  // polarization vectors of the photon
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) ++ix;
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[2],ix,Helicity::outgoing);
  }
  // total momentum of the system
  Lorentz5Momentum q(momenta[0]+momenta[1]+momenta[2]);
  // overall hadronic mass
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2(q.m2());
  unsigned int imin=0, imax = wgts_.size();
  if(ichan>0) {
    if(outgoing[0]!=outgoing[1])
      imin = ichan;
    else
      imin = ichan/2;
    imax = imin+1;
  }
  if(resonance) {
    switch(resonance->id()/1000) {
    case 0:
      imin = 0;
      break;
    case 100:
      imin = 1;
      break;
    case 30 :
      imin = 2;
      break;
    default:
      assert(false);
    }
    imax=imin+1;
  }
  vector<LorentzPolarizationVectorE> ret(3);
  // need to include exchange of identical particles for the I_3=0 case
  for(int iorder=0;iorder<2;++iorder) {
    Lorentz5Momentum pout(momenta[2]);
    if(outgoing[0]==outgoing[1]) {
      if(ichan>=0&& ichan%2!=iorder) continue;
    }
    else if(iorder==1) continue;
    // add pion momentum
    if(iorder==0) pout += momenta[1];
    else          pout += momenta[0];
    // mass of the omega
    pout.rescaleMass();
    Energy2 s2(pout.m2());
    // compute the rho width
    Energy2 mr2(sqr(rhoMasses_[0]));
    Energy grho = rhoWidths_[0]*mr2/q2*pow(max(double((q2-4.*sqr(mpi_))/(mr2-4.*sqr(mpi_))),0.),1.5);
    Energy qw = Kinematics::pstarTwoBodyDecay(q.mass(),pout.mass(),mpi_);
    grho += pow<3,1>(qw)*sqr(gRhoOmegaPi_)/12./Constants::pi;
    // compute the prefactor
    complex<InvEnergy4> pre = gRhoOmegaPi_*gGammaOmegaPi_/fRho_*
      Resonance::BreitWignerFW(s2,omegaMass_,omegaWidth_)/sqr(omegaMass_);
    if(imode==0) pre *=sqrt(2.);
    Complex bw(0.);
    for(unsigned int ix=imin;ix<imax;++ix) {
      Energy wid = ix==0 ? grho : rhoWidths_[ix];
      Energy2 mR2 = sqr(rhoMasses_[ix]);
      bw += mR2*wgts_[ix]/(mR2-q2-Complex(0.,1.)*q.mass()*wid);
    }
    pre = pre * bw;
    for(unsigned int ix=0;ix<3;++ix) {
      if(ix==1) continue;
      LorentzVector<complex<Energy2> > v2 = Helicity::epsilon(pout,temp[ix],momenta[2]);
      ret[ix] += pre*scale*Helicity::epsilon(q,v2,pout);
    }
  }
  return ret;
}

bool TwoPionPhotonSNDCurrent::accept(vector<int> id) {
  if(id.size()!=3){return false;}
  unsigned int npiplus(0),npi0(0),ngamma(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::piplus) ++npiplus;
    else if(id[ix]==ParticleID::gamma)  ++ngamma;
    else if(id[ix]==ParticleID::pi0)    ++npi0;
  }
  return (npiplus==1&&ngamma==1&&npi0==1) ||
    (npi0==2&&ngamma==1);
}

unsigned int TwoPionPhotonSNDCurrent::decayMode(vector<int> id) {
  int npip(0),npim(0),npi0(0),ngamma(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)         ++npip;
    else if(id[ix]==ParticleID::piminus)   ++npim;
    else if(id[ix]==ParticleID::pi0)       ++npi0;
    else if(id[ix]==ParticleID::gamma)   ++ngamma;
  }
  if((npip==1 || npim == 1) && npi0==1 && ngamma==1)
    return 0;
  else
    return 1;
}

// output the information for the database
void TwoPionPhotonSNDCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::TwoPionPhotonSNDCurrent " << name()
		    << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<3) output << "newdef " << name() << ":RhoMasses " << ix
		    << " " << rhoMasses_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":RhoMasses " << ix
		    << " " << rhoMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<3) output << "newdef " << name() << ":RhoWidths " << ix
		    << " " << rhoWidths_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":RhoWidths " << ix
		    << " " << rhoWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    if(ix<3) output << "newdef " << name() << ":Amplitudes " << ix
		    << " " << amp_[ix] << "\n";
    else     output << "insert " << name() << ":Amplitudes " << ix
		    << " " << amp_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<phase_.size();++ix) {
    if(ix<3) output << "newdef " << name() << ":Phases " << ix
		    << " " << phase_[ix] << "\n";
    else     output << "insert " << name() << ":Phases " << ix
		    << " " << phase_[ix] << "\n";
  }
  output << "newdef " << name() << ":fRho "    << fRho_ << "\n";
  output << "newdef " << name() << ":gRhoOmegaPi "    << gRhoOmegaPi_*GeV << "\n";
  output << "newdef " << name() << ":gGammaOmegaPi "    << gGammaOmegaPi_*GeV << "\n";
  output << "newdef " << name() << ":OmegaMass "    << omegaMass_/GeV << "\n";
  output << "newdef " << name() << ":OmegaWidth "    << omegaWidth_/GeV << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\""
		    << fullName() << "\";" << endl;
}
#line 1 "./OmegaPionSNDCurrent.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OmegaPionSNDCurrent class.
//

#include "OmegaPionSNDCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"

using namespace Herwig;

OmegaPionSNDCurrent::OmegaPionSNDCurrent() {
  // modes handled
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
  // amplitudes for the weights in the current
  amp_   = {1.,0.175,0.014};
  phase_ = {0.,124.,-63.};
  // rho masses and widths
  rhoMasses_ = {0.77526*GeV,1.510*GeV,1.720*GeV};
  rhoWidths_ = {0.1491 *GeV,0.44 *GeV,0.25 *GeV};
  // coupling
  gRhoOmegaPi_   = 15.9/GeV;
  //fRho_        = 4.9583;
  fRho_ = 5.06325;//evaluated with alphaEM in considered energy range
}

IBPtr OmegaPionSNDCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr OmegaPionSNDCurrent::fullclone() const {
  return new_ptr(*this);
}

void OmegaPionSNDCurrent::doinit() {
  WeakCurrent::doinit();
  assert(phase_.size()==amp_.size());
  wgts_.clear();
  Complex ii(0.,1.);
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    double phi = phase_[ix]/180.*Constants::pi;
    wgts_.push_back(amp_[ix]*(cos(phi)+ii*sin(phi)));
  }
  mpi_ = getParticleData(ParticleID::piplus)->mass();
}

void OmegaPionSNDCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV)
     << amp_ << phase_ << wgts_ << fRho_
     << ounit(gRhoOmegaPi_,1./GeV) << ounit(mpi_,GeV);
}

void OmegaPionSNDCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV)
     >> amp_ >> phase_ >> wgts_ >> fRho_
     >> iunit(gRhoOmegaPi_,1./GeV) >> iunit(mpi_,GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<OmegaPionSNDCurrent,WeakCurrent>
describeHerwigOmegaPionSNDCurrent("Herwig::OmegaPionSNDCurrent",
				  "HwWeakCurrents.so");

void OmegaPionSNDCurrent::Init() {

  static ClassDocumentation<OmegaPionSNDCurrent> documentation
    ("The OmegaPionSNDCurrent class provides a current for omega pi"
     " using the model of SND",
     "The current based on \\cite{Achasov:2016zvn} for $\\omega\\pi$ was used.\n",
     "\\bibitem{Achasov:2016zvn}"
     "M.~N.~Achasov {\\it et al.},\n"
     "%``Updated measurement of the $e^+e^- \\to \\omega \\pi^0 \\to \\pi^0\\pi^0\\gamma$ cross section with the SND detector,''\n"
     "Phys.\\ Rev.\\ D {\\bf 94} (2016) no.11,  112001\n"
     "doi:10.1103/PhysRevD.94.112001\n"
     "[arXiv:1610.00235 [hep-ex]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.94.112001;%%\n"
     "%12 citations counted in INSPIRE as of 22 Aug 2018\n");

  static ParVector<OmegaPionSNDCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the rho mesons",
     &OmegaPionSNDCurrent::rhoMasses_, GeV, -1, 775.26*MeV,
     0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<OmegaPionSNDCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the rho mesons",
     &OmegaPionSNDCurrent::rhoWidths_, GeV, -1, 0.1491*GeV,
     0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<OmegaPionSNDCurrent,double> interfaceAmplitudes
    ("Amplitudes",
     "THe amplitudes for the different rho resonances",
     &OmegaPionSNDCurrent::amp_, -1, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static ParVector<OmegaPionSNDCurrent,double> interfacePhase
    ("Phase",
     "The phases for the different rho resonances in degrees",
     &OmegaPionSNDCurrent::phase_, -1, 0., -360., 360.,
     false, false, Interface::limited);

  static Parameter<OmegaPionSNDCurrent,double> interfacefRho
    ("fRho",
     "The coupling of the photon and the rho meson",
     &OmegaPionSNDCurrent::fRho_, 4.9583, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<OmegaPionSNDCurrent,InvEnergy> interfacegRhoOmegaPi
    ("gRhoOmegaPi",
     "The coupling rho-omega-pi",
     &OmegaPionSNDCurrent::gRhoOmegaPi_, 1./GeV,
     15.9/GeV, 0./GeV, 1000./GeV,
     false, false, Interface::limited);

}

// complete the construction of the decay mode for integration
bool OmegaPionSNDCurrent::createMode(int icharge, tcPDPtr resonance,
					 FlavourInfo flavour,
					 unsigned int imode,PhaseSpaceModePtr mode,
					 unsigned int iloc,int ires,
					 PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode >= 1))
    return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown || flavour.I3!=IsoSpin::I3Unknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
    // and check I_3
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode!=1) return false;
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  // other flavour stuff
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  != Beauty::Zero     ) return false;
  // check that the mode is are kinematical allowed
  Energy min = getParticleData(ParticleID::omega)->massMin();
  if(imode==0)
    min += getParticleData(ParticleID::piplus)->mass();
  else
    min += getParticleData(ParticleID::pi0   )->mass();
  if(min>upp) return false;
  // set up the integration channels;
  vector<tPDPtr> rho;
  if(icharge==-3)
    rho = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  else if(icharge==0)
    rho = {getParticleData( 113),getParticleData( 100113),getParticleData( 30113)};
  else if(icharge==3)
    rho = {getParticleData( 213),getParticleData( 100213),getParticleData( 30213)};
  for(unsigned int ix=0;ix<3;++ix) {
    if(resonance && resonance!=rho[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,rho[ix],
		      ires+1,iloc+1,ires+1,iloc+2));
  }
  // reset the masses and widths of the resonances if needed
  for(unsigned int ix=0;ix<3;++ix) {
    mode->resetIntermediate(rho[ix],rhoMasses_[ix],rhoWidths_[ix]);
  }
  return true;
}

// the particles produced by the current
tPDVector OmegaPionSNDCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart = {tPDPtr(),
		       getParticleData(ParticleID::omega)};
  if(imode==0) {
    if(icharge==3)       extpart[0] = getParticleData(ParticleID::piplus );
    else if(icharge==-3) extpart[0] = getParticleData(ParticleID::piminus);
  }
  else {
    extpart[0] = getParticleData(ParticleID::pi0);
  }
  return extpart;
}

void OmegaPionSNDCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    temp[ix] = HelicityFunctions::polarizationVector(-decay[1]->momentum()
						     ,ix,Helicity::outgoing);
  }
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  VectorWaveFunction::constructSpinInfo(temp,decay[1],
					outgoing,true,true);
}

// the hadronic currents
vector<LorentzPolarizationVectorE>
OmegaPionSNDCurrent::current(tcPDPtr resonance,
			      FlavourInfo flavour,
			      const int imode, const int ichan, Energy & scale,
			      const tPDVector & outgoing,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode == 1))
    return vector<LorentzPolarizationVectorE>();
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown || flavour.I3!=IsoSpin::I3Unknown) {
    if(flavour.I!=IsoSpin::IOne) return vector<LorentzPolarizationVectorE>();
    // and check I_3
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode!=1) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  // other flavour stuff
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  != Beauty::Zero     ) return vector<LorentzPolarizationVectorE>();
  useMe();
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix)
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  // locate the particles
  Lorentz5Momentum q(momenta[0]+momenta[1]);
  // overall hadronic mass
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2(q.m2());
  // compute the rho width
  Energy2 mr2(sqr(rhoMasses_[0]));
  Energy grho = rhoWidths_[0]*mr2/q2*pow((q2-4.*sqr(mpi_))/(mr2-4.*sqr(mpi_)),1.5);
  Energy qw = Kinematics::pstarTwoBodyDecay(q.mass(),momenta[0].mass(),momenta[1].mass());
  grho += pow<3,1>(qw)*sqr(gRhoOmegaPi_)/12./Constants::pi;
  unsigned int imin=0, imax = wgts_.size();
  if(ichan>0) {
    imin = ichan;
    imax = ichan+1;
  }
  if(resonance) {
    switch(resonance->id()/1000) {
    case 0:
      imin = 0;
      break;
    case 100:
      imin = 1;
      break;
    case 30 :
      imin = 2;
      break;
    default:
      assert(false);
    }
    imax=imin+1;
  }
  // compute the prefactor
  complex<InvEnergy> pre = gRhoOmegaPi_/fRho_;
  Complex bw(0.);
  for(unsigned int ix=imin;ix<imax;++ix) {
    Energy wid = ix==0 ? grho : rhoWidths_[ix];
    Energy2 mR2 = sqr(rhoMasses_[ix]);
    bw += mR2*wgts_[ix]/(mR2-q2-Complex(0.,1.)*q.mass()*wid);
  }
  pre = pre * bw;
  vector<LorentzPolarizationVectorE> ret(3);
  for(unsigned int ix=0;ix<3;++ix) {
    ret[ix] = pre*Helicity::epsilon(q,temp[ix],momenta[1]);
  }
  if(imode==0) pre *=sqrt(2.);
  return ret;
}

bool OmegaPionSNDCurrent::accept(vector<int> id) {
  if(id.size()!=2){return false;}
  unsigned int npiplus(0),npi0(0),nomega(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::piplus) ++npiplus;
    else if(id[ix]==ParticleID::omega)  ++nomega;
    else if(id[ix]==ParticleID::pi0)    ++npi0;
  }
  return nomega==1 && (npiplus==1||npi0==1);
}

unsigned int OmegaPionSNDCurrent::decayMode(vector<int> id) {
  int npip(0),npim(0),npi0(0),nomega(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)         ++npip;
    else if(id[ix]==ParticleID::piminus)   ++npim;
    else if(id[ix]==ParticleID::pi0)       ++npi0;
    else if(id[ix]==ParticleID::omega)   ++nomega;
  }
  if((npip==1 || npim == 1) && nomega==1)
    return 0;
  else
    return 1;
}

// output the information for the database
void OmegaPionSNDCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::OmegaPionSNDCurrent " << name()
		    << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<3) output << "newdef " << name() << ":RhoMasses " << ix
		    << " " << rhoMasses_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":RhoMasses " << ix
		    << " " << rhoMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<3) output << "newdef " << name() << ":RhoWidths " << ix
		    << " " << rhoWidths_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":RhoWidths " << ix
		    << " " << rhoWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    if(ix<3) output << "newdef " << name() << ":Amplitudes " << ix
		    << " " << amp_[ix] << "\n";
    else     output << "insert " << name() << ":Amplitudes " << ix
		    << " " << amp_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<phase_.size();++ix) {
    if(ix<3) output << "newdef " << name() << ":Phases " << ix
		    << " " << phase_[ix] << "\n";
    else     output << "insert " << name() << ":Phases " << ix
		    << " " << phase_[ix] << "\n";
  }
  output << "newdef " << name() << ":fRho "    << fRho_ << "\n";
  output << "newdef " << name() << ":gRhoOmegaPi "    << gRhoOmegaPi_*GeV << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\""
		    << fullName() << "\";" << endl;
}
#line 1 "./PhiPiCurrent.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PhiPiCurrent class.
//

#include "PhiPiCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

PhiPiCurrent::PhiPiCurrent() {
  // modes handled
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
  // amplitudes for the weights in the current
  amp_   = {0.194/GeV,0.0214/GeV,0./GeV};
  phase_ = {0.,121.,0.};
  br4pi_     = {0.,0.33,0.};
  // rho masses and widths
  rhoMasses_ = {0.77526*GeV,1.593*GeV,1.909*GeV};
  rhoWidths_ = {0.1491 *GeV,0.203*GeV,0.048*GeV};
}

IBPtr PhiPiCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr PhiPiCurrent::fullclone() const {
  return new_ptr(*this);
}

void PhiPiCurrent::doinit() {
  WeakCurrent::doinit();
  assert(phase_.size()==amp_.size());
  wgts_.clear();
  Complex ii(0.,1.);
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    double phi = phase_[ix]/180.*Constants::pi;
    wgts_.push_back(amp_[ix]*(cos(phi)+ii*sin(phi)));
  }
  mpi_ = getParticleData(ParticleID::piplus)->mass();
}

void PhiPiCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV)
     << ounit(amp_,1./GeV) << phase_ << ounit(wgts_,1./GeV)
     << ounit(mpi_,GeV) << br4pi_;
}

void PhiPiCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV)
     >> iunit(amp_,1./GeV) >> phase_ >> iunit(wgts_,1./GeV)
     >> iunit(mpi_,GeV) >> br4pi_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PhiPiCurrent,WeakCurrent>
describeHerwigPhiPiCurrent("Herwig::PhiPiCurrent",
			   "HwWeakCurrents.so");

void PhiPiCurrent::Init() {

  static ClassDocumentation<PhiPiCurrent> documentation
    ("The PhiPiCurrent class implements a model of the current for phi pi"
     "based on the model of Phys.Rev. D77 (2008) 092002, 2008.",
     "The current for $\\phi\\pi$ based on \\cite{Aubert:2007ym} was used.",
     "\\bibitem{Aubert:2007ym}\n"
     "B.~Aubert {\\it et al.} [BaBar Collaboration],\n"
     "%``Measurements of $e^{+} e^{-} \\to K^{+} K^{-} \\eta$,"
     " $K^{+} K^{-} \\pi^0$ and $K^0_{s} K^\\pm \\pi^\\mp$ "
     "cross sections using initial state radiation events,''\n"
     "Phys.\\ Rev.\\ D {\\bf 77} (2008) 092002\n"
     "doi:10.1103/PhysRevD.77.092002\n"
     "[arXiv:0710.4451 [hep-ex]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.77.092002;%%\n"
     "%153 citations counted in INSPIRE as of 27 Aug 2018\n");

  static ParVector<PhiPiCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the rho mesons",
     &PhiPiCurrent::rhoMasses_, GeV, 3, 775.26*MeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<PhiPiCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the rho mesons",
     &PhiPiCurrent::rhoWidths_, GeV, 3, 149.1*MeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<PhiPiCurrent,InvEnergy> interfaceAmplitudes
    ("Amplitudes",
     "The amplitudes for the different resonances",
     &PhiPiCurrent::amp_, 1./GeV, 3, 0./GeV, 0./GeV, 100./GeV,
     false, false, Interface::limited);

  static ParVector<PhiPiCurrent,double> interfacePhase
    ("Phase",
     "The phases for the different rho resonances in degrees",
     &PhiPiCurrent::phase_, 3, 0., 0.0, 360.,
     false, false, Interface::limited);
  
  static ParVector<PhiPiCurrent,double> interfaceBR4Pi
    ("BR4Pi",
     "The branching ratios to 4 pi for the various resonances",
     &PhiPiCurrent::br4pi_, 3, 0., 0.0, 1.0,
     false, false, Interface::limited);


}

// complete the construction of the decay mode for integration
bool PhiPiCurrent::createMode(int icharge, tcPDPtr resonance,
			      FlavourInfo flavour,
			      unsigned int imode,PhaseSpaceModePtr mode,
			      unsigned int iloc,int ires,
			      PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode >= 1))
    return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode!=1) return false;
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  // other flavours
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
  // check that the mode is are kinematical allowed
  Energy min = getParticleData(ParticleID::phi)->massMin();
  if(imode==0)
    min += getParticleData(ParticleID::piplus)->mass();
  else
    min += getParticleData(ParticleID::pi0   )->mass();
  if(min>upp) return false;
  // set up the integration channels;
  vector<tPDPtr> rho;
  if(icharge==-3)
    rho = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  else if(icharge==0)
    rho = {getParticleData( 113),getParticleData( 100113),getParticleData( 30113)};
  else if(icharge==3)
    rho = {getParticleData( 213),getParticleData( 100213),getParticleData( 30213)};
  for(unsigned int ix=0;ix<3;++ix) {
    if(resonance && resonance!=rho[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,rho[ix],
		      ires+1,iloc+1,ires+1,iloc+2));
  }
  // reset the masses and widths of the resonances if needed
  for(unsigned int ix=0;ix<3;++ix) {
    mode->resetIntermediate(rho[ix],rhoMasses_[ix],rhoWidths_[ix]);
  }
  return true;
}

// the particles produced by the current
tPDVector PhiPiCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart = {tPDPtr(),
		       getParticleData(ParticleID::phi)};
  if(imode==0) {
    if(icharge==3)       extpart[0] = getParticleData(ParticleID::piplus );
    else if(icharge==-3) extpart[0] = getParticleData(ParticleID::piminus);
  }
  else {
    extpart[0] = getParticleData(ParticleID::pi0);
  }
  return extpart;
}

void PhiPiCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    temp[ix] = HelicityFunctions::polarizationVector(-decay[1]->momentum()
						     ,ix,Helicity::outgoing);
  }
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  VectorWaveFunction::constructSpinInfo(temp,decay[1],
					outgoing,true,true);
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
PhiPiCurrent::current(tcPDPtr resonance,
		      FlavourInfo flavour,
		      const int imode, const int ichan, Energy & scale, 
		      const tPDVector & outgoing,
		      const vector<Lorentz5Momentum> & momenta,
		      DecayIntegrator::MEOption) const {
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode == 1))
    return vector<LorentzPolarizationVectorE>();
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode!=1) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  // other flavours
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      )
    return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      )
    return vector<LorentzPolarizationVectorE>();
  useMe();
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix)
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  // locate the particles
  Lorentz5Momentum q(momenta[0]+momenta[1]);
  // overall hadronic mass
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2(q.m2());
  // work out the channel
  unsigned int imin=0, imax = wgts_.size();
  if(ichan>0) {
    imin = ichan;
    imax = ichan+1;
  }
  if(resonance) {
    switch(resonance->id()/1000) {
    case 0:
      imin = 0;
      break;
    case 100:
      imin = 1;
      break;
    case 30 :
      imin = 2;
      break;
    default:
      assert(false);
    }
    imax=imin+1;
  }
  complex<InvEnergy> pre(ZERO);
  for(unsigned int ix=imin;ix<imax;++ix) {
    Energy2 mR2 = sqr(rhoMasses_[ix]);
    Energy wid = rhoWidths_[ix]*
      (1.-br4pi_[ix]+ br4pi_[ix]*mR2/q2*pow((q2-16.*sqr(mpi_))/(mR2-16.*sqr(mpi_)),1.5));
    pre += wgts_[ix]*mR2/(mR2-q2-Complex(0.,1.)*q.mass()*wid);
  }
  vector<LorentzPolarizationVectorE> ret(3);
  if(imode==0) pre *= sqrt(2.);
  for(unsigned int ix=0;ix<3;++ix) {
    ret[ix] = pre*Helicity::epsilon(q,temp[ix],momenta[1]);
  }
  return ret;
}

bool PhiPiCurrent::accept(vector<int> id) {
  if(id.size()!=2){return false;}
  unsigned int npiplus(0),npi0(0),nphi(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::piplus) ++npiplus;
    else if(id[ix]==ParticleID::phi)  ++nphi;
    else if(id[ix]==ParticleID::pi0)    ++npi0;
  }
  return nphi==1 && (npiplus==1||npi0==1);
}

unsigned int PhiPiCurrent::decayMode(vector<int> id) {
  int npip(0),npim(0),npi0(0),nphi(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)         ++npip;
    else if(id[ix]==ParticleID::piminus)   ++npim;
    else if(id[ix]==ParticleID::pi0)       ++npi0;
    else if(id[ix]==ParticleID::phi)   ++nphi;
  }
  if((npip==1 || npim == 1) && nphi==1)
    return 0;
  else
    return 1;
}

// output the information for the database
void PhiPiCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::PhiPiCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<rhoMasses_.size();++ix) {
    output << "newdef " << name() << ":RhoMasses " << ix 
	   << " " << rhoMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoWidths_.size();++ix) {
    output << "newdef " << name() << ":RhoWidths " << ix 
	   << " " << rhoWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    output << "newdef " << name() << ":Amplitudes " << ix 
	   << " " << amp_[ix]*GeV << "\n";
  }
  for(unsigned int ix=0;ix<phase_.size();++ix) {
    output << "newdef " << name() << ":Phases " << ix 
	   << " " << phase_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<phase_.size();++ix) {
    output << "newdef " << name() << ":BR4Pi " << ix 
	   << " " << br4pi_[ix] << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./VectorMesonCurrent.cc"
// -*- C++ -*-
//
// VectorMesonCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonCurrent class.
//

#include "VectorMesonCurrent.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMesonCurrent::doinit() {
  unsigned int isize=numberOfModes();
  if(_id.size()!=isize||_decay_constant.size()!=isize)
    {throw InitException() << "Inconsistent parameters in VectorMesonCurrent::doinit()"
			   << Exception::abortnow;}
  WeakCurrent::doinit();
}
VectorMesonCurrent::VectorMesonCurrent()  {
  _id.push_back(213);_decay_constant.push_back(0.1764*GeV2);
  addDecayMode(2,-1);
  _id.push_back(113);_decay_constant.push_back(0.1764*GeV2);
  addDecayMode(1,-1);
  _id.push_back(113);_decay_constant.push_back(0.1764*GeV2);
  addDecayMode(2,-2);
  _id.push_back(223);_decay_constant.push_back(0.1764*GeV2);
  addDecayMode(1,-1);
  _id.push_back(223);_decay_constant.push_back(0.1764*GeV2);
  addDecayMode(2,-2);
  _id.push_back(333);_decay_constant.push_back(0.2380*GeV2);
  addDecayMode(3,-3);
  _id.push_back(313);_decay_constant.push_back(0.2019*GeV2);
  addDecayMode(1,-3);
  _id.push_back(323);_decay_constant.push_back(0.2019*GeV2);
  addDecayMode(2,-3);
  _id.push_back(20213);_decay_constant.push_back(0.4626*GeV2);
  addDecayMode(2,-1);
  _id.push_back(20113);_decay_constant.push_back(0.4626*GeV2);
  addDecayMode(1,-1);
  _id.push_back(20113);_decay_constant.push_back(0.4626*GeV2);
  addDecayMode(2,-2);
  _id.push_back(413);_decay_constant.push_back(0.402*GeV2);
  addDecayMode(4,-1);
  _id.push_back(423);_decay_constant.push_back(0.402*GeV2);
  addDecayMode(4,-2);
  _id.push_back(433);_decay_constant.push_back(0.509*GeV2);
  addDecayMode(4,-3);
  _id.push_back(443);_decay_constant.push_back(1.223*GeV2);
  addDecayMode(4,-4);
  _id.push_back(100443);_decay_constant.push_back(1.08*GeV2);
  addDecayMode(4,-4);
  _id.push_back(20433);_decay_constant.push_back(0.397*GeV2);
  addDecayMode(4,-3);
  // initial size of the vectors
  _initsize=_id.size();
  setInitialModes(_initsize);
}

void VectorMesonCurrent::persistentOutput(PersistentOStream & os) const {
  os << _id << ounit(_decay_constant,GeV2);
}

void VectorMesonCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _id >> iunit(_decay_constant,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorMesonCurrent,WeakCurrent>
describeHerwigVectorMesonCurrent("Herwig::VectorMesonCurrent", "HwWeakCurrents.so");

void VectorMesonCurrent::Init() {

  static ClassDocumentation<VectorMesonCurrent> documentation
    ("The VectorMesonCurrent class implements the current"
     " for the decay of the weak current into a (pseudo)vector meson.");

  static ParVector<VectorMesonCurrent,int> interfaceID
    ("ID",
     "The PDG code for the outgoing meson.",
     &VectorMesonCurrent::_id,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<VectorMesonCurrent,Energy2> interfaceDecay_Constant
    ("Decay_Constant",
     "The decay constant for the meson.",
     &VectorMesonCurrent::_decay_constant, GeV2, -1, 1.0*GeV2,-10.0*GeV2, 10.0*GeV2,
     false, false, true);

}

// create the decay phase space mode
bool VectorMesonCurrent::createMode(int icharge, tcPDPtr resonance,
				    FlavourInfo flavour,
				    unsigned int imode,PhaseSpaceModePtr mode,
				    unsigned int iloc,int ires,
				    PhaseSpaceChannel phase, Energy upp ) {
  assert(!resonance);
  assert(flavour.I==IsoSpin::IUnknown && flavour.I3==IsoSpin::I3Unknown);
  tPDPtr part(getParticleData(_id[imode]));
  // check the mode has the correct charge
  if(abs(icharge)!=abs(int(getParticleData(_id[imode])->iCharge()))) return false;
  // check if the particle is kinematically allowed
  Energy min=part->massMin();
  if(min>upp) return false;
  // construct the mode
  mode->addChannel((PhaseSpaceChannel(phase),ires,iloc+1));
  return true;
}

// outgoing particles 
tPDVector VectorMesonCurrent::particles(int icharge, unsigned int imode, int iq, int ia) {
  tPDPtr part(getParticleData(_id[imode]));
  tPDVector output;
  if(icharge==int(part->iCharge())) {
    if(icharge==0) {
      int iqb,iab;
      decayModeInfo(imode,iqb,iab);
      if(iq==iqb&&ia==iab) output.push_back(part);
      else                 output.push_back(part->CC());
    }
    else output.push_back(part);
  }
  else if(icharge==-int(part->iCharge())) output.push_back(part->CC());
  return output;
}

void VectorMesonCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp;
  VectorWaveFunction::
    calculateWaveFunctions(temp,decay[0],outgoing,false);
  VectorWaveFunction::constructSpinInfo(temp,decay[0],
					outgoing,true,false);
}

vector<LorentzPolarizationVectorE> 
VectorMesonCurrent::current(tcPDPtr resonance,
			    FlavourInfo flavour,
			    const int imode, const int , Energy & scale, 
			    const tPDVector & outgoing,
			    const vector<Lorentz5Momentum> & momenta,
			    DecayIntegrator::MEOption) const {
  assert(!resonance);
  assert(flavour.I==IsoSpin::IUnknown && flavour.I3==IsoSpin::I3Unknown);
  // set up the spin information for the particle and calculate the wavefunctions
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[0],ix,Helicity::outgoing);
  }
  scale=momenta[0].mass();
  // polarization vector
  Energy fact(_decay_constant[imode]/scale);
  // quarks in the current
  int iq,ia;
  decayModeInfo(imode,iq,ia);
  if(abs(iq)==abs(ia)&&abs(iq)<3) {
    fact *= sqrt(0.5);
    if(outgoing[0]->id()==ParticleID::rho0&&abs(iq)==1) fact=-fact;
  }
  // normalise the current
  vector<LorentzPolarizationVectorE> returnval(3);
  for(unsigned int ix=0;ix<3;++ix) {
    returnval[ix] = temp[ix] * fact;
  }
  // return the answer
  return returnval;
}

bool VectorMesonCurrent::accept(vector<int> id) {
  if(id.size()!=1) return false;
  int idtemp(abs(id[0]));
  for(unsigned int ix=0;ix<_id.size();++ix) {
    if(abs(_id[ix])==idtemp) return true;
  }
  return false;
}

unsigned int VectorMesonCurrent::decayMode(vector<int> idout) {
  int idtemp(abs(idout[0])); unsigned int ix(0);
  bool found(false);
  do {
    if(idtemp==abs(_id[ix])) found=true;
    else                     ++ix;
  }
  while(!found);
  return ix;
}

void VectorMesonCurrent::dataBaseOutput(ofstream & output,
					bool header,bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::VectorMesonCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<_id.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":ID " << ix 
	     << " " << _id[ix] << "\n";
      output << "newdef " << name() << ":Decay_Constant " << ix 
	     << " " << _decay_constant[ix]/GeV2 << "\n";
    }
    else {
      output << "insert " << name() << ":ID " << ix 
	     << " " << _id[ix] << "\n";
      output << "insert " << name() << ":Decay_Constant " << ix 
	     << " " << _decay_constant[ix]/GeV2 << "\n";
    }
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./FivePionCurrent.cc"
// -*- C++ -*-
//
// FivePionCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FivePionCurrent class.
//

#include "FivePionCurrent.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

FivePionCurrent::FivePionCurrent() {
  // set the number of modes
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  setInitialModes(3);
  // masses of the intermediates
  _rhomass    = 776*MeV;
  _a1mass     = 1260*MeV;
  _omegamass  = 782*MeV;
  _sigmamass  = 800*MeV;
  // widths of the intermediates
  _rhowidth   = 150*MeV;
  _a1width    = 400*MeV;
  _omegawidth = 8.5*MeV;
  _sigmawidth = 600*MeV;
  // use local values of the resonance parameters
  _localparameters=true;
  // include the rho Breit-Wigners in omega decay
  _rhoomega = true;
  // Normalisation parameters for the different currents
  _c =4.*GeV2;
  _c0=3.;
  // various meson coupling constants
  _fomegarhopi=0.07/MeV;
  _grhopipi=6.0;
  _garhopi=6.*GeV;
  _faaf=4.*GeV;
  _ffpipi=5.*GeV;
  _presigma = ZERO;
  _preomega = ZERO;
}

inline void FivePionCurrent::doinit() {
  WeakCurrent::doinit();
  if(!_localparameters) {
    _rhomass    = getParticleData(ParticleID::rhominus)->mass();
    _rhowidth   = getParticleData(ParticleID::rhominus)->width();
    _omegamass  = getParticleData(ParticleID::omega)->mass();
    _omegawidth = getParticleData(ParticleID::omega)->width();
    _sigmamass  = getParticleData(9000221)->mass();
    _sigmawidth = getParticleData(9000221)->width();
    _a1mass    = getParticleData(ParticleID::a_1minus)->mass();
    _a1width   = getParticleData(ParticleID::a_1minus)->width();
  }
  // prefactors
  _presigma =  _c/sqr(sqr(_a1mass)*_sigmamass*_rhomass)*_faaf*_ffpipi*
    _garhopi*_grhopipi;      
  _preomega = _c0*_fomegarhopi*sqr(_grhopipi/(sqr(_rhomass)*_omegamass));
}

void FivePionCurrent::persistentOutput(PersistentOStream & os) const {
  static const InvEnergy7 InvGeV7 = pow<-7,1>(GeV);
  static const InvEnergy3 InvGeV3 = pow<-3,1>(GeV);
  os << ounit(_rhomass,GeV)  << ounit(_a1mass,GeV) << ounit(_omegamass,GeV) 
     << ounit(_sigmamass,GeV) << ounit(_rhowidth,GeV)  
     << ounit(_a1width,GeV) << ounit(_omegawidth,GeV) << ounit(_sigmawidth,GeV) 
     << _localparameters << ounit(_c,GeV2) << _c0 
     << ounit(_fomegarhopi,1/GeV) << _grhopipi << ounit(_garhopi,GeV) 
     << ounit(_faaf,GeV) << ounit(_ffpipi,GeV)
     << ounit(_preomega,InvGeV7) << ounit(_presigma,InvGeV3) << _rhoomega;
}

void FivePionCurrent::persistentInput(PersistentIStream & is, int) {
  static const InvEnergy7 InvGeV7 = pow<-7,1>(GeV);
  static const InvEnergy3 InvGeV3 = pow<-3,1>(GeV);
  is >> iunit(_rhomass,GeV)  >> iunit(_a1mass,GeV) >> iunit(_omegamass,GeV) 
     >> iunit(_sigmamass,GeV) >> iunit(_rhowidth,GeV)  
     >> iunit(_a1width,GeV) >> iunit(_omegawidth,GeV) >> iunit(_sigmawidth,GeV) 
     >> _localparameters >> iunit(_c,GeV2) >> _c0 
     >> iunit(_fomegarhopi,1/GeV) >> _grhopipi >> iunit(_garhopi,GeV) 
     >> iunit(_faaf,GeV) >> iunit(_ffpipi,GeV) 
     >> iunit(_preomega,InvGeV7) >> iunit(_presigma,InvGeV3) >> _rhoomega;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FivePionCurrent,WeakCurrent>
describeHerwigFivePionCurrent("Herwig::FivePionCurrent", "HwWeakCurrents.so");

void FivePionCurrent::Init() {

  static ClassDocumentation<FivePionCurrent> documentation
    ("The FivePionCurrent class implements the model of hep-ph/0602162",
     "The model of \\cite{Kuhn:2006nw} was used for the hadronic five pion current.",
     "\\bibitem{Kuhn:2006nw} J.~H.~Kuhn and Z.~Was, hep-ph/0602162, (2006).");

  static Parameter<FivePionCurrent,Energy> interfaceRhoMass
    ("RhoMass",
     "The mass of the rho meson",
     &FivePionCurrent::_rhomass, MeV, 776*MeV, 500*MeV, 1000*MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceA1Mass
    ("A1Mass",
     "The mass of the a_1 meson",
     &FivePionCurrent::_a1mass, MeV, 1260*MeV, 1000*MeV, 1500*MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceOmegaMass
    ("OmegaMass",
     "The mass of the omega meson",
     &FivePionCurrent::_omegamass, MeV, 782*MeV, 600*MeV, 900*MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceSigmaMass
    ("SigmaMass",
     "The mass of the sigma meson",
     &FivePionCurrent::_sigmamass, MeV, 800*MeV, 400*MeV, 1200*MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceRhoWidth
    ("RhoWidth",
     "The width of the rho meson",
     &FivePionCurrent::_rhowidth, MeV, 150*MeV, 100*MeV, 300*MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceA1Width
    ("A1Width",
     "The width of the a_1 meson",
     &FivePionCurrent::_a1width, MeV, 400*MeV, 100*MeV, 800*MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceOmegaWidth
    ("OmegaWidth",
     "The width of the omega meson",
     &FivePionCurrent::_omegawidth, MeV, 8.5*MeV, 1.0*MeV, 20.0*MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceSigmaWidth
    ("SigmaWidth",
     "The width of the sigma meson",
     &FivePionCurrent::_sigmawidth, MeV, 600*MeV, 100*MeV, 1200*MeV,
     false, false, Interface::limited);

  static Switch<FivePionCurrent,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the meson masses and widths or those from the"
     " ParticleData objects.",
     &FivePionCurrent::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use values from the particle data objects",
     false);

  static Switch<FivePionCurrent,bool> interfaceRhoOmega
    ("RhoOmega",
     "Option for the treatment of the rho Breit-Wigners in the omega decay",
     &FivePionCurrent::_rhoomega, true, false, false);
  static SwitchOption interfaceRhoOmegaInclude
    (interfaceRhoOmega,
     "Yes",
     "Include the rho Breit-Wigners",
     true);
  static SwitchOption interfaceRhoOmegaOmit
    (interfaceRhoOmega,
     "No",
     "Don't include the rho Breit-Wigners",
     false);

  static Parameter<FivePionCurrent,Energy2> interfaceC
    ("C",
     "The normalisation parameter for the a_1 sigma current",
     &FivePionCurrent::_c, GeV2, 4.0*GeV2, 0.1*GeV2, 20.0*GeV2,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,double> interfaceC0
    ("C0",
     "The normalisation constant for the omega-rho current",
     &FivePionCurrent::_c0, 3., 0.1, 10.0,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,InvEnergy> interfacegomegarhopi
    ("fomegarhopi",
     "The coupling of omega-rho-pi",
     &FivePionCurrent::_fomegarhopi, 1./MeV, 0.07/MeV, 0.01/MeV, 0.2/MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,double> interfacegrhopipi
    ("grhopipi",
     "The coupling for rho-pi-pi",
     &FivePionCurrent::_grhopipi, 6.0, 1.0, 20.0,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfacegarhopi
    ("garhopi",
     "The coupling of a-rho-pi",
     &FivePionCurrent::_garhopi, GeV, 6.0*GeV, 0.1*GeV, 20.0*GeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfacefaaf
    ("faaf",
     "The coupling of a-a-f",
     &FivePionCurrent::_faaf, GeV, 4.0*GeV, 0.1*GeV, 20.0*GeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceffpipi
    ("ffpipi",
     "The coupling of f-pi-pi",
     &FivePionCurrent::_ffpipi, GeV, 5.0*GeV, 0.1*GeV, 20.0*GeV,
     false, false, Interface::limited);
}

bool FivePionCurrent::accept(vector<int> id) {
  bool allowed(false);
  // check five products
  if(id.size()!=5){return false;}
  int npiminus=0,npiplus=0,npi0=0;
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID:: piplus){++npiplus;}
    else if(id[ix]==ParticleID::piminus){++npiminus;}
    else if(id[ix]==ParticleID::pi0){++npi0;}
  }
  if(npiplus>npiminus) swap(npiplus,npiminus);
  if(     npiminus==3&&npiplus==2&&npi0==0) allowed=true;
  else if(npiminus==2&&npiplus==1&&npi0==2) allowed=true;
  else if(npiminus==1&&npiplus==0&&npi0==4) allowed=true;
  return allowed;
}

bool FivePionCurrent::createMode(int icharge, tcPDPtr resonance,
				 FlavourInfo flavour,
				 unsigned int imode,PhaseSpaceModePtr mode,
				 unsigned int iloc,int ires,
				 PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(abs(icharge)!=3) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
     return false;
      break;
    case IsoSpin::I3One:
      if(icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
  // check that the modes are kinematical allowed
  Energy min(ZERO);
  // 3 pi- 2pi+
  if(imode==0) {
    min=5.*getParticleData(ParticleID::piplus)->mass();
  }
  // 2pi- pi+ 2pi0
  else if(imode==1) {
    min=3.*getParticleData(ParticleID::piplus)->mass()
      +2.*getParticleData(ParticleID::pi0)->mass();
  }
  // pi- 4pi0
  else {
    min=   getParticleData(ParticleID::piplus)->mass()
      +4.*getParticleData(ParticleID::pi0)->mass();
  }
  if(min>upp) return false;
  // intermediates for the channels
  tPDPtr omega(getParticleData(ParticleID::omega)),rhop,rhom,
    rho0(getParticleData(ParticleID::rho0)),a1m,a10(getParticleData(ParticleID::a_10)),
    sigma(getParticleData(9000221));
  if(icharge==3) {
    rhop = getParticleData(ParticleID::rhominus);
    rhom = getParticleData(ParticleID::rhoplus);
    a1m  = getParticleData(ParticleID::a_1plus);
  }
  else {
    rhop = getParticleData(ParticleID::rhoplus);
    rhom = getParticleData(ParticleID::rhominus);
    a1m  = getParticleData(ParticleID::a_1minus);
  }
  if(resonance && resonance !=a1m) return false;
  // all charged mode
  if(imode==0) {
    if(sigma) {
      // 1st two a_1 sigma channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma ,ires+2,iloc+4,ires+2,iloc+5,
			ires+1,a1m   ,ires+3, rho0 ,ires+3,iloc+1,
			ires+4,iloc+2,ires+4,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma ,ires+2,iloc+4,ires+2,iloc+5,
			ires+1,a1m   ,ires+3,rho0  ,ires+3,iloc+2,
			ires+4,iloc+1,ires+4,iloc+3));
      // 2nd two a_1 sigma channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma ,ires+2,iloc+4,ires+2,iloc+1,
			ires+1,a1m   ,ires+3, rho0 ,ires+3,iloc+5,
			ires+4,iloc+2,ires+4,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma ,ires+2,iloc+4,ires+2,iloc+1,
			ires+1,a1m   ,ires+3,rho0  ,ires+3,iloc+2,
			ires+4,iloc+5,ires+4,iloc+3));
      // 3rd two a_1 sigma channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma ,ires+2,iloc+4,ires+2,iloc+2,
			ires+1,a1m   ,ires+3, rho0 ,ires+3,iloc+1,
			ires+4,iloc+5,ires+4,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma ,ires+2,iloc+4,ires+2,iloc+2,
			ires+1,a1m   ,ires+3,rho0  ,ires+3,iloc+5,
			ires+4,iloc+1,ires+4,iloc+3));
      // 4th two a_1 sigma channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma ,ires+2,iloc+3,ires+2,iloc+5,
			ires+1,a1m   ,ires+3, rho0 ,ires+3,iloc+1,
			ires+4,iloc+2,ires+4,iloc+4));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma ,ires+2,iloc+3,ires+2,iloc+5,
			ires+1,a1m   ,ires+3, rho0 ,ires+3,iloc+2,
			ires+4,iloc+1,ires+4,iloc+4));
      // 5th two a_1 sigma channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma ,ires+2,iloc+3,ires+2,iloc+1,
			ires+1,a1m   ,ires+3, rho0 ,ires+3,iloc+5,
			ires+4,iloc+2,ires+4,iloc+4));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma ,ires+2,iloc+3,ires+2,iloc+1,
			ires+1,a1m   ,ires+3, rho0 ,ires+3,iloc+2,
			ires+4,iloc+5,ires+4,iloc+4));
      // 6th two a_1 sigma channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma ,ires+2,iloc+3,ires+2,iloc+2,
			ires+1,a1m   ,ires+3, rho0 ,ires+3,iloc+1,
			ires+4,iloc+5,ires+4,iloc+4));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma ,ires+2,iloc+3,ires+2,iloc+2,
			ires+1,a1m   ,ires+3, rho0 ,ires+3,iloc+5,
			ires+4,iloc+1,ires+4,iloc+4));
    }
  }
  // 2 neutral mode
  else if(imode==1) {
    // first three omega channels
    mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
		      ires+1, rhom,ires+2,iloc+4,ires+2,iloc+5,
		      ires+1,omega,ires+3,  rho0,ires+3,iloc+3,
		      ires+4,iloc+1,ires+4,iloc+2));
    mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
		      ires+1, rhom,ires+2,iloc+4,ires+2,iloc+5,
		      ires+1,omega,ires+3,  rhop,ires+3,iloc+2,
		      ires+4,iloc+1,ires+4,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
		      ires+1, rhom,ires+2,iloc+4,ires+2,iloc+5,
		      ires+1,omega,ires+3,  rhom,ires+3,iloc+1,
		      ires+4,iloc+2,ires+4,iloc+3));
    // second three omega channels
    mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
		      ires+1, rhom,ires+2,iloc+2,ires+2,iloc+5,
		      ires+1,omega,ires+3,  rho0,ires+3,iloc+3,
		      ires+4,iloc+1,ires+4,iloc+4));
    mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
		      ires+1, rhom,ires+2,iloc+2,ires+2,iloc+5,
		      ires+1,omega,ires+3,  rhop,ires+3,iloc+4,
		      ires+4,iloc+1,ires+4,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
		      ires+1, rhom,ires+2,iloc+2,ires+2,iloc+5,
		      ires+1,omega,ires+3,  rhom,ires+3,iloc+1,
		      ires+4,iloc+4,ires+4,iloc+3));
    // third  three omega channels
    mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
		      ires+1, rhom,ires+2,iloc+4,ires+2,iloc+3,
		      ires+1,omega,ires+3,  rho0,ires+3,iloc+5,
		      ires+4,iloc+1,ires+4,iloc+2));
    mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
		      ires+1, rhom,ires+2,iloc+4,ires+2,iloc+3,
		      ires+1,omega,ires+3,  rhop,ires+3,iloc+2,
		      ires+4,iloc+1,ires+4,iloc+5));
    mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
		      ires+1, rhom,ires+2,iloc+4,ires+2,iloc+3,
		      ires+1,omega,ires+3,  rhom,ires+3,iloc+1,
		      ires+4,iloc+2,ires+4,iloc+5));
    // fourth three omega channels
    mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
		      ires+1, rhom,ires+2,iloc+2,ires+2,iloc+3,
		      ires+1,omega,ires+3,  rho0,ires+3,iloc+5,
		      ires+4,iloc+1,ires+4,iloc+4));
    mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
		      ires+1, rhom,ires+2,iloc+2,ires+2,iloc+3,
		      ires+1,omega,ires+3,  rhop,ires+3,iloc+4,
		      ires+4,iloc+1,ires+4,iloc+5));
    mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
		      ires+1, rhom,ires+2,iloc+2,ires+2,iloc+3,
		      ires+1,omega,ires+3,  rhom,ires+3,iloc+1,
		      ires+4,iloc+4,ires+4,iloc+5));
    if(sigma) {
      // first two sigma a_1 channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma,ires+2,iloc+1,ires+2,iloc+2,
			ires+1,  a1m,ires+3,  rhom,ires+3,iloc+3,
			ires+4,iloc+4,ires+4,iloc+5));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma,ires+2,iloc+1,ires+2,iloc+2,
			ires+1,  a1m,ires+3,  rhom,ires+3,iloc+5,
			ires+4,iloc+4,ires+4,iloc+3));
      // second two sigma a_1 channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma,ires+2,iloc+1,ires+2,iloc+4,
			ires+1,  a1m,ires+3,  rhom,ires+3,iloc+3,
			ires+4,iloc+2,ires+4,iloc+5));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma,ires+2,iloc+1,ires+2,iloc+4,
			ires+1,  a1m,ires+3,  rhom,ires+3,iloc+5,
			ires+4,iloc+2,ires+4,iloc+3));
      // third two sigma a_1 channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma,ires+2,iloc+3,ires+2,iloc+5,
			ires+1,  a1m,ires+3,  rho0,ires+3,iloc+2,
			ires+4,iloc+1,ires+4,iloc+4));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1,sigma,ires+2,iloc+3,ires+2,iloc+5,
			ires+1,  a1m,ires+3,  rho0,ires+3,iloc+4,
			ires+4,iloc+1,ires+4,iloc+2));
    }
  }
  // 4 neutral mode
  else {
    if(sigma) {
      // 1st two sigma a_1 channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1, sigma,ires+2,iloc+4,ires+2,iloc+5,
			ires+1,   a1m,ires+3,  rhom,ires+3,iloc+2,
			ires+4,iloc+1,ires+4,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1, sigma,ires+2,iloc+4,ires+2,iloc+5,
			ires+1,   a1m,ires+3,  rhom,ires+3,iloc+1,
			ires+4,iloc+2,ires+4,iloc+3));
      //   // 2nd two sigma a_1 channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1, sigma,ires+2,iloc+4,ires+2,iloc+1,
			ires+1,   a1m,ires+3,  rhom,ires+3,iloc+2,
			ires+4,iloc+5,ires+4,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1, sigma,ires+2,iloc+4,ires+2,iloc+1,
			ires+1,   a1m,ires+3,  rhom,ires+3,iloc+5,
			ires+4,iloc+2,ires+4,iloc+3));
      // 3rd two sigma a_1 channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1, sigma,ires+2,iloc+1,ires+2,iloc+5,
			ires+1,   a1m,ires+3,  rhom,ires+3,iloc+2,
			ires+4,iloc+4,ires+4,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
			ires+1, sigma,ires+2,iloc+1,ires+2,iloc+5,
			ires+1,   a1m,ires+3,  rhom,ires+3,iloc+4,
			ires+4,iloc+2,ires+4,iloc+3));
      // 4th two sigma a_1 channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
      			ires+1, sigma,ires+2,iloc+2,ires+2,iloc+5,
      			ires+1,   a1m,ires+3,  rhom,ires+3,iloc+4,
      			ires+4,iloc+1,ires+4,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
      			ires+1, sigma,ires+2,iloc+2,ires+2,iloc+5,
      			ires+1,   a1m,ires+3,  rhom,ires+3,iloc+1,
      			ires+4,iloc+4,ires+4,iloc+3));
      // 5th two sigma a_1 channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
      			ires+1, sigma,ires+2,iloc+4,ires+2,iloc+2,
      			ires+1,   a1m,ires+3,  rhom,ires+3,iloc+5,
      			ires+4,iloc+1,ires+4,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
      			ires+1, sigma,ires+2,iloc+4,ires+2,iloc+2,
      			ires+1,   a1m,ires+3,  rhom,ires+3,iloc+1,
      			ires+4,iloc+5,ires+4,iloc+3));
      // 6th two sigma a_1 channels
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
      			ires+1, sigma,ires+2,iloc+1,ires+2,iloc+2,
      			ires+1,   a1m,ires+3,  rhom,ires+3,iloc+4,
      			ires+4,iloc+5,ires+4,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1m,
      			ires+1, sigma,ires+2,iloc+1,ires+2,iloc+2,
      			ires+1,   a1m,ires+3,  rhom,ires+3,iloc+5,
      			ires+4,iloc+4,ires+4,iloc+3));
    }
  }
  // reset the parameters of the resonances if using local values
  if(_localparameters) {
    mode->resetIntermediate(rhom,_rhomass,_rhowidth);
    mode->resetIntermediate(rhop,_rhomass,_rhowidth);
    mode->resetIntermediate(rho0,_rhomass,_rhowidth);
    mode->resetIntermediate(omega,_omegamass,_omegawidth);
    mode->resetIntermediate(a1m,_a1mass,_a1width);
    mode->resetIntermediate(a10,_a1mass,_a1width);
    if(sigma) mode->resetIntermediate(sigma,_sigmamass,_sigmawidth);
  }
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector FivePionCurrent::particles(int icharge, unsigned int imode,int,int) {
  // particle data objects for the pions
  tPDPtr piplus (getParticleData(ParticleID::piplus ));
  tPDPtr pi0    (getParticleData(ParticleID::pi0    ));
  tPDPtr piminus(getParticleData(ParticleID::piminus));
  if(icharge==3) swap(piplus,piminus);
  tPDVector output(5);
  // all charged
  if(imode==0) {
    output[0]=piminus;
    output[1]=piminus;
    output[2]=piplus;;
    output[3]=piplus;
    output[4]=piminus;
  }
  // two neutral
  else if(imode==1) {
    output[0]=piplus;
    output[1]=piminus;
    output[2]=pi0;;
    output[3]=piminus;
    output[4]=pi0;
  }
  // four neutral
  else {
    output[0]=pi0;
    output[1]=pi0;
    output[2]=piminus;;
    output[3]=pi0;
    output[4]=pi0;
  }
  return output;
}

// the decay mode
unsigned int FivePionCurrent::decayMode(vector<int> idout) {
  unsigned int npi(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::pi0) ++npi;
  }
  return npi/2;
}

// output the information for the database
void FivePionCurrent::dataBaseOutput(ofstream & output,bool header,bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::FivePionCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  output << "newdef " << name() << ":RhoMass "    << _rhomass/MeV << "\n";
  output << "newdef " << name() << ":A1Mass  "    << _a1mass/MeV  << "\n";
  output << "newdef " << name() << ":SigmaMass  " << _sigmamass/MeV  << "\n";
  output << "newdef " << name() << ":OmegaMass  " << _omegamass/MeV  << "\n";
  output << "newdef " << name() << ":RhoWidth "    << _rhowidth/MeV << "\n";
  output << "newdef " << name() << ":A1Width  "    << _a1width/MeV  << "\n";
  output << "newdef " << name() << ":SigmaWidth  " << _sigmawidth/MeV  << "\n";
  output << "newdef " << name() << ":OmegaWidth  " << _omegawidth/MeV  << "\n";
  output << "newdef " << name() << ":LocalParameters " <<  _localparameters << "\n";
  output << "newdef " << name() << ":RhoOmega " << _rhoomega << "\n";
  output << "newdef " << name() << ":C " << _c/GeV2 << "\n";
  output << "newdef " << name() << ":C0 " << _c0 << "\n";
  output << "newdef " << name() << ":fomegarhopi " <<_fomegarhopi*MeV << "\n";
  output << "newdef " << name() << ":grhopipi " <<_grhopipi << "\n";
  output << "newdef " << name() << ":garhopi " << _garhopi/GeV << "\n";
  output << "newdef " << name() << ":faaf " <<_faaf/GeV << "\n";
  output << "newdef " << name() << ":ffpipi " << _ffpipi/GeV << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";\n";
}
 
vector<LorentzPolarizationVectorE> 
FivePionCurrent::current(tcPDPtr,
			 FlavourInfo flavour,
			 const int imode, const int ichan,Energy & scale, 
			 const tPDVector & outgoing,
			 const vector<Lorentz5Momentum> & momenta,
			 DecayIntegrator::MEOption) const {
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IOne)
    return vector<LorentzPolarizationVectorE>();
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return vector<LorentzPolarizationVectorE>();
  useMe();
  LorentzVector<complex<InvEnergy2> > output;
  Lorentz5Momentum q1(momenta[0]);
  Lorentz5Momentum q2(momenta[1]);
  Lorentz5Momentum q3(momenta[2]);
  Lorentz5Momentum q4(momenta[3]);
  Lorentz5Momentum q5(momenta[4]);
  // total momentum
  Lorentz5Momentum Q(q1+q2+q3+q4+q5);
  Q.rescaleMass();
  scale=Q.mass();
  // decide which decay mode
  if(imode==0) {
    if(ichan<0) {
      output=
	a1SigmaCurrent(0,Q,q1,q2,q3,q4,q5)+
	a1SigmaCurrent(0,Q,q5,q2,q3,q4,q1)+
	a1SigmaCurrent(0,Q,q1,q5,q3,q4,q2)+
	a1SigmaCurrent(0,Q,q1,q2,q4,q3,q5)+
	a1SigmaCurrent(0,Q,q5,q2,q4,q3,q1)+
	a1SigmaCurrent(0,Q,q1,q5,q4,q3,q2);
    }
    else if(ichan==0 ) output=a1SigmaCurrent(2,Q,q1,q2,q3,q4,q5);
    else if(ichan==1 ) output=a1SigmaCurrent(1,Q,q1,q2,q3,q4,q5);
    else if(ichan==2 ) output=a1SigmaCurrent(2,Q,q5,q2,q3,q4,q1);
    else if(ichan==3 ) output=a1SigmaCurrent(1,Q,q5,q2,q3,q4,q1);
    else if(ichan==4 ) output=a1SigmaCurrent(2,Q,q1,q5,q3,q4,q2);
    else if(ichan==5 ) output=a1SigmaCurrent(1,Q,q1,q5,q3,q4,q2);
    else if(ichan==6 ) output=a1SigmaCurrent(2,Q,q1,q2,q4,q3,q5);
    else if(ichan==7 ) output=a1SigmaCurrent(1,Q,q1,q2,q4,q3,q5);
    else if(ichan==8 ) output=a1SigmaCurrent(2,Q,q5,q2,q4,q3,q1);
    else if(ichan==9 ) output=a1SigmaCurrent(1,Q,q5,q2,q4,q3,q1);
    else if(ichan==10) output=a1SigmaCurrent(2,Q,q1,q5,q4,q3,q2);
    else if(ichan==11) output=a1SigmaCurrent(1,Q,q1,q5,q4,q3,q2);
    // identical particle symmetry factor
    output/=sqrt(12.);
  }
  else if(imode==1) {
    if(ichan<0) {
      output=
	rhoOmegaCurrent(0,Q,q1,q2,q3,q4,q5)
	+rhoOmegaCurrent(0,Q,q1,q4,q3,q2,q5)
	+rhoOmegaCurrent(0,Q,q1,q2,q5,q4,q3)
	+rhoOmegaCurrent(0,Q,q1,q4,q5,q2,q3)
	+a1SigmaCurrent(0,Q,q2,q4,q1,q3,q5)
	+a1SigmaCurrent(0,Q,q3,q5,q2,q1,q4)
	+a1SigmaCurrent(0,Q,q3,q5,q4,q1,q2);
    }
    else if(ichan==0 ) output=rhoOmegaCurrent(3,Q,q1,q2,q3,q4,q5);
    else if(ichan==1 ) output=rhoOmegaCurrent(2,Q,q1,q2,q3,q4,q5);
    else if(ichan==2 ) output=rhoOmegaCurrent(1,Q,q1,q2,q3,q4,q5);
    else if(ichan==3 ) output=rhoOmegaCurrent(3,Q,q1,q4,q3,q2,q5);
    else if(ichan==4 ) output=rhoOmegaCurrent(2,Q,q1,q4,q3,q2,q5);
    else if(ichan==5 ) output=rhoOmegaCurrent(1,Q,q1,q4,q3,q2,q5);
    else if(ichan==6 ) output=rhoOmegaCurrent(3,Q,q1,q2,q5,q4,q3);
    else if(ichan==7 ) output=rhoOmegaCurrent(2,Q,q1,q2,q5,q4,q3);
    else if(ichan==8 ) output=rhoOmegaCurrent(1,Q,q1,q2,q5,q4,q3);
    else if(ichan==9 ) output=rhoOmegaCurrent(3,Q,q1,q4,q5,q2,q3);
    else if(ichan==10) output=rhoOmegaCurrent(2,Q,q1,q4,q5,q2,q3);
    else if(ichan==11) output=rhoOmegaCurrent(1,Q,q1,q4,q5,q2,q3);
    else if(ichan==12) output=a1SigmaCurrent(2,Q,q3,q5,q4,q1,q2);
    else if(ichan==13) output=a1SigmaCurrent(1,Q,q3,q5,q4,q1,q2);
    else if(ichan==14) output=a1SigmaCurrent(2,Q,q3,q5,q2,q1,q4);
    else if(ichan==15) output=a1SigmaCurrent(1,Q,q3,q5,q2,q1,q4);
    else if(ichan==16) output=a1SigmaCurrent(2,Q,q2,q4,q1,q3,q5);
    else if(ichan==17) output=a1SigmaCurrent(1,Q,q2,q4,q1,q3,q5);
    // identical particle symmetry factor
    output/=2.;
  }
  else if(imode==2) {
    if(ichan<0) {
      output=
	a1SigmaCurrent(0,Q,q1,q2,q3,q4,q5)+
	a1SigmaCurrent(0,Q,q5,q2,q3,q4,q1)+
	a1SigmaCurrent(0,Q,q2,q4,q3,q1,q5)+
	a1SigmaCurrent(0,Q,q1,q4,q3,q2,q5)+
	a1SigmaCurrent(0,Q,q1,q5,q3,q4,q2)+
	a1SigmaCurrent(0,Q,q4,q5,q3,q1,q2);
    }
    else if(ichan==0 ) output=a1SigmaCurrent(1,Q,q1,q2,q3,q4,q5);
    else if(ichan==1 ) output=a1SigmaCurrent(2,Q,q1,q2,q3,q4,q5);
    else if(ichan==2 ) output=a1SigmaCurrent(1,Q,q5,q2,q3,q4,q1);
    else if(ichan==3 ) output=a1SigmaCurrent(2,Q,q5,q2,q3,q4,q1);
    else if(ichan==4 ) output=a1SigmaCurrent(2,Q,q2,q4,q3,q1,q5);
    else if(ichan==5 ) output=a1SigmaCurrent(1,Q,q2,q4,q3,q1,q5);
    else if(ichan==6 ) output=a1SigmaCurrent(1,Q,q1,q4,q3,q2,q5);
    else if(ichan==7 ) output=a1SigmaCurrent(2,Q,q1,q4,q3,q2,q5);
    else if(ichan==8 ) output=a1SigmaCurrent(1,Q,q1,q5,q3,q4,q2);
    else if(ichan==9 ) output=a1SigmaCurrent(2,Q,q1,q5,q3,q4,q2);
    else if(ichan==10) output=a1SigmaCurrent(2,Q,q4,q5,q3,q1,q2);
    else if(ichan==11) output=a1SigmaCurrent(1,Q,q4,q5,q3,q1,q2);
    // identical particle symmetry factor
    output/=sqrt(24.);
  }
  else {
    throw Exception() << "Unknown decay mode in the " 
				 << "FivePionCurrent::"
				 << "hadronCurrent()" << Exception::abortnow;
  }
  // normalise and return the current
  return vector<LorentzPolarizationVectorE>(1, output * pow<3,1>(scale));
}
#line 1 "./KPiCurrent.cc"
// -*- C++ -*-
//
// KPiCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KPiCurrent class.
//

#include "KPiCurrent.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::outgoing;

KPiCurrent::KPiCurrent() :
  _localparameters(true),_transverse(false), _cV(1.),_cS(0.2),
  _mpi(ZERO), _mK(ZERO) {
  // set up for the modes in the base class
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  setInitialModes(2);
  // parameters for the vector resonances
  _vecmag  .push_back(1.);_vecmag  .push_back(-0.135);
  _vecphase.push_back(0.);_vecphase.push_back(180.  );
  _vecmass .push_back(891.6*MeV);_vecmass .push_back(1412.*MeV);
  _vecwidth.push_back( 50. *MeV);_vecwidth.push_back( 227.*MeV);
  // parameters for the scalar resonances
  _scamag  .push_back(0.);_scamag  .push_back(1.);
  _scaphase.push_back(0.);_scaphase.push_back(0.);
  _scamass .push_back(841.*MeV);_scamass .push_back(1429.*MeV);
  _scawidth.push_back(618.*MeV);_scawidth.push_back( 287.*MeV);
}

void KPiCurrent::persistentOutput(PersistentOStream & os) const {
  os << _cV << _cS << _localparameters 
     << ounit(_mpi,GeV) << ounit(_mK,GeV) 
     << _resmap
     << _vecmag << _vecphase << _vecwgt 
     << ounit(_vecmass,MeV) << ounit(_vecwidth,MeV)
     << _scamag << _scaphase << _scawgt 
     << ounit(_scamass,MeV) << ounit(_scawidth,MeV)
     << _transverse;
}

void KPiCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _cV >> _cS >> _localparameters 
     >> iunit(_mpi,GeV) >> iunit(_mK,GeV) 
     >> _resmap
     >> _vecmag >> _vecphase >> _vecwgt 
     >> iunit(_vecmass,MeV) >> iunit(_vecwidth,MeV) 
     >> _scamag >> _scaphase >> _scawgt 
     >> iunit(_scamass,MeV) >> iunit(_scawidth,MeV)
     >> _transverse;
}

void KPiCurrent::doinit() {
  WeakCurrent::doinit();
  // check consistency of parametrers
  if(_vecmass.size()!=_vecwidth.size()||
     _scamass.size()!=_scawidth.size()) {
    throw InitException() << "Inconsistent parameters in KPiCurrent"
			  << "doinit()" << Exception::abortnow;
  }
  // the resonances
  tPDPtr vec[3]={getParticleData(-323   ),getParticleData(-100323),
		 getParticleData(-30323 )};
  tPDPtr sca[3]={getParticleData(-9000321),getParticleData(-10321)};
  // reset the masses in the form-factors if needed
  if(_localparameters) {
    if(_vecmass.size()<3) {
      for(unsigned int ix=_vecmass.size();ix<3;++ix) {
	if(vec[ix]) {
	  _vecmass.push_back( vec[ix]->mass() );
	  _vecwidth.push_back(vec[ix]->width());
	}
      }
    }
    if(_scamass.size()<2) {
      for(unsigned int ix=_scamass.size();ix<2;++ix) {
	if(sca[ix]) {
	  _scamass.push_back( sca[ix]->mass() );
	  _scawidth.push_back(sca[ix]->width());
	}
      }
    }
  }
  else {
    _vecmass.clear();_vecwidth.clear();
    for(unsigned int ix=0;ix<3;++ix) {
      if(vec[ix]) {
	_vecmass .push_back(vec[ix]->mass() );
	_vecwidth.push_back(vec[ix]->width());
      }
    }
    _scamass.clear();_scawidth.clear();
    for(unsigned int ix=0;ix<2;++ix) {
      if(sca[ix]) {
	_scamass .push_back(sca[ix]->mass() );
	_scawidth.push_back(sca[ix]->width());
      }
    }
  }
  _mpi=getParticleData(ParticleID::piplus)->mass();
  _mK =getParticleData(ParticleID::K0    )->mass();
  // weight for the vector channels
  if(_vecmag.size()!=_vecphase.size())
    throw InitException() << "The vectors containing the weights and phase for the"
			  << "vector channel must be the same size in"
			  << "KPiCurrent::doinit()" 
			  << Exception::runerror;
  _vecwgt.resize(_vecmag.size());
  for(unsigned int ix=0;ix<_vecwgt.size();++ix) {
    double angle = _vecphase[ix]/180.*Constants::pi;
    _vecwgt[ix] = _vecmag[ix]*(cos(angle)+Complex(0.,1.)*sin(angle));
  }
  // weight for the scalar channels
  if(_scamag.size()!=_scaphase.size())
    throw InitException() << "The vectors containing the weights and phase for the"
			  << "scalar channel must be the same size in"
			  << "KPiCurrent::doinit()" 
			  << Exception::runerror;
  _scawgt.resize(_scamag.size());
  for(unsigned int ix=0;ix<_scawgt.size();++ix) {
    double angle = _scaphase[ix]/180.*Constants::pi;
    _scawgt[ix] = _scamag[ix]*(cos(angle)+Complex(0.,1.)*sin(angle));
  }
  // mapping for the resonaces
  int ires(-1);
  for(unsigned int ix=0;ix<3;++ix) {
    if(vec[ix]) ++ires;
    if(ires<int(_vecwgt.size())) _resmap.push_back(ires);
  }
  if(_resmap.size()<_vecwgt.size()) {
    for(unsigned int ix=_resmap.size();ix<_vecwgt.size();++ix) {
      _resmap.push_back(-1);
    }
  }
  ires=-1;
  for(unsigned int ix=0;ix<2;++ix) {
    if(sca[ix]) ++ires;
    if(ires<int(_scawgt.size())) _resmap.push_back(ires);
  }
  if(_resmap.size()<_vecwgt.size()+_scawgt.size()) {
    for(unsigned int ix=_resmap.size()-_scawgt.size();
	ix<_scawgt.size();++ix) {
      _resmap.push_back(-1);
    }
  }
}
// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KPiCurrent,WeakCurrent>
describeHerwigKPiCurrent("Herwig::KPiCurrent", "HwWeakCurrents.so");

void KPiCurrent::Init() {

  static ClassDocumentation<KPiCurrent> documentation
    ("The KPiCurrent class",
     "The K pi weak current has the form of \\cite{Finkemeier:1996dh}.",
     "%\\cite{Finkemeier:1996dh}\n"
     "\\bibitem{Finkemeier:1996dh}\n"
     "  M.~Finkemeier and E.~Mirkes,\n"
     "  %``The scalar contribution to tau --> K pi nu/tau,''\n"
     "  Z.\\ Phys.\\  C {\\bf 72}, 619 (1996)\n"
     "  [arXiv:hep-ph/9601275].\n"
     "  %%CITATION = ZEPYA,C72,619;%%\n"
     );

  static Parameter<KPiCurrent,double> interfacecV
    ("cV",
     "The weight for the vector contribution",
     &KPiCurrent::_cV, 1., 0., 10.0,
     false, false, Interface::limited);

  static Parameter<KPiCurrent,double> interfacecS
    ("cS",
     "The weight for the scalar contribution",
     &KPiCurrent::_cS, 0.2, -10.0, 10.0,
     false, false, Interface::limited);
  
  static ParVector<KPiCurrent,double> interfaceVectorMagnitude
    ("VectorMagnitude",
     "Magnitude of the weight for the different vector resonances",
     &KPiCurrent::_vecmag, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<KPiCurrent,double> interfaceVectorPhase
    ("VectorPhase",
     "Phase of the weight of the different vector resonances",
     &KPiCurrent::_vecphase, -1, 0., 0, 0,
     false, false, Interface::nolimits);
  
  static ParVector<KPiCurrent,double> interfaceScalarMagnitude
    ("ScalarMagnitude",
     "Magnitude of the weight for the different scalar resonances",
     &KPiCurrent::_scamag, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<KPiCurrent,double> interfaceScalarPhase
    ("ScalarPhase",
     "Phase of the weight of the different scalar resonances",
     &KPiCurrent::_scaphase, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static Switch<KPiCurrent,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values for the masses and widths of the resonances or those"
     " from the ParticleData objects",
     &KPiCurrent::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use the values from the particle data objects",
     false);

  static Switch<KPiCurrent,bool> interfaceTransverse
    ("Transverse",
     "Form of the vector projection operator.",
     &KPiCurrent::_transverse, false, false, false);
  static SwitchOption interfaceTransverseTransverse
    (interfaceTransverse,
     "Transverse",
     "Use 1/Q^2 in the projection operator to force it to be transverse",
     true);
  static SwitchOption interfaceTransverseMass
    (interfaceTransverse,
     "Mass",
     "Use the on-shell mass in the projection operator",
     false);

  static ParVector<KPiCurrent,Energy> interfaceVectorMass
    ("VectorMass",
     "Masses of the vector resonances",
     &KPiCurrent::_vecmass, MeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<KPiCurrent,Energy> interfaceVectorWidth
    ("VectorWidth",
     "Widths of the vector resonances",
     &KPiCurrent::_vecwidth, MeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<KPiCurrent,Energy> interfaceScalarMass
    ("ScalarMass",
     "Masses of the scalar resonances",
     &KPiCurrent::_scamass, MeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<KPiCurrent,Energy> interfaceScalarWidth
    ("ScalarWidth",
     "Widths of the scalar resonances",
     &KPiCurrent::_scawidth, MeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);
}

bool KPiCurrent::accept(vector<int> id) {
  bool allowed(false);
  // check there are only two particles
  if(id.size()!=2){return false;}
  if      ((id[0]==ParticleID::Kminus && id[1]==ParticleID::pi0)    ||
	  (id[0]==ParticleID::pi0    && id[1]==ParticleID::Kminus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::pi0)    ||
	  (id[0]==ParticleID::pi0    && id[1]==ParticleID::Kplus)) allowed=true;
  // single neutral kaon
  else if((id[0]==ParticleID::piminus && id[1]==ParticleID::Kbar0)   ||
	  (id[0]==ParticleID::Kbar0   && id[1]==ParticleID::piminus) ||
	  (id[0]==ParticleID::piplus  && id[1]==ParticleID::K0)      ||
	  (id[0]==ParticleID::K0      && id[1]==ParticleID::piplus)) allowed=true;
  return allowed;
}

tPDVector KPiCurrent::particles(int icharge, unsigned int imode, int,int) {
  if(abs(icharge)!=3) return tPDVector();
  tPDVector output(2);
  if(imode==0) {
    output[0]=getParticleData(ParticleID::Kplus);
    output[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==1) {
    output[0]=getParticleData(ParticleID::K0);
    output[1]=getParticleData(ParticleID::piplus);
  }
  if(icharge==-3) {
    for(unsigned int ix=0;ix<output.size();++ix) {
      if(output[ix]->CC()) output[ix]=output[ix]->CC();
    }
  }
  return output;
}

unsigned int KPiCurrent::decayMode(vector<int> id) {
  unsigned int imode(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::K0) imode=1;
  }
  return imode;
}

bool KPiCurrent::createMode(int icharge, tcPDPtr resonance,
			    FlavourInfo flavour,
			    unsigned int imode,PhaseSpaceModePtr mode,
			    unsigned int iloc,int ires,
			    PhaseSpaceChannel phase, Energy upp ) {
  if(abs(icharge)!=3) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IHalf) return false;
  } 
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Half:
      if(icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusHalf:
      if(icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown) {
    if(icharge== 3 and flavour.strange != Strangeness::PlusOne ) return false;
    if(icharge==-3 and flavour.strange != Strangeness::MinusOne) return false;
  }
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
  // make sure that the decays are kinematically allowed
  tPDPtr part[2];
  if(imode==0) {
    part[0]=getParticleData(ParticleID::Kplus);
    part[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==1) {
    part[0]=getParticleData(ParticleID::K0);
    part[1]=getParticleData(ParticleID::piplus);
  }
  else {
    return false;
  }
  Energy min(part[0]->massMin()+part[1]->massMin());
  if(min>upp) return false;
  // possible resonances
  tPDPtr res[5]={getParticleData(-323   ),getParticleData(-100323),
		 getParticleData(-30323 ),getParticleData(-9000321),
		 getParticleData(-10321)};
  // create the channels
  for(unsigned int ix=0;ix<5;++ix) {
    if(!res[ix]) continue;
    if(resonance && resonance != res[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],
		      ires+1,iloc+1,ires+1,iloc+2));
  }
  // reset the masses in the intergrators if needed
  if(_localparameters) {
    // for the vectors
    for(unsigned int ix=0;ix<3;++ix) {
      if(ix<_vecmass.size()&&res[ix]) {
	mode->resetIntermediate(res[ix],_vecmass[ix],_vecwidth[ix]);
      }
    }
    // for the scalars
    for(unsigned int ix=3;ix<5;++ix) {
      if(ix-3<_scamass.size()&&res[ix]) {
	mode->resetIntermediate(res[ix],_scamass[ix-3],_scawidth[ix-3]);
      }
    }
  }
  return true;
}

void KPiCurrent::dataBaseOutput(ofstream & output,bool header,
				bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::KPiCurrent " << name() 
		    << " HeWeakCurrents.so\n";
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  output << "newdef " << name() << ":Transverse "      << _transverse << "\n";
  output << "newdef " << name() << ":cV " << _cV << "\n";
  output << "newdef " << name() << ":cS " << _cS << "\n";
  for(unsigned int ix=0;ix<_vecmag.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":VectorMagnitude " << ix << " " << _vecmag[ix]   << "\n";
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":VectorPhase "     << ix << " " << _vecphase[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_scamag.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":ScalarMagnitude " << ix << " " << _scamag[ix]  << "\n";
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":ScalarPhase "     << ix << " " << _scaphase[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_vecmass.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":VectorMass "  << ix << " " << _vecmass[ix]/MeV  << "\n";
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":VectorWidth " << ix << " " << _vecwidth[ix]/MeV << "\n";
  }
  for(unsigned int ix=0;ix<_scamass.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":ScalarMass "  << ix << " " << _scamass[ix]/MeV  << "\n";
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":ScalarWidth " << ix << " " << _scawidth[ix]/MeV << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

vector<LorentzPolarizationVectorE> 
KPiCurrent::current(tcPDPtr resonance,
		    FlavourInfo flavour,
		    const int imode, const int ichan,Energy & scale, 
		    const tPDVector & outgoing,
		    const vector<Lorentz5Momentum> & momenta,
		    DecayIntegrator::MEOption) const {
  useMe();
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IHalf)
    return vector<LorentzPolarizationVectorE>();
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Half:
      if(icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusHalf:
      if(icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown) {
    if(icharge== 3 and flavour.strange != Strangeness::PlusOne ) return vector<LorentzPolarizationVectorE>();
    if(icharge==-3 and flavour.strange != Strangeness::MinusOne) return vector<LorentzPolarizationVectorE>();
  }
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return vector<LorentzPolarizationVectorE>();
  // momentum difference and sum of the mesons
  Lorentz5Momentum pdiff(momenta[0]-momenta[1]);
  Lorentz5Momentum psum (momenta[0]+momenta[1]);
  psum.rescaleMass();
  scale=psum.mass();
  // mass2 of intermediate state
  Energy2 q2 (psum.m2());
  Energy2 dot(psum*pdiff);
  // contribution of the vector resonances
  Complex vnorm(0.),gterm(0.),sterm(0.),snorm(0.);
  complex<InvEnergy2> qterm(ZERO);
  unsigned int imin=0, imax=_vecwgt.size();
  if(resonance) {
    if(abs(resonance->id())%1000==323) {
      switch(abs(resonance->id())/1000) {
      case 0:
	imin=0; break;
      case 100:
	imin=1; break;
      case  30:
	imin=2; break;
      default:
	assert(false);
      }
      imax = imin+1;
    }
    else {
      imax=0;
    }
  }
  for(unsigned int ix=imin;ix<imax;++ix) {
    vnorm += _vecwgt[ix];
    if(ichan<0||_resmap[ix]==ichan) {
      Complex bw=_vecwgt[ix]*pWaveBreitWigner(q2,ix);
      gterm +=bw;
      qterm += _transverse ? bw/sqr(scale) : bw/sqr(_vecmass[ix]);
    }
  }
  // contribution of the scalar resonances
  imin=0, imax=_scawgt.size();
  if(resonance) {
    if(abs(resonance->id())%1000==321) {
      switch(abs(resonance->id())/1000) {
      case 9000:
	imin=0; break;
      case 10:
	imin=1; break;
      default:
	assert(false);
      }
      imax = imin+1;
    }
    else {
      imax=0;
    }
  }
  for(unsigned int ix=imin;ix<imax;++ix) {
    snorm += _scawgt[ix];
    if(ichan<0||_resmap[ix+_vecwgt.size()]==ichan) {
      sterm+=_scawgt[ix]*sWaveBreitWigner(q2,ix);
    }
  }
  // compute the current
  gterm *=_cV/vnorm;
  Complex qtermnew = qterm*_cV*dot/vnorm;
  sterm *= _cS/snorm;
  LorentzPolarizationVectorE output=gterm*pdiff+(-qtermnew+sterm)*psum;
  // return the answer
  if(imode==0) output *= sqrt(0.5);
  return vector<LorentzPolarizationVectorE>(1,output);
}
#line 1 "./OneKaonTwoPionCurrent.cc"
// -*- C++ -*-
//
// OneKaonTwoPionCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneKaonTwoPionCurrent class.
//

#include "OneKaonTwoPionCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace Herwig;

DescribeClass<OneKaonTwoPionCurrent,WeakCurrent>
describeHerwigOneKaonTwoPionCurrent("Herwig::OneKaonTwoPionCurrent",
				    "HwWeakCurrents.so");
HERWIG_INTERPOLATOR_CLASSDESC(OneKaonTwoPionCurrent,Energy,Energy2)

IBPtr OneKaonTwoPionCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr OneKaonTwoPionCurrent::fullclone() const {
  return new_ptr(*this);
}

OneKaonTwoPionCurrent::OneKaonTwoPionCurrent() {
  // the quarks for the different modes
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  setInitialModes(3);
  // rho parameters
  // rho parameters for axial-vector pieces
  _rho1wgts  = {1.0,-0.145,0.};
  _rho1mass  = {0.773*GeV,1.370*GeV,1.750*GeV};
  _rho1width = {0.145*GeV,0.510*GeV,0.120*GeV};
  // K* parameters for the axial-vector pieces
  _kstar1wgts  = {1.0,-0.135,0.};
  _kstar1mass  = {0.892*GeV,1.412*GeV,1.714*GeV};
  _kstar1width = {0.050*GeV,0.227*GeV,0.323*GeV};
  // K* parameters for vector pieces
  _kstar2wgts  = {1.0,-0.25 ,-0.038};
  _kstar2mass  = {0.892*GeV,1.412*GeV,1.714*GeV};
  _kstar2width = {0.050*GeV,0.227*GeV,0.323*GeV};
  // K_1 parameters
  _k1mass  = {1.270*GeV,1.402*GeV};
  _k1width = {0.090*GeV,0.174*GeV};
  _k1wgta = {0.33,1.};
  _k1wgtb = {1.00,0.};
  // the pion decay constant
  _fpi = 130.7*MeV/sqrt(2.);
  _mpi = ZERO;
  _mK  = ZERO;
}


void OneKaonTwoPionCurrent::persistentOutput(PersistentOStream & os) const {
  os << _rho1wgts << ounit(_rho1mass,GeV) << ounit(_rho1width,GeV) 
     << _kstar1wgts << ounit(_kstar1mass,GeV) << ounit(_kstar1width,GeV) 
     << _kstar2wgts << ounit(_kstar2mass,GeV) << ounit(_kstar2width,GeV) 
     << ounit(_k1mass,GeV) << ounit(_k1width,GeV) << _k1wgta << _k1wgtb 
     << ounit(_fpi,GeV) << ounit(_mpi,GeV) << ounit(_mK,GeV);
}

void OneKaonTwoPionCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _rho1wgts >> iunit(_rho1mass,GeV) >> iunit(_rho1width,GeV)
     >> _kstar1wgts >> iunit(_kstar1mass,GeV) >> iunit(_kstar1width,GeV) 
     >> _kstar2wgts >> iunit(_kstar2mass,GeV) >> iunit(_kstar2width,GeV) 
     >> iunit(_k1mass,GeV) >> iunit(_k1width,GeV) >> _k1wgta >> _k1wgtb 
     >> iunit(_fpi,GeV) >> iunit(_mpi,GeV) >> iunit(_mK,GeV);
}


void OneKaonTwoPionCurrent::Init() {

  static ClassDocumentation<OneKaonTwoPionCurrent> documentation
    ("The OneKaonTwoPionCurrent class implements the model of "
     "Z. Phys.  C 69 (1996) 243 [arXiv:hep-ph/9503474]"
     " for the weak current with three "
     "mesons, at least one of which is a kaon",
     "The OneKaonTwoPionCurrent class implements the model of "
     "\\cite{Finkemeier:1995sr} for the weak current with three "
     "mesons, at least one of which is a kaon.",
     "\\bibitem{Finkemeier:1995sr}\n"
     "M.~Finkemeier and E.~Mirkes,\n"
     "Z.\\ Phys.\\  C {\\bf 69} (1996) 243 [arXiv:hep-ph/9503474].\n"
     " %%CITATION = ZEPYA,C69,243;%%\n");

  static Parameter<OneKaonTwoPionCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &OneKaonTwoPionCurrent::_fpi, MeV, 92.4*MeV, ZERO, 200.0*MeV,
     false, false, true);

  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceRhoAxialMasses
    ("RhoAxialMasses",
     "The masses for the rho resonances if used local values",
     &OneKaonTwoPionCurrent::_rho1mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceRhoAxialWidths
    ("RhoAxialWidths",
     "The widths for the rho resonances if used local values",
     &OneKaonTwoPionCurrent::_rho1width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceKstarAxialMasses
    ("KstarAxialMasses",
     "The masses for the Kstar resonances if used local values",
     &OneKaonTwoPionCurrent::_kstar1mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceKstarAxialWidths
    ("KstarAxialWidths",
     "The widths for the Kstar resonances if used local values",
     &OneKaonTwoPionCurrent::_kstar1width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceKstarVectorMasses
    ("KstarVectorMasses",
     "The masses for the Kstar resonances if used local values",
     &OneKaonTwoPionCurrent::_kstar2mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceKstarVectorWidths
    ("KstarVectorWidths",
     "The widths for the Kstar resonances if used local values",
     &OneKaonTwoPionCurrent::_kstar2width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<OneKaonTwoPionCurrent,double> interfaceAxialRhoWeight
    ("AxialRhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &OneKaonTwoPionCurrent::_rho1wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<OneKaonTwoPionCurrent,double> interfaceAxialKStarWeight
    ("AxialKStarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &OneKaonTwoPionCurrent::_kstar1wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<OneKaonTwoPionCurrent,double> interfaceVectorKStarWeight
    ("VectorKStarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &OneKaonTwoPionCurrent::_kstar2wgts,
     0, 0, 0, -1000, 1000, false, false, true);

  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceK1Masses
    ("K1Masses",
     "Masses of the K_1 mesons",
     &OneKaonTwoPionCurrent::_k1mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceK1Widths
    ("K1Widths",
     "Widths of the K_1 mesons",
     &OneKaonTwoPionCurrent::_k1width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<OneKaonTwoPionCurrent,double> interfaceK1WeightKStarPi
    ("K1WeightKStarPi",
     "The relative weights for the K_1 resonances in the K* pi final-state",
     &OneKaonTwoPionCurrent::_k1wgta, -1, 1.0, 0, 10.0,
     false, false, Interface::limited);

  static ParVector<OneKaonTwoPionCurrent,double> interfaceK1WeightRhoK
    ("K1WeightRhoK",
     "The relative weights for the K_1 resonances in the rho K final-state",
     &OneKaonTwoPionCurrent::_k1wgtb, -1, 1.0, 0, 10.0,
     false, false, Interface::limited);
}

// complete the construction of the decay mode for integration
bool OneKaonTwoPionCurrent::createMode(int icharge, tcPDPtr resonance,
				       FlavourInfo flavour,
				       unsigned int imode,PhaseSpaceModePtr mode,
				       unsigned int iloc,int ires,
				       PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(abs(icharge)!=3) return false; 
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IHalf) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Half:
      if(icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusHalf:
      if(icharge == 3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown) {
    if(icharge== 3 and flavour.strange != Strangeness::PlusOne ) return false;
    if(icharge==-3 and flavour.strange != Strangeness::MinusOne) return false;
  }
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
  // get the external particles and check the mass
  int iq(0),ia(0);
  tPDVector extpart(particles(1,imode,iq,ia));
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  // the particles we will use a lot
  tPDPtr a1    = getParticleData(ParticleID::a_1minus);
  tPDPtr k1[2] = {getParticleData(ParticleID::K_1minus),
		  getParticleData(ParticleID::Kstar_1minus)};
  // the rho0 resonances
  tPDPtr rho0[3]  ={getParticleData( 113),getParticleData( 100113),
		    getParticleData( 30113)};
  // the charged rho resonances
  tPDPtr rhoc[3]  ={getParticleData(-213),getParticleData(-100213),
		    getParticleData(-30213)};
  // the K*0 resonances
  tPDPtr Kstar0[3]={getParticleData( 313),getParticleData( 100313),
		    getParticleData( 30313)};
  // the charged K* resonances
  tPDPtr Kstarc[3]={getParticleData(-323),getParticleData(-100323),
		    getParticleData(-30323)};
  if(icharge==3) {
    a1    = a1->CC();
    k1[0] = k1[0]->CC();
    k1[1] = k1[1]->CC();
    for(unsigned int ix=0;ix<3;++ix) {
      if(rhoc[ix]) rhoc[ix]=rhoc[ix]->CC();
      if(Kstar0[ix]) Kstar0[ix]=Kstar0[ix]->CC();
      if(Kstarc[ix]) Kstarc[ix]=Kstarc[ix]->CC();
    }
  }
  if(imode==0) {  
    // channels for pi0 pi0 K-
    for(unsigned int ix=0;ix<3;++ix) {
      for(unsigned int ik=0;ik<2;++ik) {
	if(resonance && resonance != k1[ik]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstarc[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstarc[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=Kstarc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstarc[iy],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstarc[iy],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
    }
  }
  else if(imode==1) {
    // channels for K- pi- pi+
    for(unsigned int ix=0;ix<3;++ix) {
      for(unsigned int ik=0;ik<2;++ik) {
	if(resonance && resonance != k1[ik]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,rho0[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstar0[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=Kstarc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,rho0[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstar0[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
    }
  }
  else if(imode==2) {
    // channels for pi- kbar0 pi0
    for(unsigned int ix=0;ix<3;++ix) {
      for(unsigned int ik=0;ik<2;++ik) {
	if(resonance && resonance != k1[ik]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,rhoc[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstar0[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstarc[ix],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=Kstarc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstar0[iy],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,  rhoc[iy],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstarc[iy],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
    }
  }
  for(unsigned int ix=0;ix<_rho1mass.size();++ix) {
    mode->resetIntermediate(rhoc[ix],_rho1mass[ix],
			    _rho1width[ix]);
    mode->resetIntermediate(rho0[ix],_rho1mass[ix],
			    _rho1width[ix]);
  }
  // K star parameters in the base class
  for(unsigned int ix=0;ix<_kstar1mass.size();++ix) {
    mode->resetIntermediate(Kstarc[ix],_kstar1mass[ix],
			    _kstar1width[ix]);
    mode->resetIntermediate(Kstar0[ix],_kstar1mass[ix],
			    _kstar1width[ix]);
  }
  return true;
}

void OneKaonTwoPionCurrent::dataBaseOutput(ofstream & os,
					   bool header,bool create) const {
  if(header) os << "update decayers set parameters=\"";
  if(create) os << "create Herwig::OneKaonTwoPionCurrent " 
		<< name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<_rho1wgts.size();++ix) {
    if(ix<3) {
      os << "newdef " << name() << ":AxialRhoWeight " << ix 
	 << " " << _rho1wgts[ix] << "\n";
    }
    else {
      os << "insert " << name() << ":AxialRhoWeight " << ix 
	 << " " << _rho1wgts[ix] << "\n";
    }
  }
  for(unsigned int ix=0;ix<_kstar1wgts.size();++ix) {
    if(ix<3) {
      os << "newdef " << name() << ":AxialKStarWeight " << ix 
	 << " " << _kstar1wgts[ix] << "\n";}
    else {
      os << "insert " << name() << ":AxialKStarWeight " << ix 
	 << " " << _kstar1wgts[ix] << "\n";
    }
  }
  for(unsigned int ix=0;ix<_kstar2wgts.size();++ix) {
    if(ix<3) {
      os << "newdef " << name() << ":VectorKStarWeight " << ix 
	 << " " << _kstar2wgts[ix] << "\n";}
    else {
      os << "insert " << name() << ":VectorKStarWeight " << ix 
	 << " " << _kstar2wgts[ix] << "\n";
    }
  }
  os << "newdef " << name() << ":FPi " << _fpi/MeV << "\n";
  for(unsigned int ix=0;ix<_k1mass.size();++ix) {
    if(ix<2) {
      os << "newdef " << name() << ":K1Masses " << ix 
	 << " " << _k1mass[ix]/GeV << "\n";
    }
    else {
      os << "insert " << name() << ":K1Masses " << ix 
	 << " " << _k1mass[ix]/GeV << "\n";
    }
  }
  for(unsigned int ix=0;ix<_k1width.size();++ix) {
    if(ix<2) {
      os << "newdef " << name() << ":K1Widths " << ix 
	 << " " << _k1width[ix]/GeV << "\n";
    }
    else {
      os << "insert " << name() << ":K1Widths " << ix 
	 << " " << _k1width[ix]/GeV << "\n";
    }
  }
  for(unsigned int ix=0;ix<_rho1mass.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":RhoAxialMasses " << ix 
		<< " " << _rho1mass[ix]/GeV << "\n";
    else     os << "insert " << name() << ": RhoAxialMasses" << ix 
		<< " " << _rho1mass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rho1width.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":RhoAxialWidths " << ix 
		    << " " << _rho1width[ix]/GeV << "\n";
    else     os << "insert " << name() << ":RhoAxialWidths " << ix 
		    << " " << _rho1width[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstar1mass.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":KstarAxialMasses " << ix 
		<< " " << _kstar1mass[ix]/GeV << "\n";
    else     os << "insert " << name() << ": KstarAxialMasses" << ix 
		<< " " << _kstar1mass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstar1width.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":KstarAxialWidths " << ix 
		    << " " << _kstar1width[ix]/GeV << "\n";
    else     os << "insert " << name() << ":KstarAxialWidths " << ix 
		    << " " << _kstar1width[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstar2mass.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":KstarVectorMasses " << ix 
		<< " " << _kstar2mass[ix]/GeV << "\n";
    else     os << "insert " << name() << ": KstarVectorMasses" << ix 
		<< " " << _kstar2mass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstar2width.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":KstarVectorWidths " << ix 
		    << " " << _kstar2width[ix]/GeV << "\n";
    else     os << "insert " << name() << ":KstarVectorWidths " << ix 
		    << " " << _kstar2width[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_k1wgta.size();++ix) {
    if(ix<2) os << "newdef " << name() << ":K1WeightKStarPi " << ix
		<< " " << _k1wgta[ix] << "\n";
    else     os << "insert " << name() << ":K1WeightKStarPi " << ix
		<< " " << _k1wgta[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_k1wgtb.size();++ix) {
    if(ix<2) os << "newdef " << name() << ":K1WeightRhoK " << ix
		<< " " << _k1wgtb[ix] << "\n";
    else     os << "insert " << name() << ":K1WeightRhoK " << ix
		<< " " << _k1wgtb[ix] << "\n";
  }
  WeakCurrent::dataBaseOutput(os,false,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" 
		<< fullName() << "\";" << endl;
}  

void OneKaonTwoPionCurrent::doinit() {
  WeakCurrent::doinit();
  // masses for the running widths
  _mpi = getParticleData(ParticleID::piplus)->mass();
  _mK  = getParticleData(ParticleID::K0    )->mass();
}
  
Complex OneKaonTwoPionCurrent::TK1(Energy2 q2,unsigned int iopt,int ires) const {
  double denom(0);
  Complex num(0.);
  if(iopt==0) {
    if(ires>=int(_k1wgta.size())) return 0.;
    denom = std::accumulate(_k1wgta.begin(),_k1wgta.end(),0.0);
    unsigned int imin=0,imax=_k1wgta.size();
    if(ires>0) {
      imin=ires;
      imax=imin+1;
    }
    for(unsigned int ix=imin;ix<imax;++ix)
      num+=_k1wgta[ix]*Resonance::BreitWignerFW_GN(q2,_k1mass[ix],_k1width[ix]);
  }
  else if(iopt==1) {
    if(ires>=int(_k1wgtb.size())) return 0.;
    denom = std::accumulate(_k1wgtb.begin(),_k1wgtb.end(),0.0);
    unsigned int imin=0,imax=_k1wgtb.size();
    if(ires>0) {
      imin=ires;
      imax=imin+1;
    }
    for(unsigned int ix=imin;ix<imax;++ix)
      num+=_k1wgtb[ix]*Resonance::BreitWignerFW_GN(q2,_k1mass[ix],_k1width[ix]);
  }
  else assert(false);
  return num/denom;
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
OneKaonTwoPionCurrent::current(tcPDPtr resonance,
			      FlavourInfo flavour,
			      const int imode, const int ichan, Energy & scale, 
			      const tPDVector & outgoing,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IHalf)
    return vector<LorentzPolarizationVectorE>();
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge()+outgoing[2]->iCharge();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Half:
      if(icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusHalf:
      if(icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown) {
    if(icharge== 3 and flavour.strange != Strangeness::PlusOne ) return vector<LorentzPolarizationVectorE>();
    if(icharge==-3 and flavour.strange != Strangeness::MinusOne) return vector<LorentzPolarizationVectorE>();
  }
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return vector<LorentzPolarizationVectorE>();
  useMe();
  // check the resonance
  int ires1=-1;
  if(resonance) {
    switch(abs(resonance->id())/1000) {
    case 0:
      ires1=0; break;
    case 100:
      ires1=1; break;
    case  30:
      ires1=2; break;
    case  10:
      ires1=3; break;
    case  20:
      ires1=4; break;
    default:
      assert(false);
    }
  }
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q;
  for(unsigned int ix=0;ix<momenta.size();++ix)
    q+=momenta[ix];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 s1 = (momenta[1]+momenta[2]).m2();
  Energy2 s2 = (momenta[0]+momenta[2]).m2();
  Energy2 s3 = (momenta[0]+momenta[1]).m2();
  Complex F1(0.), F2(0.), F5(0.);
  // calculate the pi0 pi0 K-
  if(imode==0) {
    if(ichan<0) {
      Complex K1fact;
      if(ires1<0)
	K1fact = TK1(q2,0,-1);
      else if(ires1<3)
	K1fact = 0.;
      else
	K1fact = TK1(q2,0,ires1-3);
      K1fact /= 6.;
      F1 = K1fact*TKstar1(s1,-1);
      F2 =-K1fact*TKstar1(s2,-1);
      if(ires1<0||ires1>2)
	F5 =-0.25*TKstar2(q2,   -1)*(TKstar1(s1,-1)-TKstar1(s2,-1));
      else
	F5 =-0.25*TKstar2(q2,ires1)*(TKstar1(s1,-1)-TKstar1(s2,-1));
    }
    else if(ichan%10==0) F1=  TK1(q2,0,0)/6.*TKstar1(s1,ichan/10);
    else if(ichan%10==1) F2= -TK1(q2,0,0)/6.*TKstar1(s2,ichan/10);
    else if(ichan%10==2) F1=  TK1(q2,0,1)/6.*TKstar1(s1,ichan/10);
    else if(ichan%10==3) F2= -TK1(q2,0,1)/6.*TKstar1(s2,ichan/10);
    else if(ichan%10<7 ) F5 =-sqrt(2.)/4*TKstar2(q2,ichan/10)*TKstar1(s1,(ichan-4)%10);
    else                 F5 = sqrt(2.)/4*TKstar2(q2,ichan/10)*TKstar1(s2,(ichan-7)%10);
  }
  // calculate the K- pi- pi+
  else if(imode==1) {
    double fact=sqrt(2.)/3.;
    if(ichan<0) {
      Complex K1facta(0.),K1factb(0.);
      if(ires1<0) {
	K1facta = TK1(q2,0,-1);
	K1factb = TK1(q2,1,-1);
      }
      else if(ires1<3) {
	K1facta = 0.;
      }
      else {
	K1facta = TK1(q2,0,ires1-3);
	K1factb = TK1(q2,1,ires1-3);
      }
      F1 = -fact*K1factb*Trho1(s1,-1);
      F2 =  fact*K1facta*TKstar1(s2,-1);
      if(ires1<0||ires1>2)
	F5 = -sqrt(0.5)*TKstar2(q2,-1)*(Trho1(s1,-1)+TKstar1(s2,-1));
      else
	F5 = -sqrt(0.5)*TKstar2(q2,ires1)*(Trho1(s1,-1)+TKstar1(s2,-1));
    }
    else if(ichan%10==0) F1 = -fact*TK1(q2,1,0)*Trho1  (s1,ichan/10);
    else if(ichan%10==1) F2 =  fact*TK1(q2,0,0)*TKstar1(s2,ichan/10);
    else if(ichan%10==2) F1 = -fact*TK1(q2,1,1)*Trho1(  s1,ichan/10);
    else if(ichan%10==3) F2 =  fact*TK1(q2,0,1)*TKstar1(s2,ichan/10);
    else if(ichan%10<7)  F5 = -sqrt(0.5)*TKstar2(q2,ichan/10)*Trho1(  s1,(ichan-4)%10);
    else                 F5 = -sqrt(0.5)*TKstar2(q2,ichan/10)*TKstar1(s2,(ichan-7)%10);
  }
  // calculate the pi- K0bar pi0
  else if(imode==2) {
    if(ichan<0) {
      Complex K1facta(0.),K1factb(0.);
      if(ires1<0) {
	K1facta = TK1(q2,0,-1);
	K1factb = TK1(q2,1,-1);
      }
      else if(ires1<3) {
	K1facta = 0.;
      }
      else {
	K1facta = TK1(q2,0,ires1-3);
	K1factb = TK1(q2,1,ires1-3);
      }
      F1 = K1facta*(TKstar1(s1,-1)-TKstar1(s3,-1))/3.;
      F2 =-(2.*K1factb*Trho1(s2,-1)+K1facta*TKstar1(s3,-1))/3.;
      if(ires1<0||ires1>2)
	F5 = -0.5*TKstar2(q2,-1)*(2.*Trho1(s2,-1)+TKstar1(s1,-1)+TKstar1(s3,-1));
      else
	F5 = -0.5*TKstar2(q2,ires1)*(2.*Trho1(s2,-1)+TKstar1(s1,-1)+TKstar1(s3,-1));
    }
    else if(ichan%15==0) F2 =-2.*TK1(q2,0,0)*Trho1  (s2,ichan/15)/3.;
    else if(ichan%15==1) F1 =    TK1(q2,1,0)*TKstar1(s1,ichan/15)/3.;
    else if(ichan%15==2) {
      F1 =-TK1(q2,1,0)*TKstar1(s3,ichan/15)/3.;
      F2 =-TK1(q2,1,0)*TKstar1(s3,ichan/15)/3.;
    }
    else if(ichan%15==3) F2 =-2.*TK1(q2,0,1)*Trho1  (s2,ichan/15)/3.;
    else if(ichan%15==4) F1 =    TK1(q2,1,1)*TKstar1(s1,ichan/15)/3.;
    else if(ichan%15==5) {
      F1 =-TK1(q2,1,1)*TKstar1(s3,ichan/15)/3.;
      F2 =-TK1(q2,1,1)*TKstar1(s3,ichan/15)/3.;
    }
    else if(ichan%15<9 ) F5 = -0.5*TKstar2(q2,ichan/15)*TKstar1(s1,(ichan- 6)%15);
    else if(ichan%15<12) F5 = -    TKstar2(q2,ichan/15)*Trho1  (s2,(ichan- 9)%15);
    else                 F5 = -0.5*TKstar2(q2,ichan/15)*TKstar1(s3,(ichan-12)%15);
  }
  // the first three form-factors
  LorentzPolarizationVectorE vect = (F2-F1)*momenta[2] + F1*momenta[1] - F2*momenta[0];
  // multiply by the transverse projection operator
  Complex dot=(vect*q)/q2;
  // scalar and parity violating terms
  vect -= dot*q;
  if(F5!=0.) 
    vect -= Complex(0.,1.)*F5/ sqr(Constants::twopi) / sqr(_fpi)*
      Helicity::epsilon(momenta[0],momenta[1],momenta[2]);
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,q.mass()/_fpi*vect);
}

bool OneKaonTwoPionCurrent::accept(vector<int> id) {
  if(id.size()!=3) return false;
  int npip(0),npim(0),nkp(0),nkm(0);
  int npi0(0),nk0(0),nk0bar(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
  }
  if     ( (nkp==1&&npi0==2) || (npi0==2&&nkm==1) )   return true;
  else if( (npip==1&&npim==1&&nkp==1) ||
	   (nkm==1&&npim==1&&npip==1) )               return true;
  else if( (nk0==1&&npip==1&&npi0==1)  ||
	   (npim==1&&nk0bar==1&&npi0==1))             return true;
  else return false;
}

unsigned int OneKaonTwoPionCurrent::decayMode(vector<int> id) {
  assert(id.size()==3);
  int npip(0),npim(0),nkp(0),nkm(0);
  int npi0(0),nk0(0),nk0bar(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
  }
  if     ( (nkp==1&&npi0==2) || (npi0==2&&nkm==1) ) return 0;
  else if( (npip==1&&npim==1&&nkp==1) ||
	   (nkm==1&&npim==1&&npip==1) )             return 1;
  else if( (nk0==1&&npip==1&&npi0==1)  ||
	   (npim==1&&nk0bar==1&&npi0==1))           return 2;
  assert(false);
}


tPDVector OneKaonTwoPionCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart(3);
  if(imode==0) {
    extpart[0]=getParticleData(ParticleID::pi0);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::Kminus);
  }
  else if(imode==1) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::piplus);
  }
  else if(imode==2) {
    extpart[0]=getParticleData(ParticleID::piminus);
    extpart[1]=getParticleData(ParticleID::Kbar0);
    extpart[2]=getParticleData(ParticleID::pi0);
  }
  else
    assert(false);
  // conjugate the particles if needed
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(extpart[ix]->CC()) extpart[ix]=extpart[ix]->CC();
    }
  }
  // return the answer
  return extpart;
}
#line 1 "./TwoKaonOnePionCurrent.cc"
// -*- C++ -*-
//
// TwoKaonOnePionCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoKaonOnePionCurrent class.
//

#include "TwoKaonOnePionCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace Herwig;

DescribeClass<TwoKaonOnePionCurrent,WeakCurrent>
describeHerwigTwoKaonOnePionCurrent("Herwig::TwoKaonOnePionCurrent",
				    "HwWeakCurrents.so");
HERWIG_INTERPOLATOR_CLASSDESC(TwoKaonOnePionCurrent,Energy,Energy2)

IBPtr TwoKaonOnePionCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr TwoKaonOnePionCurrent::fullclone() const {
  return new_ptr(*this);
}

TwoKaonOnePionCurrent::TwoKaonOnePionCurrent() {
  // the quarks for the different modes
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  setInitialModes(7);
  // rho parameters
  // rho parameters for axial-vector pieces
  _rho1wgts  = {1.0,-0.145,0.};
  _rho1mass  = {0.773*GeV,1.370*GeV,1.750*GeV};
  _rho1width = {0.145*GeV,0.510*GeV,0.120*GeV};
  // rho parameters for vector pieces
  _rho2wgts  = {1.0,-0.25,-0.038};
  _rho2mass  = {0.773*GeV,1.500*GeV,1.750*GeV};
  _rho2width = {0.145*GeV,0.220*GeV,0.120*GeV};
  // K* parameters
  // K* parameters for the axial-vector pieces
  _kstar1wgts  = {1.0,-0.135,0.};
  _kstar1mass  = {0.892*GeV,1.412*GeV,1.714*GeV};
  _kstar1width = {0.050*GeV,0.227*GeV,0.323*GeV};
  // a_1 parameters
  _initializea1 = false;
  _a1opt        = true;
  _a1mass  = 1.251*GeV;
  _a1width = 0.475*GeV;
  _a1runwidth = {0*GeV, 0*GeV, 0*GeV, 0*GeV, 0*GeV,
		 0*GeV, 0*GeV, 0*GeV, 0*GeV, 0*GeV,
		 0*GeV, 0*GeV, 1.47729e-06*GeV, 1.19209e-05*GeV, 3.884e-05*GeV,
		 8.83255e-05*GeV, 0.00016561*GeV, 0.000275439*GeV, 0.000422332*GeV,
		 0.000610773*GeV, 0.000845357*GeV, 0.00113092*GeV, 0.00147264*GeV,
		 0.00187616*GeV, 0.0023477*GeV, 0.00289413*GeV, 0.00352315*GeV,
		 0.00424342*GeV, 0.0050647*GeV, 0.00599808*GeV, 0.00705616*GeV,
		 0.00825335*GeV, 0.0096062*GeV, 0.0111337*GeV, 0.0128579*GeV,
		 0.0148041*GeV, 0.017002*GeV, 0.0194858*GeV, 0.0222956*GeV,
		 0.0254781*GeV, 0.0290874*GeV, 0.0331862*GeV, 0.0378467*GeV,
		 0.0431501*GeV, 0.0491862*GeV, 0.0560496*GeV, 0.0638341*GeV,
		 0.0726215*GeV, 0.0824662*GeV, 0.0933765*GeV, 0.105297*GeV,
		 0.118103*GeV, 0.131602*GeV, 0.145564*GeV, 0.159749*GeV,
		 0.173938*GeV, 0.18795*GeV, 0.201649*GeV, 0.214943*GeV,
		 0.227773*GeV, 0.240109*GeV, 0.25194*GeV, 0.263268*GeV,
		 0.274104*GeV, 0.284466*GeV, 0.294372*GeV, 0.303845*GeV,
		 0.312905*GeV, 0.321576*GeV, 0.329878*GeV, 0.337832*GeV,
		 0.345456*GeV, 0.35277*GeV, 0.35979*GeV, 0.366532*GeV,
		 0.373012*GeV, 0.379243*GeV, 0.38524*GeV, 0.391014*GeV,
		 0.396577*GeV, 0.401939*GeV, 0.407111*GeV, 0.412102*GeV,
		 0.416923*GeV, 0.421577*GeV, 0.426078*GeV, 0.430427*GeV,
		 0.434636*GeV, 0.43871*GeV, 0.442654*GeV, 0.446475*GeV,
		 0.450177*GeV, 0.453765*GeV, 0.457245*GeV, 0.460621*GeV,
		 0.463899*GeV, 0.467077*GeV, 0.470164*GeV, 0.473162*GeV,
		 0.476076*GeV, 0.478909*GeV, 0.481658*GeV, 0.484333*GeV,
		 0.486934*GeV, 0.489465*GeV, 0.491926*GeV, 0.494321*GeV,
		 0.496651*GeV, 0.49892*GeV, 0.501128*GeV, 0.503277*GeV,
		 0.505371*GeV, 0.507409*GeV, 0.509395*GeV, 0.511328*GeV,
		 0.513212*GeV, 0.515047*GeV, 0.516846*GeV, 0.518624*GeV,
		 0.520285*GeV, 0.52194*GeV, 0.523553*GeV, 0.525124*GeV,
		 0.526646*GeV, 0.52814*GeV, 0.529638*GeV, 0.531016*GeV,
		 0.532401*GeV, 0.533751*GeV, 0.535069*GeV, 0.536354*GeV,
		 0.537608*GeV, 0.538831*GeV, 0.540039*GeV, 0.541194*GeV,
		 0.542327*GeV, 0.543438*GeV, 0.544522*GeV, 0.545582*GeV,
		 0.546616*GeV, 0.54764*GeV, 0.548615*GeV, 0.549581*GeV,
		 0.550525*GeV, 0.551449*GeV, 0.552351*GeV, 0.55324*GeV,
		 0.554101*GeV, 0.554944*GeV, 0.555772*GeV, 0.556583*GeV,
		 0.557373*GeV, 0.558155*GeV, 0.558917*GeV, 0.559664*GeV,
		 0.560396*GeV, 0.561114*GeV, 0.561849*GeV, 0.562508*GeV,
		 0.563186*GeV, 0.563851*GeV, 0.564503*GeV, 0.565145*GeV,
		 0.565774*GeV, 0.566394*GeV, 0.567001*GeV, 0.567595*GeV,
		 0.568182*GeV, 0.56876*GeV, 0.56933*GeV, 0.569886*GeV,
		 0.570433*GeV, 0.570976*GeV, 0.571504*GeV, 0.572027*GeV,
		 0.572542*GeV, 0.573114*GeV, 0.573548*GeV, 0.574108*GeV,
		 0.574524*GeV, 0.575002*GeV, 0.575473*GeV, 0.575937*GeV,
		 0.576394*GeV, 0.576845*GeV, 0.57729*GeV, 0.57773*GeV,
		 0.578173*GeV, 0.5786*GeV, 0.579013*GeV, 0.579431*GeV,
		 0.579834*GeV, 0.580246*GeV, 0.580649*GeV, 0.581045*GeV,
		 0.581437*GeV, 0.581827*GeV, 0.582208*GeV, 0.582586*GeV, 0.582959*GeV};
  _a1runq2 = {  0*GeV2       , 0.0158678*GeV2, 0.0317356*GeV2, 0.0476034*GeV2, 0.0634712*GeV2,
		0.079339*GeV2, 0.0952068*GeV2,  0.111075*GeV2, 0.126942*GeV2, 0.14281*GeV2,
	        0.158678*GeV2, 0.174546*GeV2, 0.190414*GeV2, 0.206281*GeV2, 0.222149*GeV2,
		0.238017*GeV2, 0.253885*GeV2, 0.269753*GeV2, 0.285621*GeV2, 0.301488*GeV2,
		0.317356*GeV2, 0.333224*GeV2, 0.349092*GeV2, 0.36496*GeV2, 0.380827*GeV2,
		0.396695*GeV2, 0.412563*GeV2, 0.428431*GeV2, 0.444299*GeV2, 0.460166*GeV2,
		0.476034*GeV2, 0.491902*GeV2, 0.50777*GeV2, 0.523638*GeV2, 0.539505*GeV2,
		0.555373*GeV2, 0.571241*GeV2, 0.587109*GeV2, 0.602977*GeV2, 0.618844*GeV2,
		0.634712*GeV2, 0.65058*GeV2, 0.666448*GeV2, 0.682316*GeV2, 0.698183*GeV2,
		0.714051*GeV2, 0.729919*GeV2, 0.745787*GeV2, 0.761655*GeV2, 0.777523*GeV2,
		0.79339*GeV2, 0.809258*GeV2, 0.825126*GeV2, 0.840994*GeV2, 0.856862*GeV2,
		0.872729*GeV2, 0.888597*GeV2, 0.904465*GeV2, 0.920333*GeV2, 0.936201*GeV2,
		0.952068*GeV2, 0.967936*GeV2, 0.983804*GeV2, 0.999672*GeV2, 1.01554*GeV2,
		1.03141*GeV2, 1.04728*GeV2, 1.06314*GeV2, 1.07901*GeV2, 1.09488*GeV2,
		1.11075*GeV2, 1.12661*GeV2, 1.14248*GeV2, 1.15835*GeV2, 1.17422*GeV2,
		1.19009*GeV2, 1.20595*GeV2, 1.22182*GeV2, 1.23769*GeV2, 1.25356*GeV2,
		1.26942*GeV2, 1.28529*GeV2, 1.30116*GeV2, 1.31703*GeV2, 1.3329*GeV2,
		1.34876*GeV2, 1.36463*GeV2, 1.3805*GeV2, 1.39637*GeV2, 1.41223*GeV2,
		1.4281*GeV2, 1.44397*GeV2, 1.45984*GeV2, 1.47571*GeV2, 1.49157*GeV2,
		1.50744*GeV2, 1.52331*GeV2, 1.53918*GeV2, 1.55505*GeV2, 1.57091*GeV2,
		1.58678*GeV2, 1.60265*GeV2, 1.61852*GeV2, 1.63438*GeV2, 1.65025*GeV2,
		1.66612*GeV2, 1.68199*GeV2, 1.69786*GeV2, 1.71372*GeV2, 1.72959*GeV2,
		1.74546*GeV2, 1.76133*GeV2, 1.77719*GeV2, 1.79306*GeV2, 1.80893*GeV2,
		1.8248*GeV2, 1.84067*GeV2, 1.85653*GeV2, 1.8724*GeV2, 1.88827*GeV2,
		1.90414*GeV2, 1.92*GeV2, 1.93587*GeV2, 1.95174*GeV2, 1.96761*GeV2,
		1.98348*GeV2, 1.99934*GeV2, 2.01521*GeV2, 2.03108*GeV2, 2.04695*GeV2,
		2.06281*GeV2, 2.07868*GeV2, 2.09455*GeV2, 2.11042*GeV2, 2.12629*GeV2,
		2.14215*GeV2, 2.15802*GeV2, 2.17389*GeV2, 2.18976*GeV2, 2.20563*GeV2,
		2.22149*GeV2, 2.23736*GeV2, 2.25323*GeV2, 2.2691*GeV2, 2.28496*GeV2,
		2.30083*GeV2, 2.3167*GeV2, 2.33257*GeV2, 2.34844*GeV2, 2.3643*GeV2,
		2.38017*GeV2, 2.39604*GeV2, 2.41191*GeV2, 2.42777*GeV2, 2.44364*GeV2,
		2.45951*GeV2, 2.47538*GeV2, 2.49125*GeV2, 2.50711*GeV2, 2.52298*GeV2,
		2.53885*GeV2, 2.55472*GeV2, 2.57058*GeV2, 2.58645*GeV2, 2.60232*GeV2,
		2.61819*GeV2, 2.63406*GeV2, 2.64992*GeV2, 2.66579*GeV2, 2.68166*GeV2,
		2.69753*GeV2, 2.71339*GeV2, 2.72926*GeV2, 2.74513*GeV2, 2.761*GeV2,
		2.77687*GeV2, 2.79273*GeV2, 2.8086*GeV2, 2.82447*GeV2, 2.84034*GeV2,
		2.85621*GeV2, 2.87207*GeV2, 2.88794*GeV2, 2.90381*GeV2, 2.91968*GeV2,
		2.93554*GeV2, 2.95141*GeV2, 2.96728*GeV2, 2.98315*GeV2, 2.99902*GeV2,
		3.01488*GeV2, 3.03075*GeV2, 3.04662*GeV2, 3.06249*GeV2, 3.07835*GeV2,
		3.09422*GeV2, 3.11009*GeV2, 3.12596*GeV2, 3.14183*GeV2, 3.15769*GeV2};
  // parameters for the T_omega function
  _epsomega   = 0.05;
  _omegamass  = 0.782*GeV;
  _omegawidth = 0.00843*GeV;
  _phimass    = 1.020*GeV;
  _phiwidth   = 0.00443*GeV;
  _omegaKstarwgt=1./sqrt(2.);
  // the pion decay constant
  _fpi     = 130.7*MeV/sqrt(2.);
  _mpi     = ZERO;
  _mK      = ZERO;
  _maxmass = ZERO;
  _maxcalc = ZERO;
}


void TwoKaonOnePionCurrent::persistentOutput(PersistentOStream & os) const {
  os << _a1runinter
     << _rho1wgts << ounit(_rho1mass,GeV) << ounit(_rho1width,GeV) 
     << _rho2wgts << ounit(_rho2mass,GeV) << ounit(_rho2width,GeV)
     << _kstar1wgts << ounit(_kstar1mass,GeV) << ounit(_kstar1width,GeV)
     << ounit(_a1mass,GeV) << ounit(_a1width,GeV)
     << ounit(_a1runwidth,GeV) << ounit(_a1runq2,GeV2) << _epsomega 
     << ounit(_omegamass,GeV) << ounit(_omegawidth,GeV) 
     << ounit(_phimass,GeV) << ounit(_phiwidth,GeV) << _omegaKstarwgt 
     << ounit(_fpi,GeV) << ounit(_mpi,GeV) << ounit(_mK,GeV) 
     << _initializea1 << _a1opt
     << ounit(_maxmass,GeV) << ounit(_maxcalc,GeV);
}

void TwoKaonOnePionCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _a1runinter
     >> _rho1wgts >> iunit(_rho1mass,GeV) >> iunit(_rho1width,GeV) 
     >> _rho2wgts >> iunit(_rho2mass,GeV) >> iunit(_rho2width,GeV) 
     >> _kstar1wgts >> iunit(_kstar1mass,GeV) >> iunit(_kstar1width,GeV)
     >> iunit(_a1mass,GeV) >> iunit(_a1width,GeV)
     >> iunit(_a1runwidth,GeV) >> iunit(_a1runq2,GeV2) >> _epsomega 
     >> iunit(_omegamass,GeV) >> iunit(_omegawidth,GeV) 
     >> iunit(_phimass,GeV) >> iunit(_phiwidth,GeV) >> _omegaKstarwgt 
     >> iunit(_fpi,GeV) >> iunit(_mpi,GeV) >> iunit(_mK,GeV) 
     >> _initializea1 >> _a1opt
     >> iunit(_maxmass,GeV) >> iunit(_maxcalc,GeV);
}


void TwoKaonOnePionCurrent::Init() {

  static ClassDocumentation<TwoKaonOnePionCurrent> documentation
    ("The TwoKaonOnePionCurrent class implements the model of "
     "Z. Phys.  C 69 (1996) 243 [arXiv:hep-ph/9503474]"
     " for the weak current with three "
     "mesons, at least one of which is a kaon",
     "The TwoKaonOnePionCurrent class implements the model of "
     "\\cite{Finkemeier:1995sr} for the weak current with three "
     "mesons, at least one of which is a kaon.",
     "\\bibitem{Finkemeier:1995sr}\n"
     "M.~Finkemeier and E.~Mirkes,\n"
     "Z.\\ Phys.\\  C {\\bf 69} (1996) 243 [arXiv:hep-ph/9503474].\n"
     " %%CITATION = ZEPYA,C69,243;%%\n"

);

  static Switch<TwoKaonOnePionCurrent,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &TwoKaonOnePionCurrent::_initializea1, false, false, false);
  static SwitchOption interfaceInitializea1Initialization
    (interfaceInitializea1,
     "Yes",
     "Initialize the calculation",
     true);
  static SwitchOption interfaceInitializea1NoInitialization
    (interfaceInitializea1,
     "No",
     "Use the default values",
     false);

  static Parameter<TwoKaonOnePionCurrent,Energy> interfaceA1Width
    ("A1Width",
     "The a_1 width if using local values.",
     &TwoKaonOnePionCurrent::_a1width, GeV, 0.599*GeV, ZERO, 10.0*GeV,
     false, false, false);
  
  static Parameter<TwoKaonOnePionCurrent,Energy> interfaceA1Mass
    ("A1Mass",
     "The a_1 mass if using local values.",
     &TwoKaonOnePionCurrent::_a1mass, GeV, 1.251*GeV, ZERO, 10.0*GeV,
     false, false, false);

  static Parameter<TwoKaonOnePionCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &TwoKaonOnePionCurrent::_fpi, MeV, 92.4*MeV, ZERO, 200.0*MeV,
     false, false, true);

  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceRhoAxialMasses
    ("RhoAxialMasses",
     "The masses for the rho resonances if used local values",
     &TwoKaonOnePionCurrent::_rho1mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceRhoAxialWidths
    ("RhoAxialWidths",
     "The widths for the rho resonances if used local values",
     &TwoKaonOnePionCurrent::_rho1width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceRhoVectorMasses
    ("RhoVectorMasses",
     "The masses for the rho resonances if used local values",
     &TwoKaonOnePionCurrent::_rho2mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceRhoVectorWidths
    ("RhoVectorWidths",
     "The widths for the rho resonances if used local values",
     &TwoKaonOnePionCurrent::_rho2width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceKstarAxialMasses
    ("KstarAxialMasses",
     "The masses for the Kstar resonances if used local values",
     &TwoKaonOnePionCurrent::_kstar1mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceKstarAxialWidths
    ("KstarAxialWidths",
     "The widths for the Kstar resonances if used local values",
     &TwoKaonOnePionCurrent::_kstar1width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<TwoKaonOnePionCurrent,double> interfaceAxialRhoWeight
    ("AxialRhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &TwoKaonOnePionCurrent::_rho1wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<TwoKaonOnePionCurrent,double> interfaceAxialKStarWeight
    ("AxialKStarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &TwoKaonOnePionCurrent::_kstar1wgts,
     0, 0, 0, -1000, 1000, false, false, true);

  static ParVector<TwoKaonOnePionCurrent,double> interfaceVectorRhoWeight
    ("VectorRhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &TwoKaonOnePionCurrent::_rho2wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static Switch<TwoKaonOnePionCurrent,bool> interfacea1WidthOption
    ("a1WidthOption",
     "Option for the treatment of the a1 width",
     &TwoKaonOnePionCurrent::_a1opt, true, false, false);
  static SwitchOption interfacea1WidthOptionLocal
    (interfacea1WidthOption,
     "Local",
     "Use a calculation of the running width based on the parameters as"
     " interpolation table.",
     true);
  static SwitchOption interfacea1WidthOptionParam
    (interfacea1WidthOption,
     "Kuhn",
     "Use the parameterization of Kuhn and Santamaria for default parameters."
     " This should only be used for testing vs TAUOLA",
     false);

  static ParVector<TwoKaonOnePionCurrent,Energy> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &TwoKaonOnePionCurrent::_a1runwidth, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<TwoKaonOnePionCurrent,Energy2> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &TwoKaonOnePionCurrent::_a1runq2, GeV2, -1, 1.0*GeV2, ZERO, 10.0*GeV2,
     false, false, true);

  static Parameter<TwoKaonOnePionCurrent,double> interfaceEpsOmega
    ("EpsOmega",
     "The omega-phi mixing ",
     &TwoKaonOnePionCurrent::_epsomega, 0.05, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<TwoKaonOnePionCurrent,Energy> interfaceOmegaMass
    ("OmegaMass",
     "The mass of the omega meson",
     &TwoKaonOnePionCurrent::_omegamass, GeV, 0.782*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<TwoKaonOnePionCurrent,Energy> interfaceOmegaWidth
    ("OmegaWidth",
     "The width of the omega meson",
     &TwoKaonOnePionCurrent::_omegawidth, GeV, 0.00843*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<TwoKaonOnePionCurrent,Energy> interfacePhiMass
    ("PhiMass",
     "The mass of the phi meson",
     &TwoKaonOnePionCurrent::_phimass, GeV, 1.020*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<TwoKaonOnePionCurrent,Energy> interfacePhiWidth
    ("PhiWidth",
     "The width of the phi meson",
     &TwoKaonOnePionCurrent::_phiwidth, GeV, 0.00443*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<TwoKaonOnePionCurrent,double> interfaceOmegaKStarWeight
    ("OmegaKStarWeight",
     "The relative weight of the omega-phi and K* terms",
     &TwoKaonOnePionCurrent::_omegaKstarwgt, 1./sqrt(2.), 0.0, 100.0,
     false, false, Interface::limited);

}

void TwoKaonOnePionCurrent::inita1Width(int iopt) {
  if(iopt==-1) {
    _maxcalc=_maxmass;
    if(!_initializea1||_maxmass==ZERO) return; 
    // parameters for the table of values
    Energy2 step(sqr(_maxmass)/199.);
    // integrator to perform the integral
    vector<double> inweights;inweights.push_back(0.5);inweights.push_back(0.5);
    vector<int> intype;intype.push_back(2);intype.push_back(3);
    Energy mrho(getParticleData(ParticleID::rhoplus)->mass()),
      wrho(getParticleData(ParticleID::rhoplus)->width());
    vector<Energy> inmass(2,mrho),inwidth(2,wrho);
    vector<double> inpow(2,0.0);
    ThreeBodyAllOnCalculator<TwoKaonOnePionCurrent> 
      widthgen(inweights,intype,inmass,inwidth,inpow,*this,0,_mpi,_mpi,_mpi);
    // normalisation constant to give physical width if on shell
    double a1const(_a1width/(widthgen.partialWidth(sqr(_a1mass))));
    // loop to give the values
    _a1runq2.clear();_a1runwidth.clear();
    for(Energy2 moff2 = ZERO; moff2<=sqr(_maxmass); moff2+=step) {
      _a1runwidth.push_back(widthgen.partialWidth(moff2)*a1const);
      _a1runq2.push_back(moff2);
    }
  }
  // set up the interpolator
  else if(iopt==0) {
    _a1runinter = make_InterpolatorPtr(_a1runwidth,_a1runq2,3);
  }
}

// complete the construction of the decay mode for integration
bool TwoKaonOnePionCurrent::createMode(int icharge, tcPDPtr resonance,
				       FlavourInfo flavour,
				       unsigned int imode,PhaseSpaceModePtr mode,
				       unsigned int iloc,int ires,
				       PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(abs(icharge)!=3) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode<=1) return false;
      break;
    case IsoSpin::I3One:
      if( imode>1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if( imode>1 || icharge == 3) return false;
      break;
    default:
      return false;
    }
  }
  // other flavour stuff
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  != Beauty::Zero     ) return false;
  // get the particles and check the mass
  int iq(0),ia(0);
  tPDVector extpart(particles(1,imode,iq,ia));
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  // the particles we will use a lot
  tPDPtr a1    = getParticleData(ParticleID::a_1minus);
  _maxmass=max(_maxmass,upp);
  // the rho0 resonances
  tPDPtr rho0[3]  ={getParticleData( 113),getParticleData( 100113),
		    getParticleData( 30113)};
  // the charged rho resonances
  tPDPtr rhoc[3]  ={getParticleData(-213),getParticleData(-100213),
		    getParticleData(-30213)};
  // the K*0 resonances
  tPDPtr Kstar0[3]={getParticleData( 313),getParticleData( 100313),
		    getParticleData( 30313)};
  // the charged K* resonances
  tPDPtr Kstarc[3]={getParticleData(-323),getParticleData(-100323),
		    getParticleData(-30323)};
  if(icharge==3) {
    a1    = a1->CC();
    for(unsigned int ix=0;ix<3;++ix) {
      if(rhoc[ix]) rhoc[ix]=rhoc[ix]->CC();
      if(Kstar0[ix]) Kstar0[ix]=Kstar0[ix]->CC();
      if(Kstarc[ix]) Kstarc[ix]=Kstarc[ix]->CC();
    }
  }
  if(imode==0) {
    // channels for K- pi- K+
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==a1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstar0[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,rho0[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=rhoc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstar0[iy],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
      }
    }
  }
  else if(imode==1) {
    // channels for K0 pi- K0bar
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==a1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,rho0[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=rhoc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[iy],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
      }
    }
  }
  else if(imode==2) {
    // channels for K- pi0 K0
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==a1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstar0[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,rhoc[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=rhoc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstar0[iy],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[iy],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
    }
  }
  else if(imode==3||imode==4) {
    // channels for K_S0 pi- K_S0 and K_L0 pi- K_L0 
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==a1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=rhoc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[iy],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[iy],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
    }
  }
  else if(imode==5) {
    // channels for K_S0 pi- K_L0
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==a1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,rho0[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=rhoc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[ix],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
    }
  }
  // update the integration parameters
  for(unsigned int ix=0;ix<_rho1mass.size();++ix) {
    mode->resetIntermediate(rhoc[ix],_rho1mass[ix],
			    _rho1width[ix]);
    mode->resetIntermediate(rho0[ix],_rho1mass[ix],
			    _rho1width[ix]);
  }
  for(unsigned int ix=0;ix<_kstar1mass.size();++ix) {
    mode->resetIntermediate(Kstarc[ix],_kstar1mass[ix],
			    _kstar1width[ix]);
    mode->resetIntermediate(Kstar0[ix],_kstar1mass[ix],
			    _kstar1width[ix]);
  }
  return true;
}

void TwoKaonOnePionCurrent::dataBaseOutput(ofstream & os,
					   bool header,bool create) const {
  if(header) os << "update decayers set parameters=\"";
  if(create) os << "create Herwig::TwoKaonOnePionCurrent " 
		<< name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<_rho1wgts.size();++ix) {
    if(ix<3) {
      os << "newdef " << name() << ":AxialRhoWeight " << ix 
	 << " " << _rho1wgts[ix] << "\n";
    }
    else {
      os << "insert " << name() << ":AxialRhoWeight " << ix 
	 << " " << _rho1wgts[ix] << "\n";
    }
  }
  for(unsigned int ix=0;ix<_kstar1wgts.size();++ix) {
    if(ix<3) {
      os << "newdef " << name() << ":AxialKStarWeight " << ix 
	 << " " << _kstar1wgts[ix] << "\n";}
    else {
      os << "insert " << name() << ":AxialKStarWeight " << ix 
	 << " " << _kstar1wgts[ix] << "\n";
    }
  }
  for(unsigned int ix=0;ix<_rho2wgts.size();++ix) {
    if(ix<3) {
      os << "newdef " << name() << ":VectorRhoWeight " << ix 
	 << " " << _rho2wgts[ix] << "\n";
    }
    else {
      os << "insert " << name() << ":VectorRhoWeight " << ix 
	 << " " << _rho2wgts[ix] << "\n";
    }
  }
  os << "newdef " << name() << ":OmegaKStarWeight " << _omegaKstarwgt << "\n";
  os << "newdef " << name() << ":EpsOmega " << _epsomega << "\n";
  os << "newdef " << name() << ":Initializea1 " << _initializea1 << "\n";
  os << "newdef " << name() << ":a1WidthOption " << _a1opt << "\n";
  for(unsigned int ix=0;ix<_a1runwidth.size();++ix) {
    os << "newdef " << name() << ":a1RunningWidth " << ix 
	   << " " << _a1runwidth[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_a1runq2.size();++ix) {
    os << "newdef " << name() << ":a1RunningQ2 " << ix 
	   << " " << _a1runq2[ix]/GeV2 << "\n";
  }
  os << "newdef " << name() << ":A1Width " << _a1width/GeV << "\n";
  os << "newdef " << name() << ":A1Mass " << _a1mass/GeV << "\n";
  os << "newdef " << name() << ":OmegaWidth " << _omegawidth/GeV << "\n";
  os << "newdef " << name() << ":OmegaMass " << _omegamass/GeV << "\n";
  os << "newdef " << name() << ":PhiWidth " << _phiwidth/GeV << "\n";
  os << "newdef " << name() << ":PhiMass " << _phimass/GeV << "\n";
  os << "newdef " << name() << ":FPi " << _fpi/MeV << "\n";
  for(unsigned int ix=0;ix<_rho1mass.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":RhoAxialMasses " << ix 
		<< " " << _rho1mass[ix]/GeV << "\n";
    else     os << "insert " << name() << ": RhoAxialMasses" << ix 
		<< " " << _rho1mass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rho1width.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":RhoAxialWidths " << ix 
		    << " " << _rho1width[ix]/GeV << "\n";
    else     os << "insert " << name() << ":RhoAxialWidths " << ix 
		    << " " << _rho1width[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rho2mass.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":RhoVectorMasses " << ix 
		<< " " << _rho2mass[ix]/GeV << "\n";
    else     os << "insert " << name() << ": RhoVectorMasses" << ix 
		<< " " << _rho2mass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rho2width.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":RhoVectorWidths " << ix 
		    << " " << _rho2width[ix]/GeV << "\n";
    else     os << "insert " << name() << ":RhoVectorWidths " << ix 
		    << " " << _rho2width[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstar1mass.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":KstarAxialMasses " << ix 
		<< " " << _kstar1mass[ix]/GeV << "\n";
    else     os << "insert " << name() << ": KstarAxialMasses" << ix 
		<< " " << _kstar1mass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstar1width.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":KstarAxialWidths " << ix 
		    << " " << _kstar1width[ix]/GeV << "\n";
    else     os << "insert " << name() << ":KstarAxialWidths " << ix 
		    << " " << _kstar1width[ix]/GeV << "\n";
  }
  WeakCurrent::dataBaseOutput(os,false,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" 
		<< fullName() << "\";" << endl;
}  

void TwoKaonOnePionCurrent::doinit() {
  WeakCurrent::doinit();
  // masses for the running widths
  _mpi = getParticleData(ParticleID::piplus)->mass();
  _mK  = getParticleData(ParticleID::K0)    ->mass();
  // initialise the a_1 running width calculation
  inita1Width(-1);
  inita1Width(0);
}

void TwoKaonOnePionCurrent::doinitrun() {
  // set up the running a_1 width
  inita1Width(0);
  WeakCurrent::doinitrun();
}

void TwoKaonOnePionCurrent::doupdate() {
  WeakCurrent::doupdate();
  // update running width if needed
  if ( !touched() ) return;
  if(_maxmass!=_maxcalc) inita1Width(-1);
}

double TwoKaonOnePionCurrent::
threeBodyMatrixElement(const int       , const Energy2 q2,
		       const Energy2 s3, const Energy2 s2, 
		       const Energy2 s1, const Energy    , 
		       const Energy    , const Energy    ) const {
  Energy2 mpi2(sqr(_mpi));
  Complex propb(Trho1(s1,-1)),propa(Trho1(s2,-1)); 
  // the matrix element
  Energy2 output(ZERO); 
  // first resonance
  output+= ((s1-4.*mpi2)+0.25*(s3-s2)*(s3-s2)/q2)*real(propb*conj(propb)); 
  // second resonance
  output+= ((s2-4.*mpi2)+0.25*(s3-s1)*(s3-s1)/q2)*real(propa*conj(propa)); 
  // the interference term 
  output+= (0.5*q2-s3-0.5*mpi2+0.25*(s3-s2)*(s3-s1)/q2)*real(propa*conj(propb)+
							     propb*conj(propa)); 
  return output / sqr(_rho1mass[0]);
}
  
Complex TwoKaonOnePionCurrent::Tomega(Energy2 q2, int ires) const {
  double denom=(1.+_epsomega);
  Complex num(0.);
  if(ires<0) num=OmegaPhiBreitWigner(q2,0)+_epsomega*OmegaPhiBreitWigner(q2,1);
  else if(ires==0) num=OmegaPhiBreitWigner(q2,0);
  else             num=OmegaPhiBreitWigner(q2,1);
  return num/denom;
}

Complex TwoKaonOnePionCurrent::TOmegaKStar(Energy2 s1,Energy2 s2,int ires) const {
  Complex output;
  if(ires<0)         output = _omegaKstarwgt*TKstar1(s1,-1)+Tomega(s2,-1);
  else if(ires%2==0) output = _omegaKstarwgt*TKstar1(s1,ires/2);
  else if(ires%2==1) output = Tomega(s2,ires/2);
  return output/(1.+_omegaKstarwgt);
}


// the hadronic currents    
vector<LorentzPolarizationVectorE> 
TwoKaonOnePionCurrent::current(tcPDPtr resonance,
			      FlavourInfo flavour,
			      const int imode, const int ichan, Energy & scale, 
			      const tPDVector & ,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IOne)
    return vector<LorentzPolarizationVectorE>();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One: case IsoSpin::I3MinusOne:
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  // other flavour stuff
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  != Beauty::Zero     ) return vector<LorentzPolarizationVectorE>();
  // check the resonance
  int ires1=-1;
  if(resonance) {
    switch(abs(resonance->id())/1000) {
    case 0:
      ires1=0; break;
    case 100:
      ires1=1; break;
    case  30:
      ires1=2; break;
    case  10:
      ires1=3; break;
    default:
      assert(false);
    }
  }
  useMe();
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q;
  for(unsigned int ix=0;ix<momenta.size();++ix)
    q+=momenta[ix];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 s1 = (momenta[1]+momenta[2]).m2();
  Energy2 s2 = (momenta[0]+momenta[2]).m2();
  Energy2 s3 = (momenta[0]+momenta[1]).m2();
  // calculate the form factors
  useMe();
  Complex F1(0.), F2(0.), F5(0.);
  Complex a1fact = ires1<0 || ires1==3 ? a1BreitWigner(q2) : 0.;
  // calculate the K- pi - K+ factor
  if(imode==0) {
    a1fact *= sqrt(2.)/3.;
    if(ichan<0) {
      F1 = -a1fact*TKstar1(s1,-1);
      F2 =  a1fact*Trho1(s2,-1);
      if(ires1<0)
	F5 = Trho2(q2,   -1)*TOmegaKStar(s1,s2,-1)*sqrt(2.);
      else if(ires1<3)
	F5 = Trho2(q2,ires1)*TOmegaKStar(s1,s2,-1)*sqrt(2.);
      else
	F5 = 0.;
    }
    else if(ichan%5==0) F1 = -a1fact*TKstar1(s1,    ichan/5);
    else if(ichan%5==1) F2 =  a1fact*Trho1(  s2,(ichan-1)/5);
    else if(ichan%5>=2) F5 = Trho2(q2,ichan/5)*TOmegaKStar(s1,s2,2*((ichan-2)%5))
      *sqrt(2.);
  }
  // calculate the K0 pi- K0bar
  else if(imode==1) {
    a1fact *= sqrt(2.)/3.;
    if(ichan<0) {
      F1 =-a1fact*TKstar1(s1,-1);
      F2 = a1fact*Trho1  (s2,-1);
      if(ires1<0)
	F5 =-Trho2(q2,   -1)*TOmegaKStar(s1,s2,-1)*sqrt(2.);
      else if(ires1<3)
	F5 =-Trho2(q2,ires1)*TOmegaKStar(s1,s2,-1)*sqrt(2.);
      else
	F5 = 0.;
    }
    else if(ichan%5==0) F1 = -a1fact*TKstar1(s1,    ichan/5);
    else if(ichan%5==1) F2 =  a1fact*Trho1  (s2,(ichan-1)/5);
    else if(ichan%5>=2) F5 = -Trho2(q2,ichan/5)*TOmegaKStar(s1,s2,2*((ichan-2)%5))
      *sqrt(2.);
  }
  // calculate the K- pi0 k0
  else if(imode==2) {
    a1fact /= 3.;
    if(ichan<0) {
      F1 =  a1fact*( TKstar1(s1,-1)-TKstar1(s3,-1));
      F2 = -a1fact*(2.*Trho1(s2,-1)+TKstar1(s3,-1));
      if(ires1<0)
	F5 = Trho2(q2,   -1)*(TKstar1(s3,-1)-TKstar1(s1,-1))/(1.+_omegaKstarwgt)/sqrt(2.);
      else if(ires1<3)
	F5 = Trho2(q2,ires1)*(TKstar1(s3,-1)-TKstar1(s1,-1))/(1.+_omegaKstarwgt)/sqrt(2.);
      else
	F5 = 0.;
    }
    else if(ichan%9==0) F1 =  a1fact*TKstar1(s1,ichan/9)/3.;
    else if(ichan%9==1) {
      F1 = +a1fact*TKstar1(s3,(ichan-1)/9)/3.;
      F2 = -a1fact*TKstar1(s3,(ichan-1)/9)/3.;
    }
    else if(ichan%9==2) F2 = -a1fact*2.*Trho1(s2,(ichan-2)/9)/3.;
    else if(ichan%9<6)  F5 =-Trho2(q2,ichan/9)*TKstar1(s1,(ichan-3)%9)
      /(1.+_omegaKstarwgt)/sqrt(2.);
    else                F5 = Trho2(q2,ichan/9)*TKstar1(s3,(ichan-6)%9)
      /(1.+_omegaKstarwgt)/sqrt(2.);
  }
  // calculate the K_S0 pi- K_S0 or K_L0 pi- K_L0
  else if(imode==3||imode==4) {
    a1fact /=6;
    if(ichan<0) {
      F1 = a1fact*(TKstar1(s1,-1)+TKstar1(s3,-1));
      F2 = a1fact*TKstar1(s3,-1);
      if(ires1<0)
	F5 = 0.5*Trho2(q2,   -1)*(TOmegaKStar(s1,s2,-1)-TOmegaKStar(s3,s2,-1));
      else if(ires1<3)
	F5 = 0.5*Trho2(q2,ires1)*(TOmegaKStar(s1,s2,-1)-TOmegaKStar(s3,s2,-1));
      else
	F5 = 0.;
    }
    else if(ichan%8==0) F1=a1fact*TKstar1(s1,ichan/8);
    else if(ichan%8==1) {
      F1 = a1fact*TKstar1(s3,ichan/8);
      F2 = a1fact*TKstar1(s3,ichan/8);
    }
    else if(ichan%8<5 ) F5 = -Trho2(q2,ichan/8)*TKstar1(s1,(ichan-2)%8)
      /(1.+_omegaKstarwgt)/2.;
    else                F5 =  Trho2(q2,ichan/8)*TKstar1(s3,(ichan-5)%8)
      /(1.+_omegaKstarwgt)/2.;
  }
  else if(imode==5) {
    a1fact *= 1./3./sqrt(2.);
    if(ichan<0) {
      F1 = -a1fact*(TKstar1(s1,-1)-TKstar1(s3,-1));
      F2 =  a1fact*(2.*Trho1(s2,-1)+TKstar1(s3,-1));
      if(ires1<0)
	F5 = -Trho2(q2,   -1)*(TOmegaKStar(s1,s2,-1)+TOmegaKStar(s3,s2,-1))/sqrt(2.);
      else if(ires1<3)
	F5 = -Trho2(q2,ires1)*(TOmegaKStar(s1,s2,-1)+TOmegaKStar(s3,s2,-1))/sqrt(2.);
      else
	F5 = 0.;
    }
    else if(ichan%9==0) F1 =-   a1fact*TKstar1(s1,ichan/9);
    else if(ichan%9==1) {
      F1 = a1fact*TKstar1(s3,ichan/9);
      F2 = a1fact*TKstar1(s3,ichan/9);
    }
    else if(ichan%9==2) F2 = 2.*a1fact*Trho1(  s2,ichan/9);
    else if(ichan%9<6 ) F5 = -sqrt(0.5)*Trho2(q2,ichan/9)*
      TOmegaKStar(s1,s2,2*((ichan-3)%9))/sqrt(2.);
    else                F5 = -sqrt(0.5)*Trho2(q2,ichan/9)*
      TOmegaKStar(s3,s2,2*((ichan-6)%9))/sqrt(2.);
  }
  // the first three form-factors
  LorentzPolarizationVectorE vect = (F2-F1)*momenta[2] + F1*momenta[1] - F2*momenta[0];
  // multiply by the transverse projection operator
  Complex dot=(vect*q)/q2;
  // scalar and parity violating terms
  vect -= dot*q;
  if(F5!=0.) 
    vect -= Complex(0.,1.)*F5/sqr(Constants::twopi)/sqr(_fpi)*
      Helicity::epsilon(momenta[0],momenta[1],momenta[2]);
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,q.mass()/_fpi*vect);
}

bool TwoKaonOnePionCurrent::accept(vector<int> id) {
  if(id.size()!=3) return false;
  int npip(0),npim(0),nkp(0),nkm(0);
  int npi0(0),nk0(0),nk0bar(0),nks(0),nkl(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
    else if(id[ix]==ParticleID::K_S0)    ++nks;
    else if(id[ix]==ParticleID::K_L0)    ++nkl;
  }
  if     ( (nkp==1&&nkm==1&&npip==1) ||
	   (nkp==1&&nkm==1&&npim==1))             return true;
  else if( (nk0==1&&nk0bar==1&&npip==1) ||
	   (nk0==1&&nk0bar==1&&npim==1))          return true;
  else if( (nkp==1&&nk0bar==1&&npi0==1) ||
	   (nkm==1&&npi0==1&&nk0==1))             return true;
  else if( nks==2 && (npip==1||npim==1) )         return true;
  else if( nkl==2 && (npip==1||npim==1) )         return true;
  else if( nks==1&&nkl==1 && (npip==1||npim==1) ) return true;
  return false;
}

unsigned int TwoKaonOnePionCurrent::decayMode(vector<int> id) {
  assert(id.size()==3);
  int npip(0),npim(0),nkp(0),nkm(0),
    npi0(0),nk0(0),nk0bar(0),neta(0),nks(0),nkl(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
    else if(id[ix]==ParticleID::eta)     ++neta;
    else if(id[ix]==ParticleID::K_S0)    ++nks;
    else if(id[ix]==ParticleID::K_L0)    ++nkl;
  }
  if     ( (nkp==1&&nkm==1&&npip==1) ||
	   (nkp==1&&nkm==1&&npim==1))                 return 0;
  else if( (nk0==1&&nk0bar==1&&npip==1) ||
	   (nk0==1&&nk0bar==1&&npim==1))              return 1;
  else if( (nkp==1&&nk0bar==1&&npi0==1) ||
	   (nkm==1&&npi0==1&&nk0==1))                 return 2;
  else if( nks==2 && (npip==1||npim==1) )             return 3;
  else if( nkl==2 && (npip==1||npim==1) )             return 4;
  else if( nks==1&&nkl==1 && (npip==1||npim==1) )     return 5;
  assert(false);
}


tPDVector TwoKaonOnePionCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart(3);
  if(imode==0) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::Kplus);
  }
  else if(imode==1) {
    extpart[0]=getParticleData(ParticleID::K0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::Kbar0);
  }
  else if(imode==2) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::K0);
  }
  else if(imode==3) {
    extpart[0]=getParticleData(ParticleID::K_S0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::K_S0);
  }
  else if(imode==4) {
    extpart[0]=getParticleData(ParticleID::K_L0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::K_L0);
  }
  else if(imode==5) {
    extpart[0]=getParticleData(ParticleID::K_S0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::K_L0);
  }
  // conjugate the particles if needed
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(extpart[ix]->CC()) extpart[ix]=extpart[ix]->CC();
    }
  }
  // return the answer
  return extpart;
}
#line 1 "./TwoPionCzyzCurrent.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoPionCzyzCurrent class.
//

#include "TwoPionCzyzCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/ResonanceHelpers.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HERWIG_INTERPOLATOR_CLASSDESC(TwoPionCzyzCurrent,double,Energy2)

TwoPionCzyzCurrent::TwoPionCzyzCurrent()
  : omegaMag_(18.7e-4), omegaPhase_(0.106),
    omegaMass_(782.4*MeV),omegaWidth_(8.33*MeV), beta_(2.148),
    nMax_(2000), eMax_(-GeV) {
  // various parameters
  rhoMag_  =  {1.,1.,0.59,0.048,0.40,0.43};
  rhoPhase_ = {0.,0.,-2.20,-2.0,-2.9,1.19}; 
  rhoMasses_ = {773.37*MeV,1490*MeV, 1870*MeV,2120*MeV,2321*MeV,2567*MeV};
  rhoWidths_ = { 147.1*MeV, 429*MeV,  357*MeV, 300*MeV, 444*MeV, 491*MeV};
  // set up for the modes in the base class
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
}

IBPtr TwoPionCzyzCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr TwoPionCzyzCurrent::fullclone() const {
  return new_ptr(*this);
}

void TwoPionCzyzCurrent::persistentOutput(PersistentOStream & os) const {
  os << beta_ << omegaWgt_ << omegaMag_ << omegaPhase_
     << ounit(omegaMass_,GeV) << ounit(omegaWidth_,GeV)
     << rhoWgt_ << rhoMag_ << rhoPhase_
     << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV)
     << ounit(mass_,GeV) << ounit(width_,GeV) << coup_
     << dh_ << ounit(hres_,GeV2) << ounit(h0_,GeV2) << nMax_
     << ounit(eMax_,GeV) << fpiRe_ << fpiIm_;
}

void TwoPionCzyzCurrent::persistentInput(PersistentIStream & is, int) {
  is >> beta_ >> omegaWgt_ >> omegaMag_ >> omegaPhase_
     >> iunit(omegaMass_,GeV) >> iunit(omegaWidth_,GeV)
     >> rhoWgt_ >> rhoMag_ >> rhoPhase_
     >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV)
     >> iunit(mass_,GeV) >> iunit(width_,GeV) >> coup_
     >> dh_ >> iunit(hres_,GeV2) >> iunit(h0_,GeV2) >> nMax_
     >> iunit(eMax_,GeV) >> fpiRe_ >> fpiIm_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TwoPionCzyzCurrent,WeakCurrent>
describeHerwigTwoPionCzyzCurrent("Herwig::TwoPionCzyzCurrent", "HwWeakCurrents.so");

void TwoPionCzyzCurrent::Init() {

  static ClassDocumentation<TwoPionCzyzCurrent> documentation
    ("The TwoPionCzyzCurrent class uses the currents from "
     "PRD 81 094014 for the weak current with two pions",
     "The current for two pions from \\cite{Czyz:2010hj} was used.",
     "%\\cite{Czyz:2010hj}\n"
     "\\bibitem{Czyz:2010hj}\n"
     "H.~Czyz, A.~Grzelinska and J.~H.~Kuhn,\n"
     "%``Narrow resonances studies with the radiative return method,''\n"
     "Phys.\\ Rev.\\ D {\\bf 81} (2010) 094014\n"
     "doi:10.1103/PhysRevD.81.094014\n"
     "[arXiv:1002.0279 [hep-ph]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.81.094014;%%\n"
     "%28 citations counted in INSPIRE as of 30 Jul 2018\n");

  static ParVector<TwoPionCzyzCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the pi pi channel",
     &TwoPionCzyzCurrent::rhoMasses_, MeV, -1, 775.8*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<TwoPionCzyzCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the pi pi channel",
     &TwoPionCzyzCurrent::rhoWidths_, MeV, -1, 150.3*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static ParVector<TwoPionCzyzCurrent,double> interfaceRhoMagnitude
    ("RhoMagnitude",
     "Magnitude of the weight of the different resonances for the pi pi channel",
     &TwoPionCzyzCurrent::rhoMag_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoPionCzyzCurrent,double> interfaceRhoPhase
    ("RhoPhase",
     "Phase of the weight of the different resonances for the pi pi channel",
     &TwoPionCzyzCurrent::rhoPhase_, -1, 0., 0, 0,
     false, false, Interface::nolimits);
  
  static Parameter<TwoPionCzyzCurrent,unsigned int> interfacenMax
    ("nMax",
     "The maximum number of resonances to include in the sum,"
     " should be approx infinity",
     &TwoPionCzyzCurrent::nMax_, 1000, 10, 10000,
     false, false, Interface::limited);
  
  static Parameter<TwoPionCzyzCurrent,double> interfacebeta
    ("beta",
     "The beta parameter for the couplings",
     &TwoPionCzyzCurrent::beta_, 2.148, 0.0, 100.,
     false, false, Interface::limited);
  
  static Parameter<TwoPionCzyzCurrent,Energy> interfaceOmegaMass
    ("OmegaMass",
     "The mass of the omega meson",
     &TwoPionCzyzCurrent::omegaMass_, MeV,782.4*MeV, 0.0*MeV, 1000.0*MeV,
     false, false, Interface::limited);
  
  static Parameter<TwoPionCzyzCurrent,Energy> interfaceOmegaWidth
    ("OmegaWidth",
     "The mass of the omega meson",
     &TwoPionCzyzCurrent::omegaWidth_, MeV, 8.33*MeV, 0.0*MeV, 500.0*MeV,
     false, false, Interface::limited);

  static Parameter<TwoPionCzyzCurrent,double> interfaceOmegaMagnitude
    ("OmegaMagnitude",
     "The magnitude of the omega couplings",
     &TwoPionCzyzCurrent::omegaMag_, 18.7e-4, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TwoPionCzyzCurrent,double> interfaceOmegaPhase
    ("OmegaPhase",
     "The magnitude of the omega couplings",
     &TwoPionCzyzCurrent::omegaPhase_, 0.106, 0.0, 2.*Constants::pi,
     false, false, Interface::limited);

}

void TwoPionCzyzCurrent::doinit() {
  WeakCurrent::doinit();
  // check consistency of parametrers
  if(rhoMasses_.size()!=rhoWidths_.size())
    throw InitException() << "Inconsistent parameters in TwoPionCzyzCurrent"
			  << "::doinit()" << Exception::abortnow;
  // weights for the rho channels
  if(rhoMag_.size()!=rhoPhase_.size()) 
    throw InitException() << "The vectors containing the weights and phase for the"
  			  << " rho channel must be the same size in "
  			  << "TwoPionCzyzCurrent::doinit()" << Exception::runerror;
  Complex rhoSum(0.);
  for(unsigned int ix=0;ix<rhoMag_.size();++ix) {
    rhoWgt_.push_back(rhoMag_[ix]*(cos(rhoPhase_[ix])+Complex(0.,1.)*sin(rhoPhase_[ix])));
    if(ix>0) rhoSum +=rhoWgt_.back();
  }
  omegaWgt_ = omegaMag_*(cos(omegaPhase_)+Complex(0.,1.)*sin(omegaPhase_));
  // set up the masses and widths of the rho resonances
  double gamB(tgamma(2.-beta_));
  Complex cwgt(0.);
  Energy mpi(getParticleData(ParticleID::piplus)->mass());
  for(unsigned int ix=0;ix<nMax_;++ix) {
    // this is gam(2-beta+n)/gam(n+1)
    if(ix>0) {
      gamB *= ((1.-beta_+double(ix)))/double(ix);
    }
    Complex c_n = tgamma(beta_-0.5) /(0.5+double(ix)) / sqrt(Constants::pi) *
      sin(Constants::pi*(beta_-1.-double(ix)))/Constants::pi*gamB;
    if(ix%2!=0) c_n *= -1.;
    // set the masses and widths
    // calc for higher resonances
    if(ix>=rhoMasses_.size()) {
      mass_ .push_back(rhoMasses_[0]*sqrt(1.+2.*double(ix)));
      width_.push_back(rhoWidths_[0]/rhoMasses_[0]*mass_.back());
    }
    // input for lower ones
    else {
      mass_ .push_back(rhoMasses_[ix]);
      width_.push_back(rhoWidths_[ix]);
      if(ix>0) cwgt += c_n;
    }
    // parameters for the gs propagators
    hres_.push_back(Resonance::Hhat(sqr(mass_.back()),mass_.back(),width_.back(),mpi,mpi));
    dh_  .push_back(Resonance::dHhatds(mass_.back(),width_.back(),mpi,mpi));
    h0_.push_back(Resonance::H(ZERO,mass_.back(),width_.back(),mpi,mpi,dh_.back(),hres_.back()));
    coup_.push_back(c_n);
  }
  // fix up the early weights
  for(unsigned int ix=1;ix<rhoMasses_.size();++ix) {
    coup_[ix] = rhoWgt_[ix]*cwgt/rhoSum;
  }
}

void TwoPionCzyzCurrent::constructInterpolators() const {
  // construct the interpolators
  Energy mpi(getParticleData(ParticleID::piplus)->mass());
  vector<Energy2> en;
  vector<double> re,im;
  Energy maxE = eMax_>ZERO ? eMax_ : 10.*GeV;
  Energy step = (maxE-2.*mpi)/nMax_;
  Energy Q = 2.*mpi;
  for(unsigned int ix=0;ix<nMax_+1;++ix) {
    Complex value = FpiRemainder(sqr(Q),mpi,mpi);
    en.push_back(sqr(Q));
    re.push_back(value.real());
    im.push_back(value.imag());
    Q+=step;
  }
  fpiRe_ = make_InterpolatorPtr(re,en,3);
  fpiIm_ = make_InterpolatorPtr(im,en,3);
}

// complete the construction of the decay mode for integration
bool TwoPionCzyzCurrent::createMode(int icharge, tcPDPtr resonance,
				    FlavourInfo flavour,
				    unsigned int imode,PhaseSpaceModePtr mode,
				    unsigned int iloc,int ires,
				    PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((imode==0 && abs(icharge)!=3) ||
     (imode>0  && icharge !=0)) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return false;
      break;
    case IsoSpin::I3One:
      if(imode==1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode==1 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero ) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       ) return false;
  // make sure that the decays are kinematically allowed
  tPDPtr part[2];
  if(imode==0) {
    part[0]=getParticleData(ParticleID::piplus);
    part[1]=getParticleData(ParticleID::pi0);
  }
  else  {
    part[0]=getParticleData(ParticleID::piplus);
    part[1]=getParticleData(ParticleID::piminus);
  }
  Energy min(part[0]->massMin()+part[1]->massMin());
  if(min>upp) return false;
  eMax_=max(upp,eMax_);
  // set up the resonances
  tPDPtr res[3];
  if(icharge==0) {
    res[0] =getParticleData(113);
    res[1] =getParticleData(100113);
    res[2] =getParticleData(30113);
  }
  else {
    res[0] =getParticleData(213);
    res[1] =getParticleData(100213);
    res[2] =getParticleData(30213);
    if(icharge==-3) {
      for(unsigned int ix=0;ix<3;++ix) {
  	if(res[ix]&&res[ix]->CC()) res[ix]=res[ix]->CC();
      }
    }
  }
  // create the channels
  for(unsigned int ix=0;ix<3;++ix) {
    if(!res[ix]) continue;
    if(resonance && resonance != res[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],ires+1,iloc+1,ires+1,iloc+2));
  }
  // reset the masses in the intergrators
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix<rhoMasses_.size()&&res[ix]) {
      mode->resetIntermediate(res[ix],rhoMasses_[ix],rhoWidths_[ix]);
    }
  }
  return true;
}

// the particles produced by the current
tPDVector TwoPionCzyzCurrent::particles(int icharge, unsigned int imode,
					     int,int) {
  tPDVector output(2);
  if(imode==0) {
    output[0]=getParticleData(ParticleID::piplus);
    output[1]=getParticleData(ParticleID::pi0);
    if(icharge==-3) {
      for(unsigned int ix=0;ix<output.size();++ix) {
	if(output[ix]->CC()) output[ix]=output[ix]->CC();
      }
    }
  }
  else {
    output[0]=getParticleData(ParticleID::piplus);
    output[1]=getParticleData(ParticleID::piminus);
  }
  return output;
}

// hadronic current
vector<LorentzPolarizationVectorE> 
TwoPionCzyzCurrent::current(tcPDPtr resonance,
			    FlavourInfo flavour,
			    const int imode, const int ichan,Energy & scale, 
			    const tPDVector & outgoing,
			    const vector<Lorentz5Momentum> & momenta,
			    DecayIntegrator::MEOption) const {
  useMe();
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IOne)
    return vector<LorentzPolarizationVectorE>();
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode==1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode==1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return vector<LorentzPolarizationVectorE>();
  // momentum difference and sum of the mesons
  Lorentz5Momentum pdiff(momenta[0]-momenta[1]);
  Lorentz5Momentum psum (momenta[0]+momenta[1]);
  psum.rescaleMass();
  scale=psum.mass();
  // mass2 of vector intermediate state
  Energy2 q2(psum.m2());
  double dot(psum*pdiff/q2);
  psum *=dot;
  // compute the form factor
  Complex FPI=Fpi(q2,imode,ichan,resonance,momenta[0].mass(),momenta[1].mass());
  // calculate the current
  pdiff -= psum;
  return vector<LorentzPolarizationVectorE>(1,FPI*pdiff);
}
   
bool TwoPionCzyzCurrent::accept(vector<int> id) {
  // check there are only two particles
  if(id.size()!=2) return false;
  // pion modes
  if((abs(id[0])==ParticleID::piplus  &&     id[1] ==ParticleID::pi0   ) ||
     (    id[0] ==ParticleID::pi0     && abs(id[1])==ParticleID::piplus))
    return true;
  else if((id[0]==ParticleID::piminus && id[1]==ParticleID::piplus) ||
	  (id[0]==ParticleID::piplus  && id[1]==ParticleID::piminus))
    return true;
  else
    return false;
}

// the decay mode
unsigned int TwoPionCzyzCurrent::decayMode(vector<int> idout) {
  unsigned int npi(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::piplus) ++npi;
  }
  if(npi==2) return 1;
  else       return 0;
}

// output the information for the database
void TwoPionCzyzCurrent::dataBaseOutput(ofstream & output,bool header,
					bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::TwoPionCzyzCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  unsigned int ix;
  for(ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<6)  output << "newdef ";
    else      output << "insert ";
    output << name() << ":RhoMasses " << ix << " " << rhoMasses_[ix]/MeV << "\n";
  }
  for(ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<6) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidths " << ix << " " << rhoWidths_[ix]/MeV << "\n";
  }
  for(ix=0;ix<rhoWgt_.size();++ix) {
    if(ix<6) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMagnitude " << ix << " " << rhoMag_[ix]   << "\n";
    if(ix<6) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoPhase "     << ix << " " << rhoPhase_[ix] << "\n";
  }
  output << "newdef " << name() << ":OmegaMass " << omegaMass_/MeV << "\n";
  output << "newdef " << name() << ":OmegaWidth " << omegaWidth_/MeV << "\n";
  output << "newdef " << name() << ":OmegaMagnitude " << omegaMag_ << "\n";
  output << "newdef " << name() << ":OmegaPhase " << omegaPhase_ << "\n";
  output << "newdef " << name() << ":nMax " << nMax_ << "\n";
  output << "newdef " << name() << ":beta " << beta_ << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

Complex TwoPionCzyzCurrent::Fpi(Energy2 q2,const int imode, const int ichan,
				tcPDPtr resonance, Energy ma, Energy mb) const {
  Complex FPI(0.);
  unsigned int imin=0, imax = 4;
  if(ichan>0) {
    imin = ichan;
    imax = ichan+1;
  }
  if(resonance) {
    switch(resonance->id()/1000) {
    case 0:
      imax = 1;
      break;
    case 100:
      imin = 1;
      imax = 2;
      break;
    case 30 :
      imin = 2;
      imax = 3;
      break;
    default:
      assert(false);
    }
  }
  for(unsigned int ix=imin;ix<imax;++ix) {
    Complex term = coup_[ix]*Resonance::BreitWignerGS(q2,mass_[ix],width_[ix],
						      ma,mb,h0_[ix],dh_[ix],hres_[ix]);
    // include rho-omega if needed
    if(ix==0&&imode!=0)
      term *= 1./(1.+omegaWgt_)*(1.+omegaWgt_*Resonance::BreitWignerFW(q2,omegaMass_,omegaWidth_));
    FPI += term;
  }
  // interpolator for the higher resonances
  if(imax==4) {
    if(!fpiRe_) constructInterpolators();
    FPI += Complex((*fpiRe_)(q2),(*fpiIm_)(q2));
  }
  // factor for cc mode
  if(imode==0)           FPI *= sqrt(2.0);
  return FPI;
}

Complex TwoPionCzyzCurrent::FpiRemainder(Energy2 q2, Energy ma, Energy mb) const {
  Complex output(0.);
  for(unsigned int ix=4;ix<coup_.size();++ix) {
    output += coup_[ix]*Resonance::BreitWignerGS(q2,mass_[ix],width_[ix],
						 ma,mb,h0_[ix],dh_[ix],hres_[ix]);
  }
  return output;
}
#line 1 "./TwoKaonCzyzCurrent.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoKaonCzyzCurrent class.
//

#include "TwoKaonCzyzCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/ResonanceHelpers.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <cmath>

using namespace Herwig;

HERWIG_INTERPOLATOR_CLASSDESC(TwoKaonCzyzCurrent,double,Energy2)

TwoKaonCzyzCurrent::TwoKaonCzyzCurrent()
// changed parameters from 1002.0279, Fit 2 
// substituted by own fit
: betaRho_(2.19680665014), betaOmega_(2.69362046884), betaPhi_(1.94518176513),
  nMax_(200), etaPhi_(1.055), gammaOmega_(0.5), gammaPhi_(0.2), mpi_(140.*MeV), eMax_(-GeV) {
  using Constants::pi;
  // rho parameter
  rhoMag_    = {1.1148916618504967,0.050374779737077324, 0.014908906283692132,0.03902475997619905,0.038341465215871416};
  rhoPhase_  = {0    ,    pi,    pi,    pi, pi};
  rhoMasses_ = {775.49*MeV,1520.6995754050117*MeV,1740.9719246639341*MeV,1992.2811314327789*MeV};
  rhoWidths_ = {149.4 *MeV,213.41728317817743*MeV, 84.12224414791908*MeV,289.9733272437917*MeV};
  // omega parameters
  omegaMag_    = {1.3653229680598022, 0.02775156567495144, 0.32497165559032715,1.3993153161869765};
  omegaPhase_  = {0   ,    pi,    pi,   0,pi};
  omegaMasses_ = {782.65*MeV,1414.4344268685891*MeV,1655.375231284883*MeV};
  omegaWidths_ = {8.490000000000001*MeV, 85.4413887755723*MeV, 160.31760444832305*MeV};
  // phi parameters
  phiMag_    = {0.965842498579515,0.002379766320723148,0.1956211640216197,0.16527771485190898};
  phiPhase_  = {0.   ,pi    ,pi   ,0. ,0.};
  phiMasses_ = {1019.4209171596993*MeV,1594.759278457624*MeV,2156.971341201067*MeV};
  phiWidths_ = {4.252653332329334*MeV, 28.741821847408196*MeV,673.7556174184005*MeV};
  // set up for the modes in the base class
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(5);
}

IBPtr TwoKaonCzyzCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr TwoKaonCzyzCurrent::fullclone() const {
  return new_ptr(*this);
}

void TwoKaonCzyzCurrent::persistentOutput(PersistentOStream & os) const {
  os << betaRho_ << betaOmega_ << betaPhi_
     << rhoWgt_ << rhoMag_ << rhoPhase_
     << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV)
     << phiWgt_ << phiMag_ << phiPhase_
     << ounit(phiMasses_,GeV) << ounit(phiWidths_,GeV)
     << ounit(mass_,GeV) << ounit(width_,GeV) << coup_
     << dh_ << ounit(hres_,GeV2) << ounit(h0_,GeV2)
     << nMax_ << etaPhi_ << gammaOmega_ << gammaPhi_ << ounit(mpi_,GeV)
     << ounit(eMax_,GeV) << fKI0Re_ << fKI0Im_ << fKI1Re_ << fKI1Im_;
}

void TwoKaonCzyzCurrent::persistentInput(PersistentIStream & is, int) {
  is >> betaRho_ >> betaOmega_ >> betaPhi_
     >> rhoWgt_ >> rhoMag_ >> rhoPhase_
     >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV)
     >> phiWgt_ >> phiMag_ >> phiPhase_
     >> iunit(phiMasses_,GeV) >> iunit(phiWidths_,GeV)
     >> iunit(mass_,GeV) >> iunit(width_,GeV) >> coup_
     >> dh_ >> iunit(hres_,GeV2) >> iunit(h0_,GeV2)
     >> nMax_ >> etaPhi_ >> gammaOmega_ >> gammaPhi_ >> iunit(mpi_,GeV)
     >> iunit(eMax_,GeV) >> fKI0Re_ >> fKI0Im_ >> fKI1Re_ >> fKI1Im_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TwoKaonCzyzCurrent,WeakCurrent>
describeHerwigTwoKaonCzyzCurrent("Herwig::TwoKaonCzyzCurrent", "HwWeakCurrents.so");

void TwoKaonCzyzCurrent::Init() {

  static ClassDocumentation<TwoKaonCzyzCurrent> documentation
    ("The TwoKaonCzyzCurrent class uses the currents from "
     "PRD 81 094014 for the weak current with two kaons",
     "The current for two kaons from \\cite{Czyz:2010hj} was used.",
     "%\\cite{Czyz:2010hj}\n"
     "\\bibitem{Czyz:2010hj}\n"
     "H.~Czyz, A.~Grzelinska and J.~H.~Kuhn,\n"
     "%``Narrow resonances studies with the radiative return method,''\n"
     "Phys.\\ Rev.\\ D {\\bf 81} (2010) 094014\n"
     "doi:10.1103/PhysRevD.81.094014\n"
     "[arXiv:1002.0279 [hep-ph]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.81.094014;%%\n"
     "%28 citations counted in INSPIRE as of 30 Jul 2018\n");

  static ParVector<TwoKaonCzyzCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::rhoMasses_, MeV, -1, 775.8*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<TwoKaonCzyzCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::rhoWidths_, MeV, -1, 150.3*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static ParVector<TwoKaonCzyzCurrent,double> interfaceRhoMagnitude
    ("RhoMagnitude",
     "Magnitude of the weight of the different resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::rhoMag_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoKaonCzyzCurrent,double> interfaceRhoPhase
    ("RhoPhase",
     "Phase of the weight of the different resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::rhoPhase_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoKaonCzyzCurrent,Energy> interfaceOmegaMasses
    ("OmegaMasses",
     "The masses of the different omega resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::omegaMasses_, MeV, -1, 775.8*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<TwoKaonCzyzCurrent,Energy> interfaceOmegaWidths
    ("OmegaWidths",
     "The widths of the different omega resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::omegaWidths_, MeV, -1, 150.3*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static ParVector<TwoKaonCzyzCurrent,double> interfaceOmegaMagnitude
    ("OmegaMagnitude",
     "Magnitude of the weight of the different resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::omegaMag_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoKaonCzyzCurrent,double> interfaceOmegaPhase
    ("OmegaPhase",
     "Phase of the weight of the different resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::omegaPhase_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoKaonCzyzCurrent,Energy> interfacePhiMasses
    ("PhiMasses",
     "The masses of the different phi resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::phiMasses_, MeV, -1, 775.8*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<TwoKaonCzyzCurrent,Energy> interfacePhiWidths
    ("PhiWidths",
     "The widths of the different phi resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::phiWidths_, MeV, -1, 150.3*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static ParVector<TwoKaonCzyzCurrent,double> interfacePhiMagnitude
    ("PhiMagnitude",
     "Magnitude of the weight of the different resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::phiMag_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoKaonCzyzCurrent,double> interfacePhiPhase
    ("PhiPhase",
     "Phase of the weight of the different resonances for the pi pi channel",
     &TwoKaonCzyzCurrent::phiPhase_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static Parameter<TwoKaonCzyzCurrent,unsigned int> interfacenMax
    ("nMax",
     "The maximum number of resonances to include in the sum,"
     " should be approx infinity",
     &TwoKaonCzyzCurrent::nMax_, 200, 10, 10000,
     false, false, Interface::limited);

  static Parameter<TwoKaonCzyzCurrent,double> interfacebetaRho
    ("betaRho",
     "The beta parameter for the rho couplings",
     &TwoKaonCzyzCurrent::betaRho_, 2.23, 0.0, 100.,
     false, false, Interface::limited);

  static Parameter<TwoKaonCzyzCurrent,double> interfacebetaOmega
    ("betaOmega",
     "The beta parameter for the rho couplings",
     &TwoKaonCzyzCurrent::betaOmega_, 2.23, 0.0, 100.,
     false, false, Interface::limited);

  static Parameter<TwoKaonCzyzCurrent,double> interfacebetaPhi
    ("betaPhi",
     "The beta parameter for the phi couplings",
     &TwoKaonCzyzCurrent::betaPhi_, 1.97, 0.0, 100.,
     false, false, Interface::limited);

  static Parameter<TwoKaonCzyzCurrent,double> interfaceEtaPhi
    ("EtaPhi",
     "The eta_phi mixing parameter",
     &TwoKaonCzyzCurrent::etaPhi_, 1.04, 0.0, 10.0,
     false, false, Interface::limited);
  
  static Parameter<TwoKaonCzyzCurrent,double> interfacegammaOmega
    ("gammaOmega",
     "The gamma parameter for the widths of omega resonances",
     &TwoKaonCzyzCurrent::gammaOmega_, 0.5, 0.0, 1.0,
     false, false, Interface::limited);
  
  static Parameter<TwoKaonCzyzCurrent,double> interfacegammaPhi
    ("gammaPhi",
     "The gamma parameter for the widths of phi resonances",
     &TwoKaonCzyzCurrent::gammaPhi_, 0.2, 0.0, 1.0,
     false, false, Interface::limited);

}

void TwoKaonCzyzCurrent::doinit() {
  WeakCurrent::doinit();
  // check consistency of parametrers
  if(rhoMasses_.size()   != rhoWidths_.size() ||
     omegaMasses_.size() != omegaWidths_.size() ||
     phiMasses_.size()   != phiWidths_.size() )
    throw InitException() << "Inconsistent parameters in TwoKaonCzyzCurrent"
			  << "::doinit()" << Exception::abortnow;
  // weights for the rho channels
  if(rhoMag_.size()!=rhoPhase_.size()) 
    throw InitException() << "The vectors containing the weights and phase for the"
  			  << " rho channel must be the same size in "
  			  << "TwoKaonCzyzCurrent::doinit()" << Exception::runerror;
  // combine mags and phase
  for(unsigned int ix=0;ix<rhoMag_.size();++ix) {
    rhoWgt_.push_back(rhoMag_[ix]*(cos(rhoPhase_[ix])+Complex(0.,1.)*sin(rhoPhase_[ix])));
  }
  for(unsigned int ix=0;ix<omegaMag_.size();++ix) {
    omegaWgt_.push_back(omegaMag_[ix]*(cos(omegaPhase_[ix])+Complex(0.,1.)*sin(omegaPhase_[ix])));
  }
  for(unsigned int ix=0;ix<phiMag_.size();++ix) {
    phiWgt_.push_back(phiMag_[ix]*(cos(phiPhase_[ix])+Complex(0.,1.)*sin(phiPhase_[ix])));
  }
  // pion mass
  mpi_ = getParticleData(211)->mass();
  // rho masses and couplings
  double gamB(std::tgamma(2.-betaRho_));
  mass_.push_back(vector<Energy>());
  width_.push_back(vector<Energy>());
  coup_.push_back(vector<Complex>());
  Complex total(0.);
  for(unsigned int ix=0;ix<nMax_;++ix) {
    // this is gam(2-beta+n)/gam(n+1)
    if(ix>0) {
      gamB *= ((1.-betaRho_+double(ix)))/double(ix);
    }
    Complex c_n = std::tgamma(betaRho_-0.5) /(0.5+double(ix)) / sqrt(Constants::pi) *
      sin(Constants::pi*(betaRho_-1.-double(ix)))/Constants::pi*gamB;
    if(ix%2!=0) c_n *= -1.;
    // couplings
    coup_[0].push_back(c_n);
    total+=c_n;
    // set the masses and widths
    // calc for higher resonances
    if(ix>=rhoMasses_.size()) {
      mass_ [0].push_back(rhoMasses_[0]*sqrt(1.+2.*double(ix)));
      width_[0].push_back(rhoWidths_[0]/rhoMasses_[0]*mass_[0].back());
    }
    // input for lower ones
    else {
      mass_ [0].push_back(rhoMasses_[ix]);
      width_[0].push_back(rhoWidths_[ix]);
    }
    // parameters for the gs propagators
    hres_.push_back(Resonance::Hhat(sqr(mass_[0].back()),
	 			       mass_[0].back(),width_[0].back(),mpi_,mpi_));
    dh_  .push_back(Resonance::dHhatds(mass_[0].back(),width_[0].back(),mpi_,mpi_));
    h0_  .push_back(Resonance::H(ZERO,mass_[0].back(),width_[0].back(),
				 mpi_,mpi_,dh_.back(),hres_.back()));
  }
  for(unsigned int ix=0;ix<rhoWgt_.size();++ix) {
    total += rhoWgt_[ix]-coup_[0][ix];
    coup_[0][ix] = rhoWgt_[ix];
  }
  coup_[0][rhoWgt_.size()] += 1. - total;
  // omega masses and couplings
  gamB = std::tgamma(2.-betaOmega_);
  mass_.push_back(vector<Energy>());
  width_.push_back(vector<Energy>());
  coup_.push_back(vector<Complex>());
  total=0.;
  for(unsigned int ix=0;ix<nMax_;++ix) {
    // this is gam(2-beta+n)/gam(n+1)
    if(ix>0) {
      gamB *= ((1.-betaOmega_+double(ix)))/double(ix);
    }
    Complex c_n = std::tgamma(betaOmega_-0.5) /(0.5+double(ix)) / sqrt(Constants::pi) *
      sin(Constants::pi*(betaOmega_-1.-double(ix)))/Constants::pi*gamB;
    if(ix%2!=0) c_n *= -1.;
    // couplings
    coup_[1].push_back(c_n);
    total+=c_n;
    // set the masses and widths
    // calc for higher resonances
    if(ix>=omegaMasses_.size()) {
      mass_ [1].push_back(omegaMasses_[0]*sqrt(1.+2.*double(ix)));
      width_[1].push_back(gammaOmega_*mass_[1].back());
    }
    // input for lower ones
    else {
      mass_ [1].push_back(omegaMasses_[ix]);
      width_[1].push_back(omegaWidths_[ix]);
    }
  }
  for(unsigned int ix=0;ix<omegaWgt_.size();++ix) {
    total += omegaWgt_[ix]-coup_[1][ix];
    coup_[1][ix] = omegaWgt_[ix];
  }
  coup_[1][omegaWgt_.size()] += 1. - total;
  // phi masses and couplings
  gamB = std::tgamma(2.-betaPhi_);
  mass_.push_back(vector<Energy>());
  width_.push_back(vector<Energy>());
  coup_.push_back(vector<Complex>());
  total=0.;
  for(unsigned int ix=0;ix<nMax_;++ix) {
    // this is gam(2-beta+n)/gam(n+1)
    if(ix>0) {
      gamB *= ((1.-betaPhi_+double(ix)))/double(ix);
    }
    Complex c_n = std::tgamma(betaPhi_-0.5) /(0.5+double(ix)) / sqrt(Constants::pi) *
      sin(Constants::pi*(betaPhi_-1.-double(ix)))/Constants::pi*gamB;
    if(ix%2!=0) c_n *= -1.;
    // couplings
    coup_[2].push_back(c_n);
    total +=c_n;
    // set the masses and widths
    // calc for higher resonances
    if(ix>=phiMasses_.size()) {
      mass_ [2].push_back(phiMasses_[0]*sqrt(1.+2.*double(ix)));
      width_[2].push_back(gammaPhi_*mass_[2].back());
    }
    // input for lower ones
    else {
      mass_ [2].push_back(phiMasses_[ix]);
      width_[2].push_back(phiWidths_[ix]);
    }
  }
  for(unsigned int ix=0;ix<phiWgt_.size();++ix) {
    total += phiWgt_[ix]-coup_[2][ix];
    coup_[2][ix] = phiWgt_[ix];
  }
  coup_[2][phiWgt_.size()] += 1. - total;
}

void TwoKaonCzyzCurrent::constructInterpolators() const {
  // construct the interpolators
  vector<Energy2> en;
  vector<double> re0,im0;
  vector<double> re1,im1;
  Energy mK = getParticleData(ParticleID::Kplus)->mass();
  Energy maxE = eMax_>ZERO ? eMax_ : 10.*GeV; 
  Energy2 step = (sqr(maxE)-sqr(2.*mK))/nMax_;
  Energy2 Q2 = sqr(2.*mK);
  for(unsigned int ix=0;ix<nMax_+1;++ix) {
    Complex value = FkaonRemainderI1(Q2);
    re1.push_back(value.real());
    im1.push_back(value.imag());
    value = FkaonRemainderI0(Q2,mK,mK);
    re0.push_back(value.real());
    im0.push_back(value.imag());
    en.push_back(Q2);
    Q2+=step;
  }
  fKI0Re_ = make_InterpolatorPtr(re0,en,3);
  fKI0Im_ = make_InterpolatorPtr(im0,en,3);
  fKI1Re_ = make_InterpolatorPtr(re1,en,3);
  fKI1Im_ = make_InterpolatorPtr(im1,en,3);
}

// complete the construction of the decay mode for integration
bool TwoKaonCzyzCurrent::createMode(int icharge, tcPDPtr resonance,
				    FlavourInfo flavour,
				    unsigned int imode,PhaseSpaceModePtr mode,
				    unsigned int iloc,int ires,
				    PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((imode==0 && abs(icharge)!=3) ||
     (imode>0  && icharge !=0)) return false; 
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne && flavour.I!=IsoSpin::IZero ) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown&&flavour.I==IsoSpin::IOne) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return false;
      break;
    case IsoSpin::I3One:
      if(imode!=0 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode!=0 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  if(imode==0 && (flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero )) return false;
  if(imode!=0 && (flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero and
		  flavour.strange != Strangeness::ssbar )) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       ) return false;
  // make sure that the decays are kinematically allowed
  tPDPtr part[2];
  if(imode==0) {
    part[0]=getParticleData(ParticleID::Kplus);
    part[1]=getParticleData(ParticleID::Kbar0);
  }
  else if(imode==1|| imode==2)  {
    part[0]=getParticleData(ParticleID::K_S0);
    part[1]=getParticleData(ParticleID::K_L0);
  }
  else {
    part[0]=getParticleData(ParticleID::Kplus);
    part[1]=getParticleData(ParticleID::Kminus);
  }
  Energy min(part[0]->massMin()+part[1]->massMin());
  if(min>upp) return false;
  eMax_=max(upp,eMax_);
  // set the resonances
  vector<tPDPtr> res;
  if(icharge==0) {
    res.push_back(getParticleData(113   ));
    res.push_back(getParticleData(100113));
    res.push_back(getParticleData(30113 ));
    res.push_back(getParticleData(   223));
    res.push_back(getParticleData(   333));
  }
  else {
    res.push_back(getParticleData(213   ));
    res.push_back(getParticleData(100213));
    res.push_back(getParticleData(30213 ));
    if(icharge==-3) {
      for(unsigned int ix=0;ix<3;++ix) {
  	if(res[ix]&&res[ix]->CC()) res[ix]=res[ix]->CC();
      }
    }
  }
  // create the channels
  for(unsigned int ix=0;ix<res.size();++ix) {
    if(!res[ix]) continue;
    if(resonance && resonance != res[ix]) continue;
    if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IOne  && ix < 3) continue;
    if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IZero && ix >=3) continue;
    PhaseSpaceChannel newChannel((PhaseSpaceChannel(phase),ires,res[ix],ires+1,iloc+1,ires+1,iloc+2));
    mode->addChannel(newChannel);
  }
  // reset the masses in the intergrators
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix<rhoMasses_.size()&&res[ix]) {
      mode->resetIntermediate(res[ix],rhoMasses_[ix],rhoWidths_[ix]);
    }
  }
  if(res.size()>3) {
    mode->resetIntermediate(res[3],omegaMasses_[0],omegaWidths_[0]);
    mode->resetIntermediate(res[4],phiMasses_  [0],  phiWidths_[0]);
  }
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector TwoKaonCzyzCurrent::particles(int icharge, unsigned int imode,
					     int,int) {
  tPDVector output(2);
  if(imode==0) {
    output[0]=getParticleData(ParticleID::Kplus);
    output[1]=getParticleData(ParticleID::K0);
  }
  else if(imode==1||imode==2) {
    output[0]=getParticleData(ParticleID::K_S0);
    output[1]=getParticleData(ParticleID::K_L0);
  }
  else {
    output[0]=getParticleData(ParticleID::Kplus );
    output[1]=getParticleData(ParticleID::Kminus);
  }
  if(icharge==-3) {
    for(unsigned int ix=0;ix<output.size();++ix) {
      if(output[ix]->CC()) output[ix]=output[ix]->CC();
    }
  }
  return output;
}

// hadronic current   
vector<LorentzPolarizationVectorE> 
TwoKaonCzyzCurrent::current(tcPDPtr resonance,
			    FlavourInfo flavour,
			    const int imode, const int ichan,Energy & scale, 
			    const tPDVector & outgoing,
			    const vector<Lorentz5Momentum> & momenta,
			    DecayIntegrator::MEOption) const {
  useMe();
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne && flavour.I!=IsoSpin::IZero )
      return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  if(flavour.I3!=IsoSpin::I3Unknown&&flavour.I==IsoSpin::IOne) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode!=0 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode!=0 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  // momentum difference and sum of the mesons
  Lorentz5Momentum pdiff(momenta[0]-momenta[1]);
  Lorentz5Momentum psum (momenta[0]+momenta[1]);
  psum.rescaleMass();
  scale=psum.mass();
  // mass2 of vector intermediate state
  Energy2 q2(psum.m2());
  double dot(psum*pdiff/q2);
  psum *=dot;
  // calculate the current
  Complex FK = Fkaon(q2,imode,ichan,
		     flavour.I,flavour.strange,resonance,
		     momenta[0].mass(),momenta[1].mass());
  // compute the current
  pdiff -= psum;
  return vector<LorentzPolarizationVectorE>(1,FK*pdiff);
}
   
bool TwoKaonCzyzCurrent::accept(vector<int> id) {
  // check there are only two particles
  if(id.size()!=2) return false;
  // pion modes
  if((id[0]==ParticleID::Kminus && id[1]==ParticleID::K0)     ||
     (id[0]==ParticleID::K0     && id[1]==ParticleID::Kminus) ||
     (id[0]==ParticleID::Kplus  && id[1]==ParticleID::Kbar0)  ||
     (id[0]==ParticleID::Kbar0  && id[1]==ParticleID::Kplus))
    return true;
  else if((id[0]==ParticleID::Kminus && id[1]==ParticleID::Kplus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::Kminus))
    return true;
  else if((id[0]==ParticleID::K_S0 && id[1]==ParticleID::K_L0) ||
	  (id[0]==ParticleID::K_L0 && id[1]==ParticleID::K_S0))
    return true;
  else
    return false;
}

// the decay mode
unsigned int TwoKaonCzyzCurrent::decayMode(vector<int> idout) {
  unsigned int nk0(0),nkp(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::Kplus) ++nkp;
    else if(abs(idout[ix])==ParticleID::K0 ||
	    idout[ix]==ParticleID::K_L0 ||idout[ix]==ParticleID::K_S0 ) ++nk0;
  }
  if(nkp==1&&nk0==1) return 0;
  else if(nkp==2)    return 3;
  else if(nk0==2)    return 1;
  else return false;
}

// output the information for the database
void TwoKaonCzyzCurrent::dataBaseOutput(ofstream & output,bool header,
					bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::TwoKaonCzyzCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  unsigned int ix;
  for(ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<3)  output << "newdef ";
    else      output << "insert ";
    output << name() << ":RhoMasses " << ix << " " << rhoMasses_[ix]/MeV << "\n";
  }
  for(ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidths " << ix << " " << rhoWidths_[ix]/MeV << "\n";
  }
  for(ix=0;ix<rhoWgt_.size();++ix) {
    if(ix<5) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMagnitude " << ix << " " << rhoMag_[ix]   << "\n";
    if(ix<5) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoPhase "     << ix << " " << rhoPhase_[ix] << "\n";
  }
  for(ix=0;ix<omegaMasses_.size();++ix) {
    if(ix<3)  output << "newdef ";
    else      output << "insert ";
    output << name() << ":OmegaMasses " << ix << " " << omegaMasses_[ix]/MeV << "\n";
  }
  for(ix=0;ix<omegaWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":OmegaWidths " << ix << " " << omegaWidths_[ix]/MeV << "\n";
  }
  for(ix=0;ix<omegaWgt_.size();++ix) {
    if(ix<5) output << "newdef ";
    else     output << "insert ";
    output << name() << ":OmegaMagnitude " << ix << " " << omegaMag_[ix]   << "\n";
    if(ix<5) output << "newdef ";
    else     output << "insert ";
    output << name() << ":OmegaPhase "     << ix << " " << omegaPhase_[ix] << "\n";
  }
  for(ix=0;ix<phiMasses_.size();++ix) {
    if(ix<2)  output << "newdef ";
    else      output << "insert ";
    output << name() << ":PhiMasses " << ix << " " << phiMasses_[ix]/MeV << "\n";
  }
  for(ix=0;ix<phiWidths_.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":PhiWidths " << ix << " " << phiWidths_[ix]/MeV << "\n";
  }
  for(ix=0;ix<phiWgt_.size();++ix) {
    if(ix<4) output << "newdef ";
    else     output << "insert ";
    output << name() << ":PhiMagnitude " << ix << " " << phiMag_[ix]   << "\n";
    if(ix<4) output << "newdef ";
    else     output << "insert ";
    output << name() << ":PhiPhase "     << ix << " " << phiPhase_[ix] << "\n";
  }
  output << "newdef " << name() << ":betaRho " << betaRho_ << "\n";
  output << "newdef " << name() << ":betaOmega " << betaOmega_ << "\n";
  output << "newdef " << name() << ":betaPhi " << betaPhi_ << "\n";
  output << "newdef " << name() << ":gammaOmega " << gammaOmega_ << "\n";
  output << "newdef " << name() << ":gammaPhi " << gammaPhi_ << "\n";
  output << "newdef " << name() << ":etaPhi " << etaPhi_ << "\n";
  output << "newdef " << name() << ":nMax " << nMax_ << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

Complex TwoKaonCzyzCurrent::Fkaon(Energy2 q2,const int imode, const int ichan,
				  IsoSpin::IsoSpin Itotal, Strangeness::Strange strange,
				  tcPDPtr resonance, Energy ma, Energy mb) const {
  unsigned int imin=0, imax = 4;
  bool on[3] = {(Itotal==IsoSpin::IUnknown || Itotal==IsoSpin::IOne),
		(Itotal==IsoSpin::IUnknown || Itotal==IsoSpin::IZero) && imode!=0,
		(Itotal==IsoSpin::IUnknown || Itotal==IsoSpin::IZero) && imode!=0};
  if(strange!=Strangeness::Unknown) {
    if(strange==Strangeness::Zero) on[2] = false;
    else if(strange==Strangeness::ssbar) on[0]=on[1]=false;
  }
  if(ichan>=0) {
    if(ichan<3) {
      on[1]=on[2]=false;
      imin = ichan;
      imax = ichan+1;
    }
    else if(ichan==3) {
      on[0]=on[2]=false;
      imin=0;
      imax=1;
    }
    else if(ichan==4) {
      on[0]=on[1]=false;
      imin=0;
      imax=1;
    }
    else
      assert(false);
  }
  if(resonance) {
    switch(resonance->id()%1000) {
    case 223:
      imin=0;
      on[0]=on[2]=false;
      break;
    case 333:
      imin=0;
      on[0]=on[1]=false;
      break;
    case 113:
      switch(resonance->id()/1000) {
      case 0:
	imin=0;
	break;
      case 100:
	imin = 1;
	break;
      case 30 :
	imin = 2;
	break;
      default :
	assert(false);
      }
      on[1]=on[2]=false;
      break;
    default:
      assert(false);
    } 
    imax = imin+1;
  }
  // calculate the form factor
  Complex FK(0.);
  for(unsigned int ix=imin;ix<imax;++ix) {
    // rho exchange
    if(on[0]) {
      Complex term = coup_[0][ix]*Resonance::BreitWignerGS(q2,mass_[0][ix],width_[0][ix],
    							   mpi_,mpi_,h0_[ix],dh_[ix],hres_[ix]);
      FK += imode!=1 ? 0.5*term : -0.5*term;
    }
    // omega exchange
    if(on[1]) {
      Complex term = coup_[1][ix]*Resonance::BreitWignerFW(q2,mass_[1][ix],width_[1][ix]);
      FK += 1./6.*term;
    }
    // phi exchange
    if(on[2]) {
      Complex term = coup_[2][ix]*Resonance::BreitWignerPWave(q2,mass_[2][ix],width_[2][ix],ma,mb);
      if(ix==0 && imode==1 ) term *=etaPhi_;
      FK += term/3.;
    }
  }
  // remainder pieces
  if(imax==4) {
    if(!fKI1Re_) constructInterpolators();
    Complex i1((*fKI1Re_)(q2),(*fKI1Im_)(q2));
    FK += imode!=1 ? i1 : -i1;
    FK += Complex((*fKI0Re_)(q2),(*fKI0Im_)(q2));
  }
  // factor for cc mode
  if(imode==0) FK *= sqrt(2.0);
  return FK;
}

Complex  TwoKaonCzyzCurrent::FkaonRemainderI1(Energy2 q2) const {
  Complex output(0.);
  for(unsigned int ix=4;ix<coup_[0].size();++ix)
    output += 0.5*coup_[0][ix]*Resonance::BreitWignerGS(q2,mass_[0][ix],width_[0][ix],
							mpi_,mpi_,h0_[ix],dh_[ix],hres_[ix]);
  return output;
}

Complex  TwoKaonCzyzCurrent::FkaonRemainderI0(Energy2 q2,Energy ma, Energy mb) const {
  Complex output(0.);
  // omega exchange
  for(unsigned int ix=4;ix<coup_[1].size();++ix)
    output += 1./6.*coup_[1][ix]*Resonance::BreitWignerFW(q2,mass_[1][ix],width_[1][ix]);
  // phi exchange
  for(unsigned int ix=4;ix<coup_[2].size();++ix)
    output += 1./3.*coup_[2][ix]*Resonance::BreitWignerPWave(q2,mass_[2][ix],width_[2][ix],ma,mb);
  return output;
}
#line 1 "./ThreePionCzyzCurrent.cc"
// -*- C++ -*-
//
// ThreePionCzyzCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreePionCzyzCurrent class.
//

#include "ThreePionCzyzCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/ResonanceHelpers.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
namespace {
static const InvEnergy3 InvGeV3 = pow<-3,1>(GeV);
}

ThreePionCzyzCurrent::ThreePionCzyzCurrent()
  : mpip_(140*MeV), mpi0_(140*MeV) {
  // parameters for I=0
  // masses and widths
  rhoMasses_ = {0.77609*GeV,1.465*GeV,1.7  *GeV};
  rhoWidths_ = {0.14446*GeV,0.31 *GeV,0.235*GeV};
  omegaMasses_ = {782.4*MeV,1375*MeV,1631*MeV};
  omegaWidths_ = {8.69 *MeV, 250*MeV, 245*MeV};
  phiMass_     = 1019.24*MeV;
  phiWidth_    = 4.14*MeV;
  // couplings
  coup_I0_ = {18.20*InvGeV3,-0.87*InvGeV3,-0.77*InvGeV3,
	      -1.12*InvGeV3,-0.72*InvGeV3,-0.59*InvGeV3};
  // parameters for I=1
  rhoMasses_I1_ = {0.77609*GeV,1.7*GeV };
  rhoWidths_I1_ = {0.14446*GeV,0.26*GeV};
  omegaMass_I1_ = 782.59*MeV;
  omegaWidth_I1_= 8.49*MeV;
  // couplings
  sigma_ = -0.1;
  GW_pre_    = 1.55/sqrt(2.)*12.924*0.266/GeV;
  g_omega_pi_pi_ = 0.185;
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(2);
}

IBPtr ThreePionCzyzCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr ThreePionCzyzCurrent::fullclone() const {
  return new_ptr(*this);
}

void ThreePionCzyzCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV)
     << ounit(mpip_,GeV) << ounit(mpi0_,GeV)
     << ounit(omegaMasses_,GeV) << ounit(omegaWidths_,GeV)
     << ounit(phiMass_,GeV) <<  ounit(phiWidth_,GeV) << ounit(coup_I0_,InvGeV3)
     << ounit(rhoMasses_I1_,GeV) << ounit(rhoWidths_I1_,GeV)
     << ounit(omegaMass_I1_,GeV) << ounit(omegaWidth_I1_,GeV)
     << sigma_ << ounit(GW_pre_,1./GeV) << g_omega_pi_pi_ << ounit(GW_,GeV);
}

void ThreePionCzyzCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV)
     >> iunit(mpip_,GeV) >> iunit(mpi0_,GeV)
     >> iunit(omegaMasses_,GeV) >> iunit(omegaWidths_,GeV)
     >> iunit(phiMass_,GeV) >>  iunit(phiWidth_,GeV) >> iunit(coup_I0_,InvGeV3)
     >> iunit(rhoMasses_I1_,GeV) >> iunit(rhoWidths_I1_,GeV)
     >> iunit(omegaMass_I1_,GeV) >> iunit(omegaWidth_I1_,GeV)
     >> sigma_ >> iunit(GW_pre_,1./GeV) >> g_omega_pi_pi_ >> iunit(GW_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ThreePionCzyzCurrent,WeakCurrent>
describeHerwigThreePionCzyzCurrent("Herwig::ThreePionCzyzCurrent",
				     "HwWeakCurrents.so");

void ThreePionCzyzCurrent::Init() {

  static ClassDocumentation<ThreePionCzyzCurrent> documentation
    ("The ThreePionCzyzCurrent class is designed to implement "
     "the three pion current for e+e- collisions from Eur.Phys.J. C47 (2006) 617-624",
     "The current from \\cite{Czyz:2005as} was used for $\\pi^+\\pi^-\\pi^0$",
     "\\bibitem{Czyz:2005as}\n"
     "H.~Czyz, A.~Grzelinska, J.~H.~Kuhn and G.~Rodrigo,\n"
     "%``Electron-positron annihilation into three pions and the radiative return,''\n"
     "Eur.\\ Phys.\\ J.\\ C {\\bf 47} (2006) 617\n"
     "doi:10.1140/epjc/s2006-02614-7\n"
     "[hep-ph/0512180].\n"
     "%%CITATION = doi:10.1140/epjc/s2006-02614-7;%%\n"
     "%32 citations counted in INSPIRE as of 01 Aug 2018\n"
     );

  static ParVector<ThreePionCzyzCurrent,Energy> interfaceRhoMassesI0
    ("RhoMassesI0",
     "The rho masses for the I=0 part of the current",
     &ThreePionCzyzCurrent::rhoMasses_, GeV, -1, 0.766*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static ParVector<ThreePionCzyzCurrent,Energy> interfaceRhoWidthsI0
    ("RhoWidthsI0",
     "The rho masses for the I=0 part of the current",
     &ThreePionCzyzCurrent::rhoWidths_, GeV, -1, 0.766*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static ParVector<ThreePionCzyzCurrent,Energy> interfaceOmegaMassesI0
    ("OmegaMassesI0",
     "The omega masses for the I=0 part of the current",
     &ThreePionCzyzCurrent::omegaMasses_, GeV, -1, 0.766*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static ParVector<ThreePionCzyzCurrent,Energy> interfaceOmegaWidthsI0
    ("OmegaWidthsI0",
     "The omega masses for the I=0 part of the current",
     &ThreePionCzyzCurrent::omegaWidths_, GeV, -1, 0.766*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
 
  static Parameter<ThreePionCzyzCurrent,Energy> interfacePhiMass
    ("PhiMass",
     "The mass of the phi meson",
     &ThreePionCzyzCurrent::phiMass_, GeV, 1.0*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static Parameter<ThreePionCzyzCurrent,Energy> interfacePhiWidth
    ("PhiWidth",
     "The width of the phi meson",
     &ThreePionCzyzCurrent::phiWidth_, GeV, 1.0*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static ParVector<ThreePionCzyzCurrent,InvEnergy3> interfaceCouplingsI0
    ("CouplingsI0",
     "The couplings for the I=0 component",
     &ThreePionCzyzCurrent::coup_I0_, InvGeV3, -1, 1.0*InvGeV3, 0*InvGeV3, 0*InvGeV3,
     false, false, Interface::nolimits);

  static ParVector<ThreePionCzyzCurrent,Energy> interfaceRhoMassesI1
    ("RhoMassesI1",
     "The rho masses for the I=1 part of the current",
     &ThreePionCzyzCurrent::rhoMasses_I1_, GeV, -1, 0.766*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static ParVector<ThreePionCzyzCurrent,Energy> interfaceRhoWidthsI1
    ("RhoWidthsI1",
     "The rho masses for the I=0 part of the current",
     &ThreePionCzyzCurrent::rhoWidths_I1_, GeV, -1, 0.766*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
   
  static Parameter<ThreePionCzyzCurrent,Energy> interfaceOmegaMass
    ("OmegaMass",
     "The mass of the omega meson",
     &ThreePionCzyzCurrent::omegaMass_I1_, GeV, 0.78259*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static Parameter<ThreePionCzyzCurrent,Energy> interfaceOmegaWidth
    ("OmegaWidth",
     "The width of the omega meson",
     &ThreePionCzyzCurrent::omegaWidth_I1_, GeV, 0.00849*GeV, 0*GeV, 0*GeV,
     false, false, Interface::nolimits);
  
  static Parameter<ThreePionCzyzCurrent,double> interfacesigma
    ("sigma",
     "The sigma parameter for the I=1 component",
     &ThreePionCzyzCurrent::sigma_, -0.1, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<ThreePionCzyzCurrent,InvEnergy> interfaceGWPrefactor
    ("GWPrefactor",
     "The prefactor for the G omega coupling",
     &ThreePionCzyzCurrent::GW_pre_, 1./GeV,  1.55/sqrt(2.)*12.924*0.266/GeV, 0./GeV, 1e5/GeV,
     false, false, Interface::limited);

  static Parameter<ThreePionCzyzCurrent,double> interfaceg_omega_pipi
    ("g_omega_pipi",
     "The coupling of the omega meson to two pions",
     &ThreePionCzyzCurrent::g_omega_pi_pi_, 0.185, 0.0, 1.0,
     false, false, Interface::limited);

}

void ThreePionCzyzCurrent::doinit() {
  WeakCurrent::doinit();
  GW_ = GW_pre_*sqr(rhoMasses_I1_[0])*g_omega_pi_pi_;
  mpip_ = getParticleData(211)->mass();
  mpi0_ = getParticleData(111)->mass();
}

// complete the construction of the decay mode for integration
bool ThreePionCzyzCurrent::createMode(int icharge, tcPDPtr resonance, FlavourInfo flavour,
				      unsigned int imode,PhaseSpaceModePtr mode,
				      unsigned int iloc,int ires,
				      PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(imode>=2 || icharge != 0) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I==IsoSpin::IZero) {
      if(flavour.I3!=IsoSpin::I3Zero) return false;
    }
    else
      return false;
  }
  if(flavour.strange != Strangeness::Unknown)
     if(flavour.strange != Strangeness::Zero and flavour.strange != Strangeness::ssbar) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
  // check the kinematics
  tPDPtr pip = getParticleData(ParticleID::piplus);
  tPDPtr pim = getParticleData(ParticleID::piminus);
  tPDPtr pi0 = getParticleData(ParticleID::pi0);
  if(2*pip->mass()+pi0->mass()>upp) return false;
  // resonaces we need
  tPDPtr omega[4] = {getParticleData( 223),getParticleData( 100223),getParticleData( 30223),
		     getParticleData( 333)};
  tPDPtr rho0[3]  = {getParticleData( 113),getParticleData( 100113),getParticleData( 30113)};
  tPDPtr rhop[3]  = {getParticleData( 213),getParticleData( 100213),getParticleData( 30213)};
  tPDPtr rhom[3]  = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  // omega/omega -> rho pi
  unsigned int imin(0),imax(4);
  if(flavour.strange != Strangeness::Unknown) {
    if     (flavour.strange == Strangeness::Zero ) imax=3;
    else if(flavour.strange == Strangeness::ssbar) imin=3;
  }
  for(unsigned int ix=imin;ix<imax;++ix) {
    if(resonance && resonance != omega[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,omega[ix],
		      ires+1,rhom[0],ires+1,iloc+1,
		      ires+2,iloc+2,ires+2,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,omega[ix],
		      ires+1,rhop[0],ires+1,iloc+2,
		      ires+2,iloc+1,ires+2,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,omega[ix],
		      ires+1,rho0[0],ires+1,iloc+3,
		      ires+2,iloc+1,ires+2,iloc+2));
  }
  // phi rho 1450
  if(!resonance || resonance ==omega[3]) {
    mode->addChannel((PhaseSpaceChannel(phase),ires,omega[3],
		      ires+1,rhom[1],ires+1,iloc+1,
		      ires+2,iloc+2,ires+2,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,omega[3],
		      ires+1,rhop[1],ires+1,iloc+2,
		      ires+2,iloc+1,ires+2,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,omega[3],
		      ires+1,rho0[1],ires+1,iloc+3,
		      ires+2,iloc+1,ires+2,iloc+2));
  }
  // // omega 1650 rho 1700
  if(!resonance || resonance ==omega[2]) {
    mode->addChannel((PhaseSpaceChannel(phase),ires,omega[2],
		      ires+1,rhom[2],ires+1,iloc+1,
		      ires+2,iloc+2,ires+2,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,omega[2],
		      ires+1,rhop[2],ires+1,iloc+2,
		      ires+2,iloc+1,ires+2,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,omega[2],
		      ires+1,rho0[2],ires+1,iloc+3,
		      ires+2,iloc+1,ires+2,iloc+2));
  }
  // reset the masses in the intergrators
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix<rhoMasses_.size()) {
      if(rho0[ix])
	mode->resetIntermediate(rho0[ix],rhoMasses_[ix],rhoWidths_[ix]);
      if(rhop[ix])
	mode->resetIntermediate(rhop[ix],rhoMasses_[ix],rhoWidths_[ix]);
      if(rhom[ix])
	mode->resetIntermediate(rhom[ix],rhoMasses_[ix],rhoWidths_[ix]);
    }
  }
  for(unsigned int ix=0;ix<omegaMasses_.size();++ix) {
    if(omega[ix])
      mode->resetIntermediate(omega[ix],omegaMasses_[ix],omegaWidths_[ix]);
  }
  if(omega[3]) 
    mode->resetIntermediate(omega[3],phiMass_,phiWidth_);
  return true;
}

// the particles produced by the current
tPDVector ThreePionCzyzCurrent::particles(int icharge, unsigned int,
					  int,int) {
  assert(icharge==0);
  // return the answer
  return {getParticleData(ParticleID::piplus),
          getParticleData(ParticleID::piminus),
          getParticleData(ParticleID::pi0)};
}

namespace {
  Complex HChannel(const int & irho,
		   const Energy & mass, const Energy & width, const Energy2 & sp, const Energy2 & sm,
		   const Energy2 & s0, const Energy & mp, const Energy & m0) {
    if(irho<0)
      return Resonance::H(mass,width,sp,sm,s0,mp,m0);
    else if(irho==0)
      return Resonance::BreitWignerPWave(sm,mass,width,mp,m0);
    else if(irho==1)
      return Resonance::BreitWignerPWave(sp,mass,width,mp,m0);
    else if(irho==2)
      return Resonance::BreitWignerPWave(s0,mass,width,mp,mp);
    else
      assert(false);
  }
}
// hadronic current   
vector<LorentzPolarizationVectorE> 
ThreePionCzyzCurrent::current(tcPDPtr resonance,
			      FlavourInfo flavour,
			      const int, const int ichan, Energy & scale, 
			      const tPDVector & ,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I==IsoSpin::IZero) {
      if(flavour.I3!=IsoSpin::I3Zero) return vector<LorentzPolarizationVectorE>();
    }
    else
      return vector<LorentzPolarizationVectorE>();
  }
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return vector<LorentzPolarizationVectorE>();
  int ssbar=0;
  if(flavour.strange != Strangeness::Unknown) {
    if(flavour.strange == Strangeness::Zero) ssbar=1;
    else if (flavour.strange == Strangeness::ssbar) ssbar=2;
    else assert(false);
  }
  useMe();
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q;
  for(unsigned int ix=0;ix<momenta.size();++ix) q+=momenta[ix];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 sm = (momenta[1]+momenta[2]).m2();
  Energy2 sp = (momenta[0]+momenta[2]).m2();
  Energy2 s0 = (momenta[0]+momenta[1]).m2();
  int irho=-1;
  if(ichan>=0) {
    irho = ichan%3;
  }
  // isospin zero part of the current
  complex<InvEnergy3> F_I0(ZERO);
  //  if(flavour.I==IsoSpin::IUnknown || flavour.I==IsoSpin::IOne) {
  if((flavour.I==IsoSpin::IUnknown || flavour.I==IsoSpin::IZero) && !resonance && ichan<0) {
    // compute H rho
    Complex Hrho = HChannel(irho,rhoMasses_[0],rhoWidths_[0],sp,sm,s0,mpip_,mpi0_);
    // terms in the current
    if((!resonance || resonance->id() == 223) && ichan<=2 && ssbar!=2) {
      F_I0 += Hrho*coup_I0_[0]*Resonance::BreitWignerFW(q2,omegaMasses_[0],omegaWidths_[0]);
    }
    if((!resonance || resonance->id() == 333) && (ichan<0 || (ichan>=9&&ichan<=11)) && ssbar!=1) {
      F_I0 += Hrho*coup_I0_[1]*Resonance::BreitWignerFW(q2,phiMass_       ,phiWidth_      );
    }
    if((!resonance || resonance->id() == 100223) && (ichan<0 || (ichan>=3&&ichan<=5)) && ssbar!=2) {
      F_I0 += Hrho*coup_I0_[2]*Resonance::BreitWignerFW(q2,omegaMasses_[1],omegaWidths_[1]);
    }
    if((!resonance || resonance->id() == 30223) && (ichan<0 || (ichan>=6&&ichan<=8)) && ssbar!=2) {
      F_I0 +=  Hrho*coup_I0_[3]*Resonance::BreitWignerFW(q2,omegaMasses_[2],omegaWidths_[2]);
    }
    if((!resonance || resonance->id() == 333) && (ichan<0 || (ichan>=12&&ichan<=14)) && ssbar!=1) {
      F_I0 += coup_I0_[4]*HChannel(irho,rhoMasses_[1],rhoWidths_[1],sp,sm,s0,mpip_,mpi0_)*
	Resonance::BreitWignerFW(q2,phiMass_,phiWidth_);
    }
    if((!resonance || resonance->id() == 100223) && (ichan<0 || (ichan>=15&&ichan<=17)) && ssbar!=2) {
      F_I0 += coup_I0_[5]*HChannel(irho,rhoMasses_[2],rhoWidths_[2],sp,sm,s0,mpip_,mpi0_)*
	Resonance::BreitWignerFW(q2,omegaMasses_[2],omegaWidths_[2]);
    }
  }
  // isospin = 1
  complex<InvEnergy3> F_I1(ZERO);
  //  if((flavour.I==IsoSpin::IUnknown || flavour.I==IsoSpin::IZero) && !resonance && ichan<0) {
  if((flavour.I==IsoSpin::IUnknown || flavour.I==IsoSpin::IOne)&&ssbar!=2) {
    F_I1 = GW_*
      Resonance::BreitWignerFW(q2,omegaMass_I1_,omegaWidth_I1_)/sqr(omegaMass_I1_)*
      (Resonance::BreitWignerPWave(s0,rhoMasses_I1_[0],
				   rhoWidths_I1_[0],mpip_,mpip_)/sqr(rhoMasses_I1_[0])+
       sigma_*Resonance::BreitWignerPWave(s0,rhoMasses_I1_[1],
					  rhoWidths_I1_[1],mpip_,mpip_)/sqr(rhoMasses_I1_[1]));
  }
  // the current
  LorentzPolarizationVector vect = (F_I0+F_I1)*
    Helicity::epsilon(momenta[0],
  		      momenta[1],
  		      momenta[2]);
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,q.mass()*vect);
}
   
bool ThreePionCzyzCurrent::accept(vector<int> id) {
  if(id.size()!=3){return false;}
  unsigned int npiplus(0),npi0(0),npiminus(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus) ++npiplus;
    else if(id[ix]==ParticleID::piminus) ++npiplus;
    else if(id[ix]==ParticleID::pi0)    ++npi0;
  }
  return (npiplus==1&&npiminus==1&&npi0==1);
}

// the decay mode
unsigned int ThreePionCzyzCurrent::decayMode(vector<int> ) {
  return 0;
}

// output the information for the database
void ThreePionCzyzCurrent::dataBaseOutput(ofstream & output,bool header,
					bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::ThreePionCzyzCurrent " 
  		    << name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMassesI0 " << ix << " " << rhoMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidthsI0 " << ix << " " << rhoWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<omegaMasses_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":OmegaMassesI0 " << ix << " " << omegaMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<omegaWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":OmegaWidthsI0 " << ix << " " << omegaWidths_[ix]/GeV << "\n";
  }
  output << "newdef " << name() << ":PhiMass "  << phiMass_/GeV  << "\n";
  output << "newdef " << name() << ":PhiWidth " << phiWidth_/GeV << "\n";
  for(unsigned int ix=0;ix<coup_I0_.size();++ix) {
    if(ix<6) output << "newdef ";
    else     output << "insert ";
    output << name() << ":CouplingsI0 " << ix << " " << coup_I0_[ix]*GeV*GeV2 << "\n";
  }
  
  for(unsigned int ix=0;ix<rhoMasses_I1_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMassesI1 " << ix << " " << rhoMasses_I1_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoWidths_I1_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidthsI1 " << ix << " " << rhoWidths_I1_[ix]/GeV << "\n";
  }
  output << "newdef " << name() << ":OmegaMass "  << omegaMass_I1_/GeV  << "\n";
  output << "newdef " << name() << ":OmegaWidth " << omegaWidth_I1_/GeV << "\n";
  output << "newdef " << name() << ":sigma "      << sigma_     << "\n";  
  output << "newdef " << name() << ":GWPrefactor "      << GW_pre_*GeV     << "\n";  
  output << "newdef " << name() << ":g_omega_pipi "      << g_omega_pi_pi_ << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
#line 1 "./FourPionCzyzCurrent.cc"
// -*- C++ -*-
//
// FourPionCzyzCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FourPionCzyzCurrent class.
//

#include "FourPionCzyzCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/ResonanceHelpers.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FourPionCzyzCurrent::FourPionCzyzCurrent() 
  : mpip_(140*MeV), mpi0_(140*MeV), channelMap_(6,vector<int>()) {
  // Masses and widths of the particles
  // rho (PDG for most of current)
  rhoMasses_ = {0.7755*GeV,1.459*GeV,1.72*GeV};
  rhoWidths_ = {0.1494*GeV,0.4  *GeV,0.25*GeV};
  // fitted for F_rho
  rhoMasses_Frho_ = {0.7755*GeV,1.437*GeV             ,1.738*GeV             ,2.12*GeV};
  rhoWidths_Frho_ = {0.1494*GeV,0.6784258847826656*GeV,0.8049153117715863*GeV,0.20924569673342333*GeV};
  // omega
  omegaMass_  = 0.78265*GeV;
  omegaWidth_ = 0.00849*GeV;
  // f0
  f0Mass_  = 1.35*GeV;
  f0Width_ = 0.2 *GeV;
  // a1  
  a1Mass_  = 1.23*GeV;
  a1Width_ = 0.2*GeV;
  // Coefficents for sum over \f$\rho\f$ resonances 
  beta_a1_    ={1.,-0.051864694215520334, -0.0416090742847935,-0.0018940213981381284 };
  beta_f0_    ={1.,    73673.55747104406, -26116.259838644528,     332.84652898870786};
  beta_omega_ ={1., -0.3668543468911129,0.03624673983924276 , -0.0047186283455050064 };
  beta_B_     ={1.,-0.145};
  beta_bar_   ={1.,0.08,-0.0075};  
  // couplings for the various terms
  c_a1_    = -201.7929015686851/GeV2;
  c_f0_    =  124.09577796839454/GeV2;
  c_omega_ = -1.5792052826500913/GeV;
  c_rho_   = -2.308949096570198;
  // meson meson meson couplings
  g_rho_pi_pi_    = 5.997;
  g_omega_pi_rho_ = 42.3/GeV/GeV2/GeV2;
  g_rho_gamma_    = 0.1212*GeV2;
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(6);
}

IBPtr FourPionCzyzCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr FourPionCzyzCurrent::fullclone() const {
  return new_ptr(*this);
}

void FourPionCzyzCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(rhoMasses_,GeV) <<  ounit(rhoWidths_,GeV)
     << ounit(rhoMasses_Frho_,GeV) << ounit(rhoWidths_Frho_,GeV)
     << ounit(omegaMass_,GeV) << ounit(omegaWidth_,GeV)
     << ounit(f0Mass_,GeV) << ounit(f0Width_,GeV)
     << ounit(a1Mass_,GeV) << ounit(a1Width_,GeV)
     << beta_a1_ << beta_f0_ << beta_omega_ << beta_B_ << beta_bar_
     << ounit(c_a1_,1./GeV2) << ounit(c_f0_,1./GeV2) << ounit(c_omega_,1./GeV) << c_rho_
     << g_rho_pi_pi_ << ounit(g_omega_pi_rho_,1/GeV/GeV2/GeV2) << ounit(g_rho_gamma_,GeV2)
     << ounit(mpip_,GeV) << ounit(mpi0_,GeV) << channelMap_;
}

void FourPionCzyzCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(rhoMasses_,GeV) >>  iunit(rhoWidths_,GeV)
     >> iunit(rhoMasses_Frho_,GeV) >> iunit(rhoWidths_Frho_,GeV)
     >> iunit(omegaMass_,GeV) >> iunit(omegaWidth_,GeV)
     >> iunit(f0Mass_,GeV) >> iunit(f0Width_,GeV)
     >> iunit(a1Mass_,GeV) >> iunit(a1Width_,GeV)
     >> beta_a1_ >> beta_f0_ >> beta_omega_ >> beta_B_ >> beta_bar_
     >> iunit(c_a1_,1./GeV2) >> iunit(c_f0_,1./GeV2) >> iunit(c_omega_,1./GeV) >> c_rho_
     >> g_rho_pi_pi_ >> iunit(g_omega_pi_rho_,1/GeV/GeV2/GeV2) >> iunit(g_rho_gamma_,GeV2)
     >> iunit(mpip_,GeV) >> iunit(mpi0_,GeV) >> channelMap_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FourPionCzyzCurrent,WeakCurrent>
describeHerwigFourPionCzyzCurrent("Herwig::FourPionCzyzCurrent",
				  "HwWeakCurrents.so");

void FourPionCzyzCurrent::Init() {

  static ClassDocumentation<FourPionCzyzCurrent> documentation
    ("The FourPionCzyzCurrent class is designed to implement "
     "the four pion current for e+e- collisions from Phys.Rev. D77 (2008) 114005",
     "The current from \\cite{Czyz:2008kw} was used for four pions.",
     "\\bibitem{Czyz:2008kw}\n"
     "H.~Czyz, J.~H.~Kuhn and A.~Wapienik,\n"
     "%``Four-pion production in tau decays and e+e- annihilation: An Update,''\n"
     "Phys.\\ Rev.\\ D {\\bf 77} (2008) 114005\n"
     "doi:10.1103/PhysRevD.77.114005\n"
     "[arXiv:0804.0359 [hep-ph]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.77.114005;%%\n"
     "%35 citations counted in INSPIRE as of 02 Aug 2018\n");

  static ParVector<FourPionCzyzCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the rho mesons used by default in the current",
     &FourPionCzyzCurrent::rhoMasses_, GeV, -1, 0.7755*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static ParVector<FourPionCzyzCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the rho mesons used by default in the current",
     &FourPionCzyzCurrent::rhoWidths_, GeV, -1, 0.1494*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<FourPionCzyzCurrent,Energy> interfaceRhoMassesFrho
    ("RhoMassesFrho",
     "The masses of the rho mesons used in the F_rho piece",
     &FourPionCzyzCurrent::rhoMasses_Frho_, GeV, -1, 0.7755*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static ParVector<FourPionCzyzCurrent,Energy> interfaceRhoWidthsFrho
    ("RhoWidthsFrho",
     "The widths of the rho mesons used in the F_rho piece",
     &FourPionCzyzCurrent::rhoWidths_Frho_, GeV, -1, 0.1494*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,Energy> interfaceomegaMass
    ("omegaMass",
     "The mass of the omega meson",
     &FourPionCzyzCurrent::omegaMass_, GeV, 0.78265*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,Energy> interfaceomegaWidth
    ("omegaWidth",
     "The width of the omega meson",
     &FourPionCzyzCurrent::omegaWidth_, GeV, 0.00849*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,Energy> interfaceF0Mass
    ("f0Mass",
     "The mass of the f0 meson",
     &FourPionCzyzCurrent::f0Mass_, GeV, 1.35*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,Energy> interfaceF0Width
    ("f0Width",
     "The width of the f0 meson",
     &FourPionCzyzCurrent::f0Width_, GeV, 0.2*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,Energy> interfaceA1Mass
    ("a1Mass",
     "The mass of the a1 meson",
     &FourPionCzyzCurrent::a1Mass_, GeV, 1.23*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,Energy> interfaceA1Width
    ("a1Width",
     "The width of the a1 meson",
     &FourPionCzyzCurrent::a1Width_, GeV, 0.2*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<FourPionCzyzCurrent,double> interfacebeta_a1
    ("beta_a1",
     "The coefficients for the sum over rho resonances in the a_1 term",
     &FourPionCzyzCurrent::beta_a1_, -1, 1.0, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<FourPionCzyzCurrent,double> interfacebeta_f0
    ("beta_f0",
     "The coefficients for the sum over rho resonances in the f_0 term",
     &FourPionCzyzCurrent::beta_f0_, -1, 1.0, 0, 0,
     false, false, Interface::nolimits);
  
  static ParVector<FourPionCzyzCurrent,double> interfacebeta_omega
    ("beta_omega",
     "The coefficients for the sum over rho resonances in the omega term",
     &FourPionCzyzCurrent::beta_omega_, -1, 1.0, 0, 0,
     false, false, Interface::nolimits);
  
  static ParVector<FourPionCzyzCurrent,double> interfacebeta_B
    ("beta_B",
     "The coefficients for the sum over rho resonances in the B_rho term",
     &FourPionCzyzCurrent::beta_B_, -1, 1.0, 0, 0,
     false, false, Interface::nolimits);
  
  static ParVector<FourPionCzyzCurrent,double> interfacebeta_bar
    ("beta_bar",
     "The coefficients for the sum over rho resonances in the T_rho term",
     &FourPionCzyzCurrent::beta_bar_, -1, 1.0, 0, 0,
     false, false, Interface::nolimits);

  static Parameter<FourPionCzyzCurrent,InvEnergy2> interfacec_a1
    ("c_a1",
     "The coupling for the a_1 channel",
     &FourPionCzyzCurrent::c_a1_, 1./GeV2, -255./GeV2, -1e5/GeV2, 1e5/GeV2,
     false, false, Interface::limited);

  static Parameter<FourPionCzyzCurrent,InvEnergy2> interfacec_f0
    ("c_f0",
     "The coupling for the f_0 channel",
     &FourPionCzyzCurrent::c_f0_, 1./GeV2, 64./GeV2, -1e5/GeV2, 1e5/GeV2,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,InvEnergy> interfacec_omega
    ("c_omega",
     "The coupling for the omega channel",
     &FourPionCzyzCurrent::c_omega_, 1./GeV, -1.47/GeV, -1e5/GeV, 1e5/GeV,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,double> interfacec_rho
    ("c_rho",
     "The coupling for the rho channel",
     &FourPionCzyzCurrent::c_rho_, -2.46, 0, 0,
     false, false, Interface::nolimits);

  static Parameter<FourPionCzyzCurrent,double> interfaceg_rho_pi_pi
    ("g_rho_pi_pi",
     "The coupling of rho to two pions",
     &FourPionCzyzCurrent::g_rho_pi_pi_, 5.997, 0, 0,
     false, false, Interface::nolimits);

  static Parameter<FourPionCzyzCurrent,ThePEG::Qty<std::ratio<0,1>, std::ratio<-5,1>, std::ratio<0,1> >> interfaceg_omega_pi_rho
    ("g_omega_pi_rho",
     "The coupling of omega to rho and pi",
     &FourPionCzyzCurrent::g_omega_pi_rho_, 1./GeV/GeV2/GeV2, 42.3/GeV/GeV2/GeV2, 0.0/GeV/GeV2/GeV2, 1e5/GeV/GeV2/GeV2,
     false, false, Interface::limited);
  
  static Parameter<FourPionCzyzCurrent,Energy2> interfaceg_rho_gamma
    ("g_rho_gamma",
     "The coupling of the rho to the photon",
     &FourPionCzyzCurrent::g_rho_gamma_, GeV2, 0.1212*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);
}

void FourPionCzyzCurrent::doinit() {
  WeakCurrent::doinit();
  mpip_ = getParticleData(211)->mass();
  mpi0_ = getParticleData(111)->mass();
  // test of the current for a fixed momentum configuration
  // Lorentz5Momentum q1(0.13061870567796208*GeV,-0.21736300316234394*GeV,
  // 		      0.51725254282500699*GeV,0.59167288008090657*GeV);
  // Lorentz5Momentum q2(-1.1388573867484255 *GeV,      0.37727761929037396 *GeV,      0.31336706796993302 *GeV,    1.2472979400305677*GeV);
  // Lorentz5Momentum q3(0.11806773412672231 *GeV,      0.17873024832600765   *GeV,    0.10345508827447017   *GeV,   0.27580297667647757*GeV );
  // Lorentz5Momentum q4 (7.7487017488620830E-002*GeV, 0.16118198754624435    *GeV,    6.5813182962706746E-002*GeV, 0.23620982448991124*GeV );
  // q1.rescaleMass();
  // q2.rescaleMass();
  // q2.rescaleMass();
  // q3.rescaleMass();
  // cerr << q1/GeV << " " << q1.mass()/GeV << "\n";
  // cerr << q2/GeV << " " << q2.mass()/GeV << "\n";
  // cerr << q3/GeV << " " << q3.mass()/GeV << "\n";
  // cerr << q4/GeV << " " << q4.mass()/GeV << "\n";
  // Lorentz5Momentum Q(q1+q2+q3+q4);
  // Q.rescaleMass();
 
  // LorentzVector<complex<InvEnergy> > base = baseCurrent(Q.mass2(),tcPDPtr(),-1,Q,q1,q2,q3,q4);
  // LorentzVector<complex<InvEnergy> > test( Complex(  376.35697290395467     ,  66.446392015809550     )/GeV,
  // 					   Complex( -135.73591401998152     ,  112.36912660073307     )/GeV,
  // 					   Complex(  83.215302375273723     , -54.986430577097920     )/GeV,
  // 					   Complex( -123.56434266557559     , -22.465096431505703     )/GeV);
  // cerr << "current test X " << (base.x()-test.x())/(base.x()+test.x()) << "\n";
  // cerr << "current test Y " << (base.y()-test.y())/(base.x()+test.y()) << "\n";
  // cerr << "current test Z " << (base.z()-test.z())/(base.x()+test.z()) << "\n";
  // cerr << "current test T " << (base.t()-test.t())/(base.x()+test.t()) << "\n";
}

void FourPionCzyzCurrent::createChannels(unsigned int imode,
					 int icharge,  tcPDPtr resonance,
					 unsigned int iloc,int ires,
					 tPDVector outgoing, PhaseSpaceModePtr mode,
					 PhaseSpaceChannel phase,
					 unsigned int j1, unsigned int j2,
					 unsigned int j3, unsigned int j4,
					 int & nchan) {
  tPDPtr rho0[3]  = {getParticleData( 113),getParticleData( 100113),getParticleData( 30113)};
  tPDPtr rhop[3]  = {getParticleData( 213),getParticleData( 100213),getParticleData( 30213)};
  tPDPtr rhom[3]  = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  tPDPtr omega(getParticleData(ParticleID::omega));
  tPDPtr f0(getParticleData(10221));
  tPDPtr a1p  = getParticleData(ParticleID::a_1plus);
  tPDPtr a1m  = getParticleData(ParticleID::a_1minus);
  tPDPtr a10  = getParticleData(ParticleID::a_10);
  if(icharge==3) {
    swap(rhop,rhom);
    swap(a1p,a1m);
  }
  int rhoCharge;
  tPDPtr rho;
  tPDPtr rhoin;
  // first the a1 channels
  for(unsigned int irho=0;irho<3;++irho) {
    tPDPtr rhoin = icharge==0 ? rho0[irho] : rhom[irho];
    if(resonance && rhoin!=resonance) {
      nchan+=8;
      continue;
    }
    int a1Charge;
    tPDPtr a1;
    for(unsigned int irho2=0;irho2<2;++irho2) {
      // // find the a1
      a1Charge = outgoing[j1-1]->iCharge()+outgoing[j2-1]->iCharge()+outgoing[j4-1]->iCharge();
      a1 = a1Charge==0 ? a10 : a1p;
      if(a1->iCharge()!=a1Charge) a1=a1m;
      assert(abs(a1Charge)<=3);
      // find the rho
      rhoCharge = outgoing[j1-1]->iCharge()+outgoing[j4-1]->iCharge();
      rho = rhoCharge==0 ? rho0[irho2] : rhop[irho2];
      if(rho->iCharge()!=rhoCharge) rho = rhom[irho2];
      assert(abs(rhoCharge)<=3);
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,a1,ires+1,iloc+j3,
      			ires+2,rho,ires+2,iloc+j2,ires+3,iloc+j1,ires+3,iloc+j4));
      ++nchan; channelMap_[imode].push_back(nchan);
      // find the rho
      if(imode!=4) {
	rhoCharge = outgoing[j2-1]->iCharge()+outgoing[j4-1]->iCharge();
	rho = rhoCharge==0 ? rho0[irho2] : rhop[irho2];
	if(rho->iCharge()!=rhoCharge) rho = rhom[irho2];
	assert(abs(rhoCharge)<=3);
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,a1,ires+1,iloc+j3,
			  ires+2,rho,ires+2,iloc+j1,ires+3,iloc+j2,ires+3,iloc+j4));
	++nchan; channelMap_[imode].push_back(nchan);
      }
      else ++nchan;
      // find the second a1
      a1Charge = outgoing[j1-1]->iCharge()+outgoing[j2-1]->iCharge()+outgoing[j3-1]->iCharge();
      a1 = a1Charge==0 ? a10 : a1p;
      if(a1->iCharge()!=a1Charge) a1=a1m;
      assert(abs(a1Charge)<=3);
      // find the rho
      if(imode!=4) {
	rhoCharge = outgoing[j1-1]->iCharge()+outgoing[j3-1]->iCharge();
	rho = rhoCharge==0 ? rho0[irho2] : rhop[irho2];
	if(rho->iCharge()!=rhoCharge) rho = rhom[irho2];
	assert(abs(rhoCharge)<=3);
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,a1,ires+1,iloc+j4,
			  ires+2,rho,ires+2,iloc+j2,ires+3,iloc+j1,ires+3,iloc+j3));
	++nchan; channelMap_[imode].push_back(nchan);
      }
      else
	++nchan;
      // find the rho
      if(imode!=1) {
	rhoCharge = outgoing[j2-1]->iCharge()+outgoing[j3-1]->iCharge();
	rho = rhoCharge==0 ? rho0[irho2] : rhop[irho2];
	if(rho->iCharge()!=rhoCharge) rho = rhom[irho2];
	assert(abs(rhoCharge)<=3);
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,a1,ires+1,iloc+j4,
			  ires+2,rho,ires+2,iloc+j1,ires+3,iloc+j2,ires+3,iloc+j3));
	++nchan; channelMap_[imode].push_back(nchan);
      }
      else
	++nchan;
    }
  }
  // now the f_0 channel
  for(unsigned int irho=0;irho<3;++irho) {
    tPDPtr rhoin = icharge==0 ? rho0[irho] : rhom[irho];
    if(resonance && rhoin!=resonance) continue;
    rhoCharge = outgoing[j1-1]->iCharge()+outgoing[j2-1]->iCharge();
    for(unsigned int irho2=0;irho2<3;++irho2) {
      rho= rhoCharge==0 ? rho0[irho2] : rhop[irho2];
      if(rho->iCharge()!=rhoCharge) rho = rhom[irho2];
      assert(abs(rhoCharge)<=3);
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,rho,ires+1,f0,
  			ires+2,iloc+j1,ires+2,iloc+j2,ires+3,iloc+j3,ires+3,iloc+j4));
      ++nchan; channelMap_[imode].push_back(nchan);
    }
  }
  // now the omega channels
  if(imode>=1&&imode<=3) {
    for(unsigned int irho=0;irho<3;++irho) {
      tPDPtr rhoin = icharge==0 ? rho0[irho] : rhom[irho];
      if(resonance && rhoin!=resonance) continue;
      if(imode!=1) {
	rhoCharge = outgoing[j2-1]->iCharge()+outgoing[j3-1]->iCharge();
	rho= rhoCharge==0 ? rho0[0] : rhop[0];
	if(rho->iCharge()!=rhoCharge) rho = rhom[0];
	assert(abs(rhoCharge)<=3);
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,omega,ires+1,iloc+j1,
			  ires+2,rho,ires+2,iloc+j4,ires+3,iloc+j2,ires+3,iloc+j3));
	++nchan; channelMap_[imode].push_back(nchan);
	rhoCharge = outgoing[j2-1]->iCharge()+outgoing[j4-1]->iCharge();
	rho= rhoCharge==0 ? rho0[0] : rhop[0];
	if(rho->iCharge()!=rhoCharge) rho = rhom[0];
	assert(abs(rhoCharge)<=3);
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,omega,ires+1,iloc+j1,
			  ires+2,rho,ires+2,iloc+j3,ires+3,iloc+j2,ires+3,iloc+j4));
	++nchan; channelMap_[imode].push_back(nchan);
	rhoCharge = outgoing[j3-1]->iCharge()+outgoing[j4-1]->iCharge();
	rho= rhoCharge==0 ? rho0[0] : rhop[0];
	if(rho->iCharge()!=rhoCharge) rho = rhom[0];
	assert(abs(rhoCharge)<=3);
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,omega,ires+1,iloc+j1,
			  ires+2,rho,ires+2,iloc+j2,ires+3,iloc+j3,ires+3,iloc+j4));
	++nchan; channelMap_[imode].push_back(nchan);
      }
      else nchan+=3;
      rhoCharge = outgoing[j1-1]->iCharge()+outgoing[j3-1]->iCharge();
      rho= rhoCharge==0 ? rho0[0] : rhop[0];
      if(rho->iCharge()!=rhoCharge) rho = rhom[0];
      assert(abs(rhoCharge)<=3);
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,omega,ires+1,iloc+j2,
			ires+2,rho,ires+2,iloc+j4,ires+3,iloc+j1,ires+3,iloc+j3));
      ++nchan; channelMap_[imode].push_back(nchan);
      rhoCharge = outgoing[j1-1]->iCharge()+outgoing[j4-1]->iCharge();
      rho= rhoCharge==0 ? rho0[0] : rhop[0];
      if(rho->iCharge()!=rhoCharge) rho = rhom[0];
      assert(abs(rhoCharge)<=3);
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,omega,ires+1,iloc+j2,
			ires+2,rho,ires+2,iloc+j3,ires+3,iloc+j1,ires+3,iloc+j4));
      ++nchan; channelMap_[imode].push_back(nchan);
      rhoCharge = outgoing[j3-1]->iCharge()+outgoing[j4-1]->iCharge();
      rho= rhoCharge==0 ? rho0[0] : rhop[0];
      if(rho->iCharge()!=rhoCharge) rho = rhom[0];
      assert(abs(rhoCharge)<=3);
      mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,omega,ires+1,iloc+j2,
			ires+2,rho,ires+2,iloc+j1,ires+3,iloc+j3,ires+3,iloc+j4));
      ++nchan; channelMap_[imode].push_back(nchan);
    }
  }
  else
    nchan+=18;
  // the rho channels cancel for -000 and ++--
  if(imode==0 || imode>3) {
    nchan+=16;
    return;
  }
  // rho channels
  for(unsigned int irho=0;irho<2;++irho) {
    tPDPtr rhoin = icharge==0 ? rho0[irho] : rhom[irho];
    if(resonance && rhoin!=resonance) continue;
    for(unsigned int irho1=0;irho1<2;++irho1) {
      for(unsigned int irho2=0;irho2<2;++irho2) {
  	int rho1Charge =  outgoing[j1-1]->iCharge()+outgoing[j3-1]->iCharge();
  	tPDPtr rho1= rho1Charge==0 ? rho0[irho1] : rhop[irho1];
  	if(rho1->iCharge()!=rho1Charge) rho1 = rhom[irho1];
  	assert(abs(rho1Charge)<=3);
  	int rho2Charge =  outgoing[j2-1]->iCharge()+outgoing[j4-1]->iCharge();
  	tPDPtr rho2= rho2Charge==0 ? rho0[irho2] : rhop[irho2];
  	if(rho2->iCharge()!=rho2Charge) rho2 = rhom[irho2];
  	assert(abs(rho2Charge)<=3);
  	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,rho1,ires+1,rho2,
  			  ires+2,iloc+j1,ires+2,iloc+j3,ires+3,iloc+j2,ires+3,iloc+j4));
	++nchan; channelMap_[imode].push_back(nchan);
  	assert(icharge==rho1Charge+rho2Charge);
	if(imode!=1) {
	  rho1Charge =  outgoing[j1-1]->iCharge()+outgoing[j4-1]->iCharge();
	  rho1= rho1Charge==0 ? rho0[irho1] : rhop[irho1];
	  if(rho1->iCharge()!=rho1Charge) rho1 = rhom[irho1];
	  assert(abs(rho1Charge)<=3);
	  rho2Charge =  outgoing[j2-1]->iCharge()+outgoing[j3-1]->iCharge();
	  rho2= rho2Charge==0 ? rho0[irho2] : rhop[irho2];
	  if(rho2->iCharge()!=rho2Charge) rho2 = rhom[irho2];
	  assert(abs(rho2Charge)<=3);
	  mode->addChannel((PhaseSpaceChannel(phase),ires,rhoin,ires+1,rho1,ires+1,rho2,
			    ires+2,iloc+j1,ires+2,iloc+j4,ires+3,iloc+j2,ires+3,iloc+j3));
	  ++nchan; channelMap_[imode].push_back(nchan);
	  assert(icharge==rho1Charge+rho2Charge);
	}
	else
	  ++nchan;
      }
    }
  }
}

// complete the construction of the decay mode for integration
bool FourPionCzyzCurrent::createMode(int icharge, tcPDPtr resonance,
				     FlavourInfo flavour,
				     unsigned int imode,PhaseSpaceModePtr mode,
				     unsigned int iloc,int ires,
				     PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((abs(icharge)!=3&&imode<2) ||
     (imode>=2&&icharge!=0)) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return false;
      break;
    case IsoSpin::I3One:
      if(imode==1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode==1 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  // and other flavour
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
  // check that the modes are kinematical allowed
  Energy min(ZERO);
  if(imode==0)                min =    mpip_+3.*mpi0_;
  else if(imode==1)           min = 3.*mpip_+   mpi0_;
  else if(imode==2||imode==3) min = 2.*mpip_+2.*mpi0_;
  else                        min = 4.*mpip_;
  if(min>upp) return false;
  // resonances we need
  // get the external particles
  int iq(0),ia(0);
  tPDVector outgoing = particles(icharge,imode,iq,ia);
  int nchan(-1);
  channelMap_[imode].clear();
  if(imode==0) {
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,3,4,1,2,nchan);
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,2,4,1,3,nchan);
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,2,3,1,4,nchan);
  }
  else if(imode==1) {
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,3,2,1,4,nchan);
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,3,1,2,4,nchan);
  }
  // pi0 pi0 pi+ pi-
  else if(imode==2||imode==3) {
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,3,4,1,2,nchan);
  }
  else {
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,2,4,1,3,nchan);
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,1,4,2,3,nchan);
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,2,3,1,4,nchan);
    createChannels(imode,icharge,resonance,iloc,ires,outgoing,mode,phase,1,3,2,4,nchan);
  }
  return true;
}

// the particles produced by the current
tPDVector FourPionCzyzCurrent::particles(int icharge, unsigned int imode,
					 int,int) {
  tPDVector output(4);
  tPDPtr pi0=getParticleData(ParticleID::pi0);
  tPDPtr pip=getParticleData(ParticleID::piplus);
  tPDPtr pim=pip->CC();
  if(imode==0) {
    output[0]=pim;
    output[1]=pi0;
    output[2]=pi0;
    output[3]=pi0;
  }
  else if(imode==1) {
    output[0]=pim;
    output[1]=pim;
    output[2]=pip;
    output[3]=pi0;
  }
  else if(imode==2||imode==3) {
    output[0]=pip;
    output[1]=pim;
    output[2]=pi0;
    output[3]=pi0;
  }
  else {
    output[0]=pip;
    output[1]=pip;
    output[2]=pim;
    output[3]=pim;
  }
  if(icharge==3) {
    for(unsigned int ix=0;ix<output.size();++ix) {
      if(output[ix]->CC()) output[ix]=output[ix]->CC();
    }
  }
  // return the answer
  return output;
}


// hadronic current   
vector<LorentzPolarizationVectorE> 
FourPionCzyzCurrent::current(tcPDPtr resonance,
			     FlavourInfo flavour,
			     const int imode, const int ichan,Energy & scale,
			     const tPDVector & outgoing,
			     const vector<Lorentz5Momentum> & momenta,
			     DecayIntegrator::MEOption) const {
  int icharge(0);
  for(tPDPtr out : outgoing) icharge+=out->iCharge();
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode==1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode==1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return vector<LorentzPolarizationVectorE>();
  useMe();
  // the momenta of the particles
  Lorentz5Momentum q1(momenta[0]);
  Lorentz5Momentum q2(momenta[1]);
  Lorentz5Momentum q3(momenta[2]);
  Lorentz5Momentum q4(momenta[3]);
  Lorentz5Momentum Q(q1+q2+q3+q4);
  Q.rescaleMass();
  scale = Q.mass();
  LorentzVector<complex<InvEnergy> > output;
  assert(ichan<int(channelMap_[imode].size()));
  int ichannelB = ichan<0 ? -1 : channelMap_[imode][ichan];
  if(imode==0) {
    if(ichannelB<51)
      output  += baseCurrent(Q.mass2(),resonance,ichannelB   ,Q,q3,q4,q1,q2);
    if(ichannelB<0 || (ichannelB>=67&&ichannelB<133))
       output += baseCurrent(Q.mass2(),resonance,ichannelB-67,Q,q2,q4,q1,q3);
    if(ichannelB<0 ||  ichannelB>=133) 
       output += baseCurrent(Q.mass2(),resonance,ichannelB-133,Q,q2,q3,q1,q4);
    output *= sqrt(1./3.);
  }
  else if(imode==1) {
    if(ichannelB<117)
      output += baseCurrent(Q.mass2(),resonance,ichannelB,Q,q3,q2,q1,q4);
    if(ichannelB<0||ichannelB>=67)
      output += baseCurrent(Q.mass2(),resonance,ichannelB-67,Q,q3,q1,q2,q4);
  }
  else if(imode==2||imode==3) {
    output = baseCurrent(Q.mass2(),resonance,ichannelB,Q,q3,q4,q1,q2);
  }
  else if(imode==4||imode==5) {
    if(ichannelB<67)
      output += baseCurrent(Q.mass2(),resonance,ichannelB    ,Q,q2,q4,q1,q3);
    if(ichannelB<0 || (ichannelB>=67&&ichannelB<133))
      output += baseCurrent(Q.mass2(),resonance,ichannelB-67 ,Q,q1,q4,q2,q3);
    if(ichannelB<0 || (ichannelB>=133&&ichannelB<200)) 
      output += baseCurrent(Q.mass2(),resonance,ichannelB-133,Q,q2,q3,q1,q4);
    if(ichannelB<0 ||  ichannelB>=200)
      output += baseCurrent(Q.mass2(),resonance,ichannelB-200,Q,q1,q3,q2,q4);
  }
  return  vector<LorentzPolarizationVectorE>(1,output*Q.mass2());
}
   
bool FourPionCzyzCurrent::accept(vector<int> id) {
  bool allowed(false);
  // check four products
  if(id.size()!=4){return false;}
  int npiminus=0,npiplus=0,npi0=0;
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID:: piplus)      ++npiplus;
    else if(id[ix]==ParticleID::piminus) ++npiminus;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
  }
  if(npiminus==2&&npiplus==1&&npi0==1)      allowed=true;
  else if(npiminus==1&&npi0==3)             allowed=true;
  else if(npiplus==2&&npiminus==1&&npi0==1) allowed=true;
  else if(npiplus==1&&npi0==3)              allowed=true;
  else if(npiplus==2&&npiminus==2)          allowed=true;
  else if(npiplus==1&&npiminus==1&&npi0==2) allowed=true;
  return allowed;
}

// the decay mode
unsigned int FourPionCzyzCurrent::decayMode(vector<int> idout) {
  unsigned int npi(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::piplus) ++npi;
  }
  if(npi==1) return 0;
  else if(npi==2) return 2;
  else if(npi==3) return 1;
  else return 4;
}

// output the information for the database
void FourPionCzyzCurrent::dataBaseOutput(ofstream & output,bool header,
					bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::FourPionCzyzCurrent " 
  		    << name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMasses " << ix << " " << rhoMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidths " << ix << " " << rhoWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<4) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMassesFrho " << ix << " " << rhoMasses_Frho_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<4) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidthsFrho " << ix << " " << rhoWidths_Frho_[ix]/GeV << "\n";
  }
  output << "newdef " << name() << ":omegaMass "  << omegaMass_/GeV  << "\n";
  output << "newdef " << name() << ":omegaWidth " << omegaWidth_/GeV << "\n";
  output << "newdef " << name() << ":f0Mass "  << f0Mass_/GeV  << "\n";
  output << "newdef " << name() << ":f0Width " << f0Width_/GeV << "\n";
  output << "newdef " << name() << ":a1Mass "  << a1Mass_/GeV  << "\n";
  output << "newdef " << name() << ":a1Width " << a1Width_/GeV << "\n";
  for(unsigned int ix=0;ix<beta_a1_.size();++ix) {
    if(ix<4) output << "newdef ";
    else     output << "insert ";
    output << name() << ":beta_a1 " << ix << " " << beta_a1_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<beta_f0_.size();++ix) {
    if(ix<4) output << "newdef ";
    else     output << "insert ";
    output << name() << ":beta_f0 " << ix << " " << beta_f0_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<beta_omega_.size();++ix) {
    if(ix<4) output << "newdef ";
    else     output << "insert ";
    output << name() << ":beta_omega " << ix << " " << beta_omega_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<beta_B_.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":beta_B " << ix << " " << beta_B_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<beta_bar_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":beta_bar " << ix << " " << beta_bar_[ix] << "\n";
  }
  output << "newdef " << name() << ":c_a1    " << c_a1_*GeV2   << "\n";  
  output << "newdef " << name() << ":c_f0    " << c_f0_*GeV2   << "\n";  
  output << "newdef " << name() << ":c_omega " << c_omega_*GeV << "\n";  
  output << "newdef " << name() << ":c_rho   " << c_rho_   << "\n";  
  output << "newdef " << name() << ":g_rho_pi_pi   " <<  g_rho_pi_pi_  << "\n";  
  output << "newdef " << name() << ":g_omega_pi_rho   "
	 << g_omega_pi_rho_*GeV2*GeV2*GeV   << "\n";  
  output << "newdef " << name() << ":g_rho_gamma   "
	 <<g_rho_gamma_/GeV2    << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}

LorentzVector<complex<InvEnergy> > FourPionCzyzCurrent::
baseCurrent(Energy2 Q2, tcPDPtr resonance,const int ichan,
	    const Lorentz5Momentum & Q , const Lorentz5Momentum & q1,
	    const Lorentz5Momentum & q2, const Lorentz5Momentum & q3,
	    const Lorentz5Momentum & q4) const {
  // check the resonance
  int ires0(-1);
  if(resonance) {
    switch(resonance->id()/1000) {
    case 0:
      ires0=0;
      break;
    case 100:
      ires0=1;
      break;
    case 30 :
      ires0=2;
      break;
    default:
      assert(false);
    }
  }
  // various dot products we'll need
  Energy2 m12 = sqr(q1.mass()), m22 = sqr(q2.mass());
  Energy2 m32 = sqr(q3.mass()), m42 = sqr(q4.mass());
  Energy2 Qq1 = Q*q1, Qq2 = Q*q2, Qq3 = Q*q3, Qq4 = Q*q4;
  Energy2 Qm32= Q2-2.*Qq3+m32;
  Energy2 Qm42= Q2-2.*Qq4+m42;
  Energy2 q1q2 = q1*q2, q1q3 = q1*q3, q1q4 = q1*q4;
  Energy2 q2q3 = q2*q3, q2q4 = q2*q4, q3q4 = q3*q4;
  // coefficients of the momenta
  complex<InvEnergy2> c1(ZERO),c2(ZERO),c34(ZERO),cq(ZERO),c3(ZERO),c4(ZERO);
  // first the a_1 terms from A.3 0804.0359 (N.B. sign due defns)
  // common coefficent
  if(ichan<24) {
    int ires1(-1),ires2(-1),iterm(-1);
    if(ichan>=0) {
      ires1  = ichan/8;
      ires2  = (ichan/4)%2;
      iterm = ichan%4;
    }
    if(ires0>=0 && ires1<0) ires1=ires0;
    complex<InvEnergy2> a1_fact;
    if(ires1<0) {
      a1_fact = -c_a1_*
  	Resonance::F_rho(Q2,beta_a1_,rhoMasses_Frho_,
  			 rhoWidths_Frho_,mpip_,mpip_);
    }
    else {
      if(ires0<0 || ires0==ires1) {
  	a1_fact = -c_a1_*beta_a1_[ires1]*
  	  Resonance::BreitWignerPWave(Q2,rhoMasses_Frho_[ires1],rhoWidths_Frho_[ires1],mpip_,mpip_)
  	  /std::accumulate(beta_a1_.begin(),beta_a1_.end(),0.);
      }
      else
  	a1_fact=ZERO;
    }
    // first 2 terms
    if(iterm<2) {
      Complex bw_a1_Qm3 = Resonance::BreitWignera1(Qm32,a1Mass_,a1Width_);
      // first term
      if(iterm<=0) {
  	Complex Brhoq1q4(0.);
  	if(ires2<0) {
  	  Brhoq1q4 = bw_a1_Qm3*
  	  Resonance::F_rho(m12+m42+2.*q1q4,beta_B_,rhoMasses_,rhoWidths_,mpip_,mpip_);
  	}
  	else {
  	  Brhoq1q4 = bw_a1_Qm3*beta_B_[ires2]*
  	    Resonance::BreitWignerPWave(m12+m42+2.*q1q4,rhoMasses_[ires2],rhoWidths_[ires2],mpip_,mpip_)/
  	    std::accumulate(beta_B_.begin(),beta_B_.end(),0.);
  	}
  	double dot1  = (q1q2-q2q4)/Qm32;
  	double dot1B = (Qq1-Qq4)/Q2;
  	double dot1C = Qq3/Q2*dot1;
  	// coefficients of the momenta to construct the current
  	c1 += 0.5*a1_fact*Brhoq1q4*( 3.-dot1);
  	c2 += 0.5*a1_fact*Brhoq1q4*( 1.-dot1);
  	c34+= 0.5*a1_fact*Brhoq1q4*( 1.+dot1);
  	cq += 0.5*a1_fact*Brhoq1q4*(-1.+dot1-2.*dot1B-2.*dot1C);
      }
      // second term
      if(iterm<0||iterm==1) {
      	Complex Brhoq2q4(0.);
      	if(ires2<0) {
      	  Brhoq2q4= bw_a1_Qm3*
      	    Resonance::F_rho(m12+m42+2.*q2q4,beta_B_,rhoMasses_,rhoWidths_,mpip_,mpip_);
      	}
      	else {
      	  Brhoq2q4= bw_a1_Qm3*beta_B_[ires2]*
      	    Resonance::BreitWignerPWave(m12+m42+2.*q2q4,rhoMasses_[ires2],rhoWidths_[ires2],mpip_,mpip_)/
      	    std::accumulate(beta_B_.begin(),beta_B_.end(),0.);
      	}
      	double dot2  = (q1q2-q1q4)/Qm32;
      	double dot2B = (Qq2-Qq4)/Q2;
      	double dot2C = Qq3/Q2*dot2;
      	// coefficients of the momenta to construct the current
      	// a_1 terms
      	c1 += 0.5*a1_fact*Brhoq2q4*( 1.-dot2);
      	c2 += 0.5*a1_fact*Brhoq2q4*( 3.-dot2);
      	c34+= 0.5*a1_fact*Brhoq2q4*( 1.+dot2);
      	cq += 0.5*a1_fact*Brhoq2q4*(-1.+dot2-2.*dot2B-2.*dot2C);
      }
    }
    // second 2 terms
    if(iterm<0||iterm>=2) {
      Complex bw_a1_Qm4 = Resonance::BreitWignera1(Qm42,a1Mass_,a1Width_);
      // third term
      if(iterm<0||iterm==2) {
      	Complex Brhoq1q3(0.);
      	if(ires2<0) {
      	  Brhoq1q3 = bw_a1_Qm4*
      	    Resonance::F_rho(m12+m32+2.*q1q3,beta_B_,rhoMasses_,rhoWidths_,mpip_,mpip_);
      	}
      	else {
      	  Brhoq1q3 = bw_a1_Qm4*beta_B_[ires2]*
      	    Resonance::BreitWignerPWave(m12+m32+2.*q1q3,rhoMasses_[ires2],rhoWidths_[ires2],mpip_,mpip_)/
      	    std::accumulate(beta_B_.begin(),beta_B_.end(),0.);
      	}
      	double dot3  = (q1q2-q2q3)/Qm42;
      	double dot3B = (Qq1-Qq3)/Q2;
      	double dot3C = Qq4/Q2*dot3;
      	// coefficients of the momenta to construct the current
      	// a_1 terms
      	c1 += 0.5*a1_fact*Brhoq1q3*(-3.+dot3);
      	c2 += 0.5*a1_fact*Brhoq1q3*(-1.+dot3);
      	c34+= 0.5*a1_fact*Brhoq1q3*( 1.+dot3);
      	cq += 0.5*a1_fact*Brhoq1q3*( 1.-dot3+2.*dot3B+2.*dot3C);
      }
      // fourth term
      if(iterm<0||iterm==3) {
  	Complex Brhoq2q3(0.);
  	if(ires2<0) {
  	  Brhoq2q3 = bw_a1_Qm4*
  	  Resonance::F_rho(m12+m32+2.*q2q3,beta_B_,rhoMasses_,rhoWidths_,mpip_,mpip_);
  	}
  	else {
  	  Brhoq2q3 = bw_a1_Qm4*beta_B_[ires2]*
  	  Resonance::BreitWignerPWave(m12+m32+2.*q2q3,rhoMasses_[ires2],rhoWidths_[ires2],mpip_,mpip_)/
  	    std::accumulate(beta_B_.begin(),beta_B_.end(),0.);
  	}
  	double dot4  = (q1q2-q1q3)/Qm42;
  	double dot4B = (Qq2-Qq3)/Q2;
  	double dot4C = Qq4/Q2*dot4;
  	// coefficients of the momenta to construct the current
  	// a_1 terms
  	c1 += 0.5*a1_fact*Brhoq2q3*(-1.+dot4);
  	c2 += 0.5*a1_fact*Brhoq2q3*(-3.+dot4);
  	c34+= 0.5*a1_fact*Brhoq2q3*( 1.+dot4);
  	cq += 0.5*a1_fact*Brhoq2q3*( 1.-dot4+2.*dot4B+2.*dot4C);
      }
    }
  }
  // f_0
  if(ichan<0 || (ichan>=24 && ichan<33) ) {
    int ires1(-1),ires2(-2);
    if(ichan>0) {
      ires1 = (ichan-24)/3;
      ires2 = (ichan-24)%3;
    }
    Complex rho1(0.);
    if(ires0>=0 && ires1<0) ires1=ires0;
    if(ires1<0) {
      rho1 = Resonance::F_rho(Q2,beta_f0_,rhoMasses_Frho_,rhoWidths_Frho_,mpip_,mpip_);
    }
    else {
      if(ires0<0 || ires0==ires1) {
  	rho1 = beta_f0_[ires1]*Resonance::BreitWignerPWave(Q2,rhoMasses_Frho_[ires1],rhoWidths_Frho_[ires1],mpip_,mpip_)/
  	  std::accumulate(beta_f0_.begin(),beta_f0_.end(),0.);
      }
      else
  	rho1=0.;
    }
    Complex rho2(0.);
    if(ires2<0) {
      rho2 = Resonance::F_rho(m32+m42+2.*q3q4,beta_bar_,rhoMasses_,rhoWidths_,mpip_,mpip_);
    }
    else {
      rho2 = beta_bar_[ires2]*Resonance::BreitWignerPWave(m32+m42+2.*q3q4,rhoMasses_[ires2],rhoWidths_[ires2],mpip_,mpip_)/
  	std::accumulate(beta_bar_.begin(),beta_bar_.end(),0.);
    }
    complex<InvEnergy2> f0fact = c_f0_*rho1*rho2*
      Resonance::BreitWignerSWave(m12+m22+2.*q1q2,f0Mass_,f0Width_,mpip_,mpip_);
    // add contribution to the coefficients
    c34 -= f0fact;
    cq  += f0fact*(Qq3-Qq4)/Q2;
  }
  // omega
  if(ichan<0 || (ichan>=33&&ichan<=50) ) {
    int ires1(-1),iterm(-1);
    if(ichan>0) {
      ires1 = (ichan-33)/6;
      iterm = (ichan-33)%6;
    }
    Complex rho1(0.);
    if(ires0>=0 && ires1<0) ires1=ires0;
    if(ires1<0) {
      rho1 = Resonance::F_rho(Q2,beta_omega_,rhoMasses_Frho_,rhoWidths_Frho_,mpip_,mpip_);
    }
    else {
      if(ires0<0 || ires0==ires1) {
	rho1 = beta_omega_[ires1]*Resonance::BreitWignerPWave(Q2,rhoMasses_Frho_[ires1],rhoWidths_Frho_[ires1],mpip_,mpip_)/
	  std::accumulate(beta_omega_.begin(),beta_omega_.end(),0.);
      }
      else
	rho1=0.;
    }
    complex<ThePEG::Qty<std::ratio<0,1>, std::ratio<-6,1>, std::ratio<0,1> > >
      wfact = 2.*c_omega_*g_omega_pi_rho_*g_rho_pi_pi_*rho1;
    if(iterm<3) {
      Complex Hterm(0.);
      if(iterm<0) {
	Hterm = Resonance::H(rhoMasses_[0],rhoWidths_[0],m22+2.*q2q3+m32,m22+2.*q2q4+m42,
			     m32+2.*q3q4+m42,mpip_,mpip_);
      }
      else if(iterm==0)
	Hterm = Resonance::BreitWignerPWave(m22+2.*q2q3+m32,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_);
      else if(iterm==1)
	Hterm = Resonance::BreitWignerPWave(m22+2.*q2q4+m42,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_);
      else if(iterm==2)
	Hterm = Resonance::BreitWignerPWave(m32+2.*q3q4+m42,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_);
      else
	assert(false);
      complex<ThePEG::Qty<std::ratio<0,1>, std::ratio<-6,1>, std::ratio<0,1> > >
	bw_omega_1 = wfact*Resonance::BreitWignerFW(m12-2.*Qq1+Q2,omegaMass_,omegaWidth_)*Hterm;
      c2 -=  bw_omega_1*(q1q4*Qq3-q1q3*Qq4);
      c3 -=  bw_omega_1*(q1q2*Qq4-q1q4*Qq2);
      c4 -=  bw_omega_1*(q1q3*Qq2-q1q2*Qq3);
    }
    if(iterm<0||iterm>=3) {
      Complex Hterm(0.);
      if(iterm<0) {
	Hterm = Resonance::H(rhoMasses_[0],rhoWidths_[0],m12+2.*q1q3+m32,m12+2.*q1q4+m42,
			     m32+2.*q3q4+m42,mpip_,mpip_);
      }
      else if(iterm==3)
	Hterm = Resonance::BreitWignerPWave(m12+2.*q1q3+m32,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_);
      else if(iterm==4)
	Hterm = Resonance::BreitWignerPWave(m12+2.*q1q4+m42,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_);
      else if(iterm==5)
	Hterm = Resonance::BreitWignerPWave(m32+2.*q3q4+m42,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_);
      else
	assert(false);
      complex<ThePEG::Qty<std::ratio<0,1>, std::ratio<-6,1>, std::ratio<0,1> > >
	bw_omega_2 = wfact*Resonance::BreitWignerFW(m22-2.*Qq2+Q2,omegaMass_,omegaWidth_)*Hterm;
      c1 -= bw_omega_2*(q2q4*Qq3-q2q3*Qq4);
      c3 -= bw_omega_2*(q1q2*Qq4-q2q4*Qq1);
      c4 -= bw_omega_2*(q2q3*Qq1-q1q2*Qq3);
    }
  }
  // the rho term
  LorentzVector<complex<InvEnergy> > v_rho;
  if(ichan<0||ichan>=51) {
    int ires1(-1),ires2(-1),ires3(-1),iterm(-1);
    if(ichan>0) {
      ires1 = (ichan-51)/8;
      ires2 = ((ichan-51)/4)%2;
      ires3 = ((ichan-51)/2)%2;
      iterm = (ichan-51)%2;
    }
    // prefactor
    complex<InvEnergy2> rho1;
    if(ires0>=0 && ires1<0) ires1=ires0;
    if(ires1<0) {
      rho1 = Resonance::BreitWignerDiff(Q2,rhoMasses_[0],rhoWidths_[0],
   					rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
    }
    else {
      if(ires0<0 || ires0==ires1) {
   	if(ires1==0)
   	  rho1 =  Resonance::BreitWignerPWave(Q2,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_)/sqr(rhoMasses_[0]);
   	else if(ires1==1)
   	  rho1 = -Resonance::BreitWignerPWave(Q2,rhoMasses_[1],rhoWidths_[1],mpip_,mpip_)/sqr(rhoMasses_[1]);
   	else
   	  rho1 = ZERO;
      }
      else
   	rho1 = ZERO;
    }
    Complex pre_rho = c_rho_*pow(g_rho_pi_pi_,3)*g_rho_gamma_*rho1;
    complex<InvEnergy2> BW12_q1q3_A(ZERO),BW12_q1q4_A(ZERO),BW12_q2q3_B(ZERO),BW12_q2q4_B(ZERO);
    if(ires2<0) {
      BW12_q1q3_A = Resonance::BreitWignerDiff(m12+m32+2.*q1q3,rhoMasses_[0],rhoWidths_[0],
  					       rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
      BW12_q1q4_A = Resonance::BreitWignerDiff(m12+m42+2.*q1q4,rhoMasses_[0],rhoWidths_[0],
  					       rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
      BW12_q2q3_B = Resonance::BreitWignerDiff(m22+m32+2.*q2q3,rhoMasses_[0],rhoWidths_[0],
  					       rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
      
      BW12_q2q4_B = Resonance::BreitWignerDiff(m22+m42+2.*q2q4,rhoMasses_[0],rhoWidths_[0],
  					       rhoMasses_[1],rhoWidths_[1],mpip_,mpip_);
      
    }
    else {
      if(ires2==0) {
  	BW12_q1q3_A = Resonance::BreitWignerPWave(m12+m32+2.*q1q3,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_)/sqr(rhoMasses_[0]);
  	BW12_q1q4_A = Resonance::BreitWignerPWave(m12+m42+2.*q1q4,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_)/sqr(rhoMasses_[0]);
      }
      else {
  	BW12_q1q3_A = -Resonance::BreitWignerPWave(m12+m32+2.*q1q3,rhoMasses_[1],rhoWidths_[1],mpip_,mpip_)/sqr(rhoMasses_[1]);
  	BW12_q1q4_A = -Resonance::BreitWignerPWave(m12+m42+2.*q1q4,rhoMasses_[1],rhoWidths_[1],mpip_,mpip_)/sqr(rhoMasses_[1]);
      }
    }
    if(ires3>=0) {
      if(ires3==0) {
  	BW12_q2q3_B = Resonance::BreitWignerPWave(m22+m32+2.*q2q3,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_)/sqr(rhoMasses_[0]);
  	BW12_q2q4_B = Resonance::BreitWignerPWave(m22+m42+2.*q2q4,rhoMasses_[0],rhoWidths_[0],mpip_,mpip_)/sqr(rhoMasses_[0]);
      }
      else {
  	BW12_q2q3_B = -Resonance::BreitWignerPWave(m22+m32+2.*q2q3,rhoMasses_[1],rhoWidths_[1],mpip_,mpip_)/sqr(rhoMasses_[1]);
  	BW12_q2q4_B = -Resonance::BreitWignerPWave(m22+m42+2.*q2q4,rhoMasses_[1],rhoWidths_[1],mpip_,mpip_)/sqr(rhoMasses_[1]);
      }
    }
    // now the various terms
    if(iterm<=0) {
      Energy2 d1 = Qq2 - Qq4 + 2.*q2q3 - 2.*q3q4;
      v_rho += pre_rho*BW12_q1q3_A*(BW12_q2q4_B*d1 + 2.)*q1;
    }
    if(iterm<=0) {
      Energy2 d5 = Qq1 - Qq3 + 2.*q1q2 - 2.*q2q3;
      v_rho += pre_rho*BW12_q2q4_B*(BW12_q1q3_A*d5 + 2.)*q4;
    }
    if(iterm<0||iterm==1) {
      Energy2 d2 = Qq2 - Qq3 + 2.*q2q4 - 2.*q3q4;
      v_rho -= pre_rho*BW12_q1q4_A*(BW12_q2q3_B*d2 + 2.)*q1;
    }
    if(iterm<0||iterm==1) {
      Energy2 d7 = Qq1 - Qq4 + 2.*q1q2 - 2.*q2q4;
      v_rho -= pre_rho*BW12_q2q3_B*(BW12_q1q4_A*d7 + 2.)*q3;
    }
    if(iterm<0||iterm==1) {
      Energy2 d3 = Qq1 - Qq4 + 2.*q1q3 - 2.*q3q4;
      v_rho += pre_rho*BW12_q2q3_B*(BW12_q1q4_A*d3 + 2.)*q2;
    }
    if(iterm<0||iterm==1) {
      Energy2 d6 = Qq2 - Qq3 + 2.*q1q2 - 2.*q1q3;
      v_rho += pre_rho*BW12_q1q4_A*(BW12_q2q3_B*d6 + 2.)*q4;
    }
    if(iterm<=0) {
      Energy2 d4 = Qq1 - Qq3 + 2.*q1q4 - 2.*q3q4;
      v_rho -= pre_rho*BW12_q2q4_B*(BW12_q1q3_A*d4 + 2.)*q2;
    }
    if(iterm<=0) {
      Energy2 d8 = Qq2 - Qq4 + 2.*q1q2 - 2.*q1q4;  
      v_rho -= pre_rho*BW12_q1q3_A*(BW12_q2q4_B*d8 + 2.)*q3;
    }
    complex<InvEnergy2> vdot = (Q*v_rho)/Q2;  
    v_rho  = -v_rho + vdot*Q;
  }
  // put everything together
  return c1*q1+c2*q2+(c3+c34)*q3+(c4-c34)*q4+cq*Q+v_rho;
}
#line 1 "./EtaPiPiCurrent.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiCurrent class.
//

#include "EtaPiPiCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

EtaPiPiCurrent::EtaPiPiCurrent() : fpi_(93.3*MeV) {
  rhoMasses_ = {0.77549*GeV,1.54*GeV,1.76*GeV};
  rhoWidths_ = {0.1494 *GeV,0.356*GeV,.113*GeV};
  amp_       = {1.,0.326,0.0115};
  phase_     = {0.,Constants::pi,Constants::pi};
  // set up for the modes in the base class
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
}

IBPtr EtaPiPiCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr EtaPiPiCurrent::fullclone() const {
  return new_ptr(*this);
}

void EtaPiPiCurrent::doinit() {
  WeakCurrent::doinit();
  // check consistency of parametrers
  if(rhoMasses_.size()   != rhoWidths_.size())
    throw InitException() << "Inconsistent parameters in EtaPiPiCurrent"
			  << "::doinit()" << Exception::abortnow;
  // weights for the rho channels
  if(amp_.size()!=phase_.size()) 
    throw InitException() << "The vectors containing the weights and phase for the"
  			  << " rho channel must be the same size in "
  			  << "EtaPiPiCurrent::doinit()" << Exception::runerror;
  // combine mags and phase
  weights_.clear();
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    weights_.push_back(amp_[ix]*(cos(phase_[ix])+Complex(0.,1.)*sin(phase_[ix])));
  }
  Complex denom = std::accumulate(weights_.begin(),weights_.end(),Complex(0.));
  for(unsigned int ix=0;ix<weights_.size();++ix)
    weights_[ix] /=denom;
}

void EtaPiPiCurrent::persistentOutput(PersistentOStream & os) const {
  os << weights_ << amp_ << phase_ << ounit(fpi_,MeV)
     << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV);
}

void EtaPiPiCurrent::persistentInput(PersistentIStream & is, int) {
  is >> weights_ >> amp_ >> phase_ >> iunit(fpi_,MeV)
     >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EtaPiPiCurrent,WeakCurrent>
describeHerwigEtaPiPiCurrent("Herwig::EtaPiPiCurrent", "HwWeakCurrents.so");

void EtaPiPiCurrent::Init() {

  static ClassDocumentation<EtaPiPiCurrent> documentation
    ("There is no documentation for the EtaPiPiCurrent class");

  static ParVector<EtaPiPiCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the pi pi channel",
     &EtaPiPiCurrent::rhoMasses_, MeV, -1, 775.8*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<EtaPiPiCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the pi pi channel",
     &EtaPiPiCurrent::rhoWidths_, MeV, -1, 150.3*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static ParVector<EtaPiPiCurrent,double> interfaceRhoMagnitude
    ("RhoMagnitude",
     "Magnitude of the weight of the different resonances for the pi pi channel",
     &EtaPiPiCurrent::amp_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<EtaPiPiCurrent,double> interfaceRhoPhase
    ("RhoPhase",
     "Phase of the weight of the different resonances for the pi pi channel",
     &EtaPiPiCurrent::phase_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static Parameter<EtaPiPiCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &EtaPiPiCurrent::fpi_, MeV, 93.3*MeV, ZERO, 200.0*MeV,
     false, false, true);
}

// complete the construction of the decay mode for integration
bool EtaPiPiCurrent::createMode(int icharge, tcPDPtr resonance,
				FlavourInfo flavour,
				unsigned int imode,PhaseSpaceModePtr mode,
				unsigned int iloc,int ires,
				PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((imode==0 && abs(icharge)!=3) ||
     (imode>0  && icharge !=0)) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return false;
      break;
    case IsoSpin::I3One:
      if(imode==1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode==1 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero ) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       ) return false;
  // make sure that the decays are kinematically allowed
  int iq(0),ia(0);
  tPDVector part = particles(icharge,imode,iq,ia);
  Energy min=ZERO;
  for(tPDPtr p : part) min += p->massMin();
  if(min>upp) return false;
  // set up the resonances
  tPDPtr res[3];
  if(icharge==0) {
    res[0] =getParticleData(113);
    res[1] =getParticleData(100113);
    res[2] =getParticleData(30113);
  }
  else {
    res[0] =getParticleData(213);
    res[1] =getParticleData(100213);
    res[2] =getParticleData(30213);
    if(icharge==-3) {
      for(unsigned int ix=0;ix<3;++ix) {
  	if(res[ix]&&res[ix]->CC()) res[ix]=res[ix]->CC();
      }
    }
  }
  // create the channels
  for(unsigned int ix=0;ix<3;++ix) {
    if(!res[ix]) continue;
    if(resonance && resonance != res[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],ires+1,res[0],ires+1,iloc+3,
		      ires+2,iloc+1,ires+2,iloc+2));
  }
  // reset the masses in the intergrators
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix<rhoMasses_.size()&&res[ix]) {
      mode->resetIntermediate(res[ix],rhoMasses_[ix],rhoWidths_[ix]);
    }
  }
  return true;
}

// the particles produced by the current
tPDVector EtaPiPiCurrent::particles(int icharge, unsigned int imode,
				    int,int) {
  tPDVector output(3);
  output[0]=getParticleData(ParticleID::piplus);
  output[2]=getParticleData(ParticleID::eta);
  if(imode==0) {
    output[1]=getParticleData(ParticleID::pi0);
    if(icharge==-3) {
      for(unsigned int ix=0;ix<output.size();++ix) {
	if(output[ix]->CC()) output[ix]=output[ix]->CC();
      }
    }
  }
  else {
    output[1]=getParticleData(ParticleID::piminus);
  }
  return output;
}

// hadronic current
vector<LorentzPolarizationVectorE> 
EtaPiPiCurrent::current(tcPDPtr resonance,
			FlavourInfo flavour,
			const int imode, const int ichan,Energy & scale, 
			const tPDVector & outgoing,
			const vector<Lorentz5Momentum> & momenta,
			DecayIntegrator::MEOption) const {
  useMe();
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IOne)
    return vector<LorentzPolarizationVectorE>();
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode==1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode==1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return vector<LorentzPolarizationVectorE>();
  Lorentz5Momentum q=momenta[0]+momenta[1]+momenta[2];
  q.rescaleMass();
  Energy2 s1 = (momenta[0]+momenta[1]).m2();
  Energy2 Q2 = q.mass2();
  Energy  Q  = q.mass();
  Complex BW1 = Resonance::BreitWignerPWave(s1,rhoMasses_[0],rhoWidths_[0],
					    momenta[0].mass(),momenta[1].mass());
  vector<Complex> BWs = {Resonance::BreitWignerPWave(Q2,rhoMasses_[0],rhoWidths_[0],
						     momenta[0].mass(),momenta[1].mass()),
			 Resonance::BreitWignerFW(Q2,rhoMasses_[1],
						  rhoWidths_[1]*pow(Q/rhoMasses_[1],3)),
			 Resonance::BreitWignerFW(Q2,rhoMasses_[2],
						  rhoWidths_[2]*pow(Q/rhoMasses_[2],3))};
  unsigned int imin=0,imax=3;
  if(resonance) {
    switch(resonance->id()/1000) {
    case 0:
      imax = 1;
      break;
    case 100:
      imin = 1;
      imax = 2;
      break;
    case 30 :
      imin = 2;
      imax = 3;
      break;
    default:
      assert(false);
    }
  }
  if(ichan>0) {
    imin = ichan;
    imax = ichan+1;
  }
  // form factor
  Complex fact(0.);
  for(unsigned int ix=imin;ix<imax;++ix)
    fact += weights_[ix]*BWs[ix];
  fact *= -0.25*Complex(0.,1.)/sqr(Constants::pi)/sqrt(3.)*BW1;
  if(imode==0) fact *=sqrt(2.);
  scale=Q;
  LorentzPolarizationVectorE output = fact/pow<3,1>(fpi_)*Q*
    Helicity::epsilon(momenta[0],momenta[1],momenta[2]);
  return vector<LorentzPolarizationVectorE>(1,output);
}
   
bool EtaPiPiCurrent::accept(vector<int> id) {
  // check there are only three particles
  if(id.size()!=3) return false;
  unsigned int npip(0),npim(0),npi0(0),neta(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::eta)     ++neta;
  }
  if( (npip==1&&npim==1&&neta==1) ||
      (npi0==1&&npim+npip==1&&neta==1))
    return true;
  else
    return false;
}

// the decay mode
unsigned int EtaPiPiCurrent::decayMode(vector<int> idout) {
  unsigned int npi(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::piplus) ++npi;
  }
  if(npi==2) return 1;
  else       return 0;
}

// output the information for the database
void EtaPiPiCurrent::dataBaseOutput(ofstream & output,bool header,
				    bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::EtaPiPiCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  unsigned int ix;
  for(ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<3)  output << "newdef ";
    else      output << "insert ";
    output << name() << ":RhoMasses " << ix << " " << rhoMasses_[ix]/MeV << "\n";
  }
  for(ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidths " << ix << " " << rhoWidths_[ix]/MeV << "\n";
  }
  for(ix=0;ix<weights_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMagnitude " << ix << " " << amp_[ix]   << "\n";
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoPhase "     << ix << " " << phase_[ix] << "\n";
  }
  output << "newdef " << name() << ":FPi " << fpi_/MeV << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

#line 1 "./EtaPrimePiPiCurrent.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPrimePiPiCurrent class.
//

#include "EtaPrimePiPiCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

EtaPrimePiPiCurrent::EtaPrimePiPiCurrent() : fpi_(93.3*MeV) {
  rhoMasses_ = {0.77549*GeV,1.54*GeV ,1.76*GeV,2.11*GeV};
  rhoWidths_ = {0.1494 *GeV,0.356*GeV,.113*GeV,.176*GeV};
  amp_       = {1.,0.,0.,0.02};
  phase_     = {0.,Constants::pi,Constants::pi,Constants::pi};
  // set up for the modes in the base class
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
}

IBPtr EtaPrimePiPiCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr EtaPrimePiPiCurrent::fullclone() const {
  return new_ptr(*this);
}

void EtaPrimePiPiCurrent::doinit() {
  WeakCurrent::doinit();
  // check consistency of parametrers
  if(rhoMasses_.size()   != rhoWidths_.size())
    throw InitException() << "Inconsistent parameters in EtaPrimePiPiCurrent"
			  << "::doinit()" << Exception::abortnow;
  // weights for the rho channels
  if(amp_.size()!=phase_.size()) 
    throw InitException() << "The vectors containing the weights and phase for the"
  			  << " rho channel must be the same size in "
  			  << "EtaPrimePiPiCurrent::doinit()" << Exception::runerror;
  // combine mags and phase
  weights_.clear();
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    weights_.push_back(amp_[ix]*(cos(phase_[ix])+Complex(0.,1.)*sin(phase_[ix])));
  }
  Complex denom = std::accumulate(weights_.begin(),weights_.end(),Complex(0.));
  for(unsigned int ix=0;ix<weights_.size();++ix)
    weights_[ix] /=denom;
}

void EtaPrimePiPiCurrent::persistentOutput(PersistentOStream & os) const {
  os << weights_ << amp_ << phase_ << ounit(fpi_,MeV)
     << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV);
}

void EtaPrimePiPiCurrent::persistentInput(PersistentIStream & is, int) {
  is >> weights_ >> amp_ >> phase_ >> iunit(fpi_,MeV)
     >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EtaPrimePiPiCurrent,WeakCurrent>
describeHerwigEtaPrimePiPiCurrent("Herwig::EtaPrimePiPiCurrent", "HwWeakCurrents.so");

void EtaPrimePiPiCurrent::Init() {

  static ClassDocumentation<EtaPrimePiPiCurrent> documentation
    ("There is no documentation for the EtaPrimePiPiCurrent class");

  static ParVector<EtaPrimePiPiCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the pi pi channel",
     &EtaPrimePiPiCurrent::rhoMasses_, MeV, -1, 775.8*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<EtaPrimePiPiCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the pi pi channel",
     &EtaPrimePiPiCurrent::rhoWidths_, MeV, -1, 150.3*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static ParVector<EtaPrimePiPiCurrent,double> interfaceRhoMagnitude
    ("RhoMagnitude",
     "Magnitude of the weight of the different resonances for the pi pi channel",
     &EtaPrimePiPiCurrent::amp_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<EtaPrimePiPiCurrent,double> interfaceRhoPhase
    ("RhoPhase",
     "Phase of the weight of the different resonances for the pi pi channel",
     &EtaPrimePiPiCurrent::phase_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static Parameter<EtaPrimePiPiCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &EtaPrimePiPiCurrent::fpi_, MeV, 93.3*MeV, ZERO, 200.0*MeV,
     false, false, true);
}

// complete the construction of the decay mode for integration
bool EtaPrimePiPiCurrent::createMode(int icharge, tcPDPtr resonance,
				FlavourInfo flavour,
				unsigned int imode,PhaseSpaceModePtr mode,
				unsigned int iloc,int ires,
				PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((imode==0 && abs(icharge)!=3) ||
     (imode>0  && icharge !=0)) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return false;
      break;
    case IsoSpin::I3One:
      if(imode==1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode==1 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return false;
  // make sure that the decays are kinematically allowed
  int iq(0),ia(0);
  tPDVector part = particles(icharge,imode,iq,ia);
  Energy min=ZERO;
  for(tPDPtr p : part) min += p->massMin();
  if(min>upp) return false;
  // set up the resonances
  tPDPtr res[3];
  if(icharge==0) {
    res[0] =getParticleData(113);
    res[1] =getParticleData(100113);
    res[2] =getParticleData(30113);
  }
  else {
    res[0] =getParticleData(213);
    res[1] =getParticleData(100213);
    res[2] =getParticleData(30213);
    if(icharge==-3) {
      for(unsigned int ix=0;ix<3;++ix) {
  	if(res[ix]&&res[ix]->CC()) res[ix]=res[ix]->CC();
      }
    }
  }
  // create the channels
  for(unsigned int ix=0;ix<3;++ix) {
    if(!res[ix]) continue;
    if(resonance && resonance != res[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],ires+1,res[0],ires+1,iloc+3,
		      ires+2,iloc+1,ires+2,iloc+2));
  }
  // reset the masses in the intergrators
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix<rhoMasses_.size()&&res[ix]) {
      mode->resetIntermediate(res[ix],rhoMasses_[ix],rhoWidths_[ix]);
    }
  }
  return true;
}

// the particles produced by the current
tPDVector EtaPrimePiPiCurrent::particles(int icharge, unsigned int imode,
				    int,int) {
  tPDVector output(3);
  output[0]=getParticleData(ParticleID::piplus);
  output[2]=getParticleData(ParticleID::etaprime);
  if(imode==0) {
    output[1]=getParticleData(ParticleID::pi0);
    if(icharge==-3) {
      for(unsigned int ix=0;ix<output.size();++ix) {
	if(output[ix]->CC()) output[ix]=output[ix]->CC();
      }
    }
  }
  else {
    output[1]=getParticleData(ParticleID::piminus);
  }
  return output;
}

// hadronic current
vector<LorentzPolarizationVectorE> 
EtaPrimePiPiCurrent::current(tcPDPtr resonance,
			FlavourInfo flavour,
			const int imode, const int ichan,Energy & scale, 
			const tPDVector & outgoing,
			const vector<Lorentz5Momentum> & momenta,
			DecayIntegrator::MEOption) const {
  useMe();
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IOne)
    return vector<LorentzPolarizationVectorE>();
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode==1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode==1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return vector<LorentzPolarizationVectorE>();
  Lorentz5Momentum q=momenta[0]+momenta[1]+momenta[2];
  q.rescaleMass();
  Energy2 s1 = (momenta[0]+momenta[1]).m2();
  Energy2 Q2 = q.mass2();
  Energy  Q  = q.mass();
  Complex BW1 = Resonance::BreitWignerPWave(s1,rhoMasses_[0],rhoWidths_[0],
					    momenta[0].mass(),momenta[1].mass());
  vector<Complex> BWs = {Resonance::BreitWignerPWave(Q2,rhoMasses_[0],rhoWidths_[0],
						     momenta[0].mass(),momenta[1].mass()),
			 Resonance::BreitWignerFW(Q2,rhoMasses_[1],
						  rhoWidths_[1]*pow(Q/rhoMasses_[1],3)),
			 Resonance::BreitWignerFW(Q2,rhoMasses_[2],
						  rhoWidths_[2]*pow(Q/rhoMasses_[2],3)),
			 Resonance::BreitWignerFW(Q2,rhoMasses_[3],
						  rhoWidths_[3]*pow(Q/rhoMasses_[3],3))};
  unsigned int imin=0,imax=4;
  if(resonance) {
    switch(resonance->id()/1000) {
    case 0:
      imax = 1;
      break;
    case 100:
      imin = 1;
      imax = 2;
      break;
    case 30 :
      imin = 2;
      imax = 3;
      break;
    default:
      assert(false);
    }
  }
  if(ichan>0&&ichan!=3) {
    imin = ichan;
    imax = ichan+1;
  }
  else if(ichan==3) {
    return vector<LorentzPolarizationVectorE>(1,LorentzPolarizationVectorE());
  }
  // form factor
  Complex fact(0.);
  for(unsigned int ix=imin;ix<imax;++ix) {
    fact += weights_[ix]*BWs[ix];
  }
  fact *= -0.25*Complex(0.,1.)/sqr(Constants::pi)*sqrt(2./3.)*BW1;
  if(imode==0) fact *=sqrt(2.);
  scale=Q;
  LorentzPolarizationVectorE output = fact/pow<3,1>(fpi_)*Q*
    Helicity::epsilon(momenta[0],momenta[1],momenta[2]);
  return vector<LorentzPolarizationVectorE>(1,output);
}
   
bool EtaPrimePiPiCurrent::accept(vector<int> id) {
  // check there are only three particles
  if(id.size()!=3) return false;
  unsigned int npip(0),npim(0),npi0(0),neta(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::etaprime)     ++neta;
  }
  if( (npip==1&&npim==1&&neta==1) ||
      (npi0==1&&npim+npip==1&&neta==1))
    return true;
  else
    return false;
}

// the decay mode
unsigned int EtaPrimePiPiCurrent::decayMode(vector<int> idout) {
  unsigned int npi(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::piplus) ++npi;
  }
  if(npi==2) return 1;
  else       return 0;
}

// output the information for the database
void EtaPrimePiPiCurrent::dataBaseOutput(ofstream & output,bool header,
				    bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::EtaPrimePiPiCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  unsigned int ix;
  for(ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<3)  output << "newdef ";
    else      output << "insert ";
    output << name() << ":RhoMasses " << ix << " " << rhoMasses_[ix]/MeV << "\n";
  }
  for(ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidths " << ix << " " << rhoWidths_[ix]/MeV << "\n";
  }
  for(ix=0;ix<weights_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMagnitude " << ix << " " << amp_[ix]   << "\n";
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoPhase "     << ix << " " << phase_[ix] << "\n";
  }
  output << "newdef " << name() << ":FPi " << fpi_/MeV << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

#line 1 "./KKPiCurrent.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KKPiCurrent class.
//

#include "KKPiCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace Herwig;

KKPiCurrent::KKPiCurrent() {
  // masses for the isoscalar component
  isoScalarMasses_ = {1019.461*MeV,1630*MeV,1960*MeV};
  isoScalarWidths_ = {  4.249*MeV, 218*MeV, 267*MeV};
  // masses for the isovector component
  isoVectorMasses_ = {775.26*MeV,1465*MeV,1720*MeV};
  isoVectorWidths_ = {149.1 *MeV, 400*MeV, 250*MeV};
  // amplitude and phases for the isoscalar
  isoScalarKStarAmp_  ={ZERO, 0.233/GeV, 0.0405/GeV};
  isoScalarKStarPhase_={ 0.,  1.1E-07,  5.19};
  // amplitudes and phase for the isovector component
  isoVectorKStarAmp_  ={-2.34/GeV, 0.594/GeV, -0.0179/GeV};
  isoVectorKStarPhase_={0.,0.317, 2.57};
  // Coupling for the K* to Kpi
  gKStar_ = 5.37392360229;
  // mstar masses
  mKStarP_ = 895.6*MeV;
  mKStar0_ = 895.6*MeV;
  wKStarP_ = 47.0*MeV;
  wKStar0_ = 47.0*MeV;
  // modes
  addDecayMode(3,-3);
  addDecayMode(3,-3);
  addDecayMode(3,-3);
  addDecayMode(3,-3);
  addDecayMode(3,-3);
  addDecayMode(3,-3);
}

IBPtr KKPiCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr KKPiCurrent::fullclone() const {
  return new_ptr(*this);
}

void KKPiCurrent::doinit() {
  WeakCurrent::doinit();
  static const Complex ii(0.,1.);
  assert(isoScalarKStarAmp_.size()==isoScalarKStarPhase_.size());
  for(unsigned int ix=0;ix<isoScalarKStarAmp_.size();++ix) {
    isoScalarKStarCoup_.push_back(isoScalarKStarAmp_[ix]*(cos(isoScalarKStarPhase_[ix])
							  +ii*sin(isoScalarKStarPhase_[ix])));
  }
  assert(isoVectorKStarAmp_.size()==isoVectorKStarPhase_.size());
  for(unsigned int ix=0;ix<isoVectorKStarAmp_.size();++ix)
    isoVectorKStarCoup_.push_back(isoVectorKStarAmp_[ix]*(cos(isoVectorKStarPhase_[ix])
						+ii*sin(isoVectorKStarPhase_[ix])));
}

void KKPiCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(isoScalarMasses_,GeV) << ounit(isoScalarWidths_,GeV)
     << ounit(isoVectorMasses_,GeV) << ounit(isoVectorWidths_,GeV)
     << ounit(isoScalarKStarAmp_,1./GeV) << ounit(isoVectorKStarAmp_,1./GeV)
     << isoScalarKStarPhase_ << isoVectorKStarPhase_
     << ounit(isoScalarKStarCoup_,1./GeV) << ounit(isoVectorKStarCoup_,1./GeV)
     << gKStar_
     << ounit(mKStarP_,GeV) <<  ounit(mKStar0_,GeV)
     << ounit(wKStarP_,GeV) << ounit(wKStar0_,GeV);
}

void KKPiCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(isoScalarMasses_,GeV) >> iunit(isoScalarWidths_,GeV)
     >> iunit(isoVectorMasses_,GeV) >> iunit(isoVectorWidths_,GeV)
     >> iunit(isoScalarKStarAmp_,1./GeV) >> iunit(isoVectorKStarAmp_,1./GeV)
     >> isoScalarKStarPhase_ >> isoVectorKStarPhase_
     >> iunit(isoScalarKStarCoup_,1./GeV) >> iunit(isoVectorKStarCoup_,1./GeV)
     >> gKStar_
     >> iunit(mKStarP_,GeV) >>  iunit(mKStar0_,GeV)
     >> iunit(wKStarP_,GeV) >> iunit(wKStar0_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KKPiCurrent,WeakCurrent>
describeHerwigKKPiCurrent("Herwig::KKPiCurrent", "HwWeakCurrents.so");

void KKPiCurrent::Init() {

  static ClassDocumentation<KKPiCurrent> documentation
    ("The KKPiCurrent class implements a simple model for the production"
     " of K K pi in e+e- collisions via rho and phi resonances with a"
     " K*K intermediate state.");

  static ParVector<KKPiCurrent,Energy> interfaceIsoScalarMasses
    ("IsoScalarMasses",
     "The masses for I=0 part of the current",
     &KKPiCurrent::isoScalarMasses_, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static ParVector<KKPiCurrent,Energy> interfaceIsoVectorMasses
    ("IsoVectorMasses",
     "The masses for I=1 part of the current",
     &KKPiCurrent::isoVectorMasses_, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<KKPiCurrent,Energy> interfaceIsoScalarWidths
    ("IsoScalarWidths",
     "The widths for I=0 part of the current",
     &KKPiCurrent::isoScalarWidths_, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static ParVector<KKPiCurrent,Energy> interfaceIsoVectorWidths
    ("IsoVectorWidths",
     "The widths for I=1 part of the current",
     &KKPiCurrent::isoVectorWidths_, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<KKPiCurrent,InvEnergy> interfaceIsoScalarKStarAmp
    ("IsoScalarKStarAmp",
     "Amplitude for the I=0 K* component",
     &KKPiCurrent::isoScalarKStarAmp_, 1./GeV, -1, 0.0/GeV, -1000.0/GeV, 1000.0/GeV,
     false, false, Interface::limited);
  
  static ParVector<KKPiCurrent,InvEnergy> interfaceIsoVectorKStarAmp
    ("IsoVectorKStarAmp",
     "Amplitude for the I=1 K* component",
     &KKPiCurrent::isoVectorKStarAmp_, 1./GeV, -1, 0.0/GeV, -1000.0/GeV, 1000.0/GeV,
     false, false, Interface::limited);
  
  static ParVector<KKPiCurrent,double> interfaceIsoScalarKStarPhase
    ("IsoScalarKStarPhase",
     "The phase of the couplings for the I=0 part of the current",
     &KKPiCurrent::isoScalarKStarPhase_, -1, 0., -Constants::pi, Constants::pi,
     false, false, Interface::limited);
  
  static ParVector<KKPiCurrent,double> interfaceIsoVectorKStarPhase
    ("IsoVectorKStarPhase",
     "The phase of the couplings for the I=1 part of the current",
     &KKPiCurrent::isoVectorKStarPhase_, -1, 0., -Constants::pi, Constants::pi,
     false, false, Interface::limited);

  static Parameter<KKPiCurrent,Energy> interfacemKStarPlus
    ("mKStarPlus",
     "The mass of the charged K* resonace",
     &KKPiCurrent::mKStarP_, GeV, 0.8956*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<KKPiCurrent,Energy> interfacemKStar0
    ("mKStar0",
     "The mass of the neutral K* resonace",
     &KKPiCurrent::mKStar0_, GeV, 0.8956*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<KKPiCurrent,Energy> interfacewKStarPlus
    ("wKStarPlus",
     "The width of the charged K* resonance",
     &KKPiCurrent::wKStarP_, GeV, 0.047*GeV, 0.0*GeV, 1.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<KKPiCurrent,Energy> interfacewKStar0
    ("wKStar0",
     "The width of the neutral K* resonance",
     &KKPiCurrent::wKStar0_, GeV, 0.047*GeV, 0.0*GeV, 1.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<KKPiCurrent,double> interfacegKStar
    ("gKStar",
     "The coupling of K* K pi",
     &KKPiCurrent::gKStar_, 5.37392360229, 0.0, 10.0,
     false, false, Interface::limited);

}


// complete the construction of the decay mode for integration
bool KKPiCurrent::createMode(int icharge, tcPDPtr resonance,
			     FlavourInfo flavour,
			     unsigned int imode,PhaseSpaceModePtr mode,
			     unsigned int iloc,int ires,
			     PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(icharge!=0) return false;
  if(imode>5) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I==IsoSpin::IZero) {
      if(flavour.I3!=IsoSpin::I3Zero) return false;
    }
    else if(flavour.I==IsoSpin::IOne) {
      if(flavour.I3!=IsoSpin::I3Zero) return false;
    }
    else
      return false;
  }
  if(flavour.strange != Strangeness::Unknown)
     if(flavour.strange != Strangeness::Zero and flavour.strange != Strangeness::ssbar) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  != Beauty::Zero      ) return false;
  // get the external particles
  tPDVector out = particles(0,imode,0,0);
  // check the kinematics
  Energy mout(ZERO);
  for(unsigned int ix=0;ix<out.size();++ix)
    mout += out[ix]->mass();
  if(mout>upp) return false;
  // resonances we need
  vector<tPDPtr> resI1 = {getParticleData(   113),getParticleData(100113),getParticleData( 30113)};
  vector<tPDPtr> resI0 = {getParticleData(   333),getParticleData(100333)};
  tPDPtr res[2];
  if(imode==0) {
    res[0] = getParticleData(ParticleID::Kstar0);
    res[1] = getParticleData(ParticleID::Kstarbar0);
  }
  else if(imode==1) {
    res[0] = getParticleData(ParticleID::Kstarplus);
    res[1] = getParticleData(ParticleID::Kstarminus);
  }
  else if(imode==2||imode==4) {
    res[0] = getParticleData(ParticleID::Kstarplus);
    res[1] = getParticleData(ParticleID::Kstarbar0);
  }
  else if(imode==3||imode==5) {
    res[0] = getParticleData(ParticleID::Kstarminus);
    res[1] = getParticleData(ParticleID::Kstar0);
  }
  for(unsigned int ix=0;ix<resI0.size();++ix) {
    if(resonance && resonance != resI0[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,resI0[ix],ires+1,res[0],ires+1,iloc+2,
		      ires+2,iloc+1,ires+2,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,resI0[ix],ires+1,res[1],ires+1,iloc+1,
		      ires+2,iloc+2,ires+2,iloc+3));
  }
  for(unsigned int ix=0;ix<resI1.size();++ix) {
    if(resonance && resonance != resI1[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,resI1[ix],ires+1,res[0],ires+1,iloc+2,
		      ires+2,iloc+1,ires+2,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,resI1[ix],ires+1,res[1],ires+1,iloc+1,
		      ires+2,iloc+2,ires+2,iloc+3));
  }
  return true;
}

// the particles produced by the current
tPDVector KKPiCurrent::particles(int icharge, unsigned int imode,
				 int,int) {
  assert(icharge==0);
  if(imode==0)
    return {getParticleData(ParticleID::K_S0 ),getParticleData(ParticleID::K_L0  ),getParticleData(ParticleID::pi0)};
  else if(imode==1) 
    return {getParticleData(ParticleID::Kplus),getParticleData(ParticleID::Kminus),getParticleData(ParticleID::pi0)};
  else if(imode==2)
    return {getParticleData(ParticleID::K_S0 ),getParticleData(ParticleID::Kminus),getParticleData(ParticleID::piplus)};
  else if(imode==3)
    return {getParticleData(ParticleID::K_S0 ),getParticleData(ParticleID::Kplus ),getParticleData(ParticleID::piminus)};
  else if(imode==4)
    return {getParticleData(ParticleID::K_L0 ),getParticleData(ParticleID::Kminus),getParticleData(ParticleID::piplus)};
  else if(imode==5)
    return {getParticleData(ParticleID::K_L0 ),getParticleData(ParticleID::Kplus ),getParticleData(ParticleID::piminus)};
  else
    assert(false);
}


// hadronic current   
vector<LorentzPolarizationVectorE> 
KKPiCurrent::current(tcPDPtr resonance,
		     FlavourInfo flavour,
		     const int imode, const int ichan, Energy & scale, 
		     const tPDVector & ,
		     const vector<Lorentz5Momentum> & momenta,
		     DecayIntegrator::MEOption) const {
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I==IsoSpin::IZero) {
      if(flavour.I3!=IsoSpin::I3Zero) return vector<LorentzPolarizationVectorE>();
    }
    else if(flavour.I==IsoSpin::IOne) {
      if(flavour.I3!=IsoSpin::I3Zero) return vector<LorentzPolarizationVectorE>();
    }
    else
      return vector<LorentzPolarizationVectorE>();
  }
  int ssbar=0;
  if(flavour.strange != Strangeness::Unknown) {
    if(flavour.strange == Strangeness::Zero) ssbar=1;
    else if (flavour.strange == Strangeness::ssbar) ssbar=2;
    else assert(false);
  }
  useMe();
  // calculate q2,s1,s2
  Lorentz5Momentum q;
  for(unsigned int ix=0;ix<momenta.size();++ix) q+=momenta[ix];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 s1 = (momenta[0]+momenta[2]).m2();
  Energy2 s2 = (momenta[1]+momenta[2]).m2();
  // I=0 coefficient
  complex<InvEnergy> A0(ZERO);
  int ires=-1;
  if(ichan>=0) ires=ichan/2;
  if(resonance) {
    int ires2=-1;
    switch(abs(resonance->id())) {
    case 113:
      ires2=2;
      break;
    case 100113:
      ires2=3;
      break;
    case 30113:
      ires2=4;
      break;
    case 333:
      ires2=0;
      break;
    case 100333:
      ires2=1;
      break;
    };
    if(ires>=0 && ires!=ires2) {
      return vector<LorentzPolarizationVectorE>();
    }
    ires=ires2;
  }
  if((flavour.I==IsoSpin::IUnknown || flavour.I==IsoSpin::IZero) && ssbar!=1) {
    if(ires>=0) {
      if(ires<int(isoScalarMasses_.size()))
	A0 = isoScalarKStarCoup_[ires]*Resonance::BreitWignerFW(q2,isoScalarMasses_[ires],isoScalarWidths_[ires]);
    }
    else {
      for(unsigned int ix=0;ix<isoScalarMasses_.size();++ix) {
	A0 += isoScalarKStarCoup_[ix]*Resonance::BreitWignerFW(q2,isoScalarMasses_[ix],isoScalarWidths_[ix]);
      }
    }
  }
  ires-=2;
  // I=1 coefficient
  complex<InvEnergy> A1(ZERO);
  if((flavour.I==IsoSpin::IUnknown || flavour.I==IsoSpin::IOne) && ssbar!=2) {
    if(ires>=0) {
      if(ires<int(isoVectorMasses_.size()))
	A1  = isoVectorKStarCoup_[ires]*Resonance::BreitWignerFW(q2,isoVectorMasses_[ires],isoVectorWidths_[ires]);
    }
    else {
      for(unsigned int ix=0;ix<isoVectorMasses_.size();++ix) {
	A1 += isoVectorKStarCoup_[ix]*Resonance::BreitWignerFW(q2,isoVectorMasses_[ix],isoVectorWidths_[ix]);
      }
    }
  }
  complex<InvEnergy3> amp(ZERO);
  ires = -1;
  if(ichan>=0) ires = ichan%2;
  if(imode==0) {
    complex<InvEnergy2> r1 = (ires<0||ires==0) ?
     Resonance::BreitWignerPWave(s1,mKStar0_,wKStar0_,momenta[0].mass(),momenta[2].mass())/sqr(mKStar0_) : InvEnergy2(); 
    complex<InvEnergy2> r2 = (ires<0||ires==1) ?
     Resonance::BreitWignerPWave(s2,mKStar0_,wKStar0_,momenta[1].mass(),momenta[2].mass())/sqr(mKStar0_) : InvEnergy2();
    amp = sqrt(1./6.)*(A0+A1)*(r1+r2);
  }
  else if(imode==1) {
    complex<InvEnergy2> r1 = (ires<0||ires==0) ?
     Resonance::BreitWignerPWave(s1,mKStarP_,wKStarP_,momenta[0].mass(),momenta[2].mass())/sqr(mKStarP_) : InvEnergy2(); 
    complex<InvEnergy2> r2 = (ires<0||ires==1) ?
     Resonance::BreitWignerPWave(s2,mKStarP_,wKStarP_,momenta[1].mass(),momenta[2].mass())/sqr(mKStarP_) : InvEnergy2();
    amp = sqrt(1./6.)*(A0-A1)*(r1+r2);
  }
  else {
    complex<InvEnergy2> r1 = (ires<0||ires==0) ?
     Resonance::BreitWignerPWave(s1,mKStarP_,wKStarP_,momenta[0].mass(),momenta[2].mass())/sqr(mKStarP_) : InvEnergy2(); 
    complex<InvEnergy2> r2 = (ires<0||ires==1) ?
     Resonance::BreitWignerPWave(s2,mKStar0_,wKStar0_,momenta[1].mass(),momenta[2].mass())/sqr(mKStar0_) : InvEnergy2();
    amp = sqrt(1./6.)*((A0+A1)*r1+(A0-A1)*r2);
  }
  amp *= 2.*gKStar_;
  // the current
  LorentzPolarizationVector vect = amp*Helicity::epsilon(momenta[0],momenta[1],momenta[2]);
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,scale*vect);
}
   
bool KKPiCurrent::accept(vector<int> id) {
  if(id.size()!=3) return false;
  int npip(0),npim(0),nkp(0),nkm(0),
    npi0(0),nks(0),nkl(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K_S0)    ++nks;
    else if(id[ix]==ParticleID::K_L0)    ++nkl;
  }
  if ( (npi0==1 && (( nks==1&&nkl==1 ) ||
		    ( nkp==1&&nkm==1 )) ) ||
       ( (nkl==1||nks==1) &&
	 ( (nkm==1&&npip==1) || (nkp==1&&npim==1) ) ) ) return true;
  return false;
}

// the decay mode
unsigned int KKPiCurrent::decayMode(vector<int> id) {
  assert(id.size()==3);
  int npip(0),npim(0),nkp(0),nkm(0),
    npi0(0),nks(0),nkl(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K_S0)    ++nks;
    else if(id[ix]==ParticleID::K_L0)    ++nkl;
  }
  if     ( nks==1&&nkl==1&&npi0==1 ) return 0;
  else if( nkp==1&&nkm==1&&npi0==1 ) return 1;
  else if( nks==1&&nkm==1&&npip==1 ) return 2;
  else if( nks==1&&nkp==1&&npim==1 ) return 3;
  else if( nkl==1&&nkm==1&&npip==1 ) return 4;
  else if( nkl==1&&nkp==1&&npim==1 ) return 5;
  assert(false);
}

// output the information for the database
void KKPiCurrent::dataBaseOutput(ofstream & output,bool header,
					bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::KKPiCurrent " 
  		    << name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<isoScalarMasses_.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoScalarMasses " << ix << " " << isoScalarMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<isoVectorMasses_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoVectorMasses " << ix << " " << isoVectorMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<isoScalarWidths_.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoScalarWidths " << ix << " " << isoScalarWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<isoVectorWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoVectorWidths " << ix << " " << isoVectorWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<isoScalarKStarAmp_.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoScalarKStarAmp " << ix << " " << isoScalarKStarAmp_[ix]*GeV << "\n";
  }
  for(unsigned int ix=0;ix<isoVectorKStarAmp_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoVectorKStarAmp " << ix << " " << isoVectorKStarAmp_[ix]*GeV << "\n";
  }
  for(unsigned int ix=0;ix<isoScalarKStarPhase_.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoScalarKStarPhase " << ix << " " << isoScalarKStarPhase_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<isoVectorKStarPhase_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoVectorKStarPhase " << ix << " " << isoVectorKStarPhase_[ix] << "\n";
  }
  output << "newdef " << name() << ":mKStarPlus " << mKStarP_/GeV  << "\n";
  output << "newdef " << name() << ":mKStar0 "    << mKStar0_/GeV  << "\n";
  output << "newdef " << name() << ":wKStarPlus " << wKStarP_/GeV  << "\n";
  output << "newdef " << name() << ":wKStar0 "    << wKStar0_/GeV  << "\n";
  output << "newdef " << name() << ":gKStar "     << gKStar_ << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
#line 1 "./EtaPhiCurrent.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPhiCurrent class.
//

#include "EtaPhiCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

EtaPhiCurrent::EtaPhiCurrent() {
  addDecayMode(3,-3);
  setInitialModes(3);
  // Masses for the resonances
  resMasses_ = {1.670*GeV,2.14*GeV};
  // widths for the resonances
  resWidths_ = {122*MeV,43.5*MeV};
  // amplitudes
  amp_   = {0.175/GeV,0.00409/GeV};
  // phases
  phase_ = {0.,2.19};
}

IBPtr EtaPhiCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr EtaPhiCurrent::fullclone() const {
  return new_ptr(*this);
}

void EtaPhiCurrent::doinit() {
  WeakCurrent::doinit();
  assert(phase_.size()==amp_.size());
  couplings_.clear();
  Complex ii(0.,1.);
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    couplings_.push_back(amp_[ix]*(cos(phase_[ix])+ii*sin(phase_[ix])));
  }
}

void EtaPhiCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(resMasses_,GeV) << ounit(resWidths_,GeV)
     << ounit(amp_,1./GeV) << phase_ << ounit(couplings_,1./GeV);
}

void EtaPhiCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(resMasses_,GeV) >> iunit(resWidths_,GeV)
     >> iunit(amp_,1./GeV) >> phase_ >> iunit(couplings_,1./GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EtaPhiCurrent,WeakCurrent>
describeHerwigEtaPhiCurrent("Herwig::EtaPhiCurrent",
			    "HwWeakCurrents.so");

void EtaPhiCurrent::Init() {

  static ClassDocumentation<EtaPhiCurrent> documentation
    ("The EtaPhiCurrent class implements a current based"
     " on the model of SND for eta + phi "
     "The current based on the model of \\cite{Achasov:2018ygm}"
     " for eta and phi was used.",
     "\\bibitem{Achasov:2018ygm}\n"
     "M.~N.~Achasov {\\it et al.},\n"
     "%``Measurement of the $e^+e^−\\to\\eta K^+K^−$ Cross Section by Means of the SND Detector,''\n"
     "Phys.\\ Atom.\\ Nucl.\\  {\\bf 81} (2018) no.2,  205\n"
     " [Yad.\\ Fiz.\\  {\\bf 81} (2018) no.2,  195].\n"
     "doi:10.1134/S1063778818020023\n"
     "%%CITATION = doi:10.1134/S1063778818020023;%%\n");

  static ParVector<EtaPhiCurrent,Energy> interfaceResonanceMasses
    ("ResonanceMasses",
     "The masses of the resonances for the form factor",
     &EtaPhiCurrent::resMasses_, GeV, 1, 1680*MeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<EtaPhiCurrent,Energy> interfaceResonanceWidths
    ("ResonanceWidths",
     "The widths of the resonances for the form factor",
     &EtaPhiCurrent::resWidths_, GeV, 1, 150*MeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<EtaPhiCurrent,InvEnergy> interfaceAmplitude
    ("Amplitude",
     "The amplitudes of the couplings",
     &EtaPhiCurrent::amp_, 0.00115/GeV, 1, 1./GeV, 0.0/GeV, 10/GeV,
     false, false, Interface::limited);

  static ParVector<EtaPhiCurrent,double> interfacePhase
    ("Phase",
     "The phases of the couplings in radians",
     &EtaPhiCurrent::phase_, 1, 0., 0.0, 2.*Constants::pi,
     false, false, Interface::limited);
}

// complete the construction of the decay mode for integration
bool EtaPhiCurrent::createMode(int icharge, tcPDPtr resonance,
			       FlavourInfo flavour,
			       unsigned int, PhaseSpaceModePtr mode,
			       unsigned int iloc,int ires,
			       PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(icharge!=0) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IZero) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    if(flavour.I3!=IsoSpin::I3Zero) return false;
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::ssbar)
    return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return false;
  // check that the mode is are kinematical allowed
  Energy min = getParticleData(ParticleID::eta)->mass()+
               getParticleData(ParticleID::phi)->massMin();
  if(min>upp) return false;
  // resonances for the intermediate channels
  tPDVector res = {getParticleData(100333)};
  // set up the integration channels;
  for(unsigned int ix=0;ix<res.size();++ix) {
    if(resonance && resonance!=res[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],
 		      ires+1,iloc+1,ires+1,iloc+2));
  }
  // reset the masses and widths of the resonances if needed
  for(unsigned int ix=0;ix<res.size();++ix) {
    mode->resetIntermediate(res[ix],resMasses_[ix],resWidths_[ix]);
  }
  return true;
}

// the particles produced by the current
tPDVector EtaPhiCurrent::particles(int icharge, unsigned int imode,int,int) {
  assert(icharge==0 && imode<=1);
  return {getParticleData(ParticleID::eta),getParticleData(ParticleID::phi)};
}

void EtaPhiCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    temp[ix] = HelicityFunctions::polarizationVector(-decay[1]->momentum(),
						     ix,Helicity::outgoing);
  }
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  VectorWaveFunction::constructSpinInfo(temp,decay[1],
					outgoing,true,true);
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
EtaPhiCurrent::current(tcPDPtr resonance,
		       FlavourInfo flavour,
		       const int, const int ichan, Energy & scale, 
		       const tPDVector & ,
		       const vector<Lorentz5Momentum> & momenta,
		       DecayIntegrator::MEOption) const {
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IZero) return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    if(flavour.I3!=IsoSpin::I3Zero) return vector<LorentzPolarizationVectorE>();
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::ssbar)
    return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return vector<LorentzPolarizationVectorE>();
  useMe();
  // polarization vectors of the photon
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // total momentum of the system
  Lorentz5Momentum q(momenta[0]+momenta[1]);
  // overall hadronic mass
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2(q.m2());
  unsigned int imin = 0;
  unsigned int imax = couplings_.size();
  if(ichan>0) {
    imin = ichan;
    imax = imin+1;
  }
  if(resonance) {
    switch(abs(resonance->id())) {
    case 100333:
      imin=0;
      break;
    default:
      assert(false);
    }
    imax=imin+1;
  }
  // compute the form factor
  complex<InvEnergy> formFactor(ZERO);
  // loop over the resonances
  for(unsigned int ix=imin;ix<imax;++ix) {
    Energy2 mR2(sqr(resMasses_[ix]));
    // compute the width
    Energy width = resWidths_[ix];
    formFactor += couplings_[ix]*mR2/(mR2-q2-Complex(0.,1.)*q.mass()*width);
  }
  // calculate the current
  vector<LorentzPolarizationVectorE> ret(3);
  for(unsigned int ix=0;ix<3;++ix) {
    ret[ix] += formFactor*Helicity::epsilon(q,temp[ix],momenta[1]);
  }
  return ret;
}

bool EtaPhiCurrent::accept(vector<int> id) {
  if(id.size()!=2) return false;
  unsigned int neta(0),nphi(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::eta)   ++neta;
    else if(id[ix]==ParticleID::phi) ++nphi;
  }
  return nphi == 1 && neta==1;
}

unsigned int EtaPhiCurrent::decayMode(vector<int>) {
  return 0;
}

// output the information for the database
void EtaPhiCurrent::dataBaseOutput(ofstream & output,bool header,
				   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::EtaPhiCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<resMasses_.size();++ix) {
    if(ix<1) output << "newdef " << name() << ":ResonanceMasses " << ix 
		    << " " << resMasses_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":ResonanceMasses " << ix 
		    << " " << resMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<resWidths_.size();++ix) {
    if(ix<1) output << "newdef " << name() << ":ResonanceWidths " << ix 
		    << " " << resWidths_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":ResonanceWidths " << ix 
		    << " " << resWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    if(ix<1) output << "newdef " << name() << ":Amplitude " << ix 
		    << " " << amp_[ix]*GeV << "\n";
    else     output << "insert " << name() << ":Amplitude " << ix 
		    << " " << amp_[ix]*GeV << "\n";
  }
  for(unsigned int ix=0;ix<phase_.size();++ix) {
    if(ix<1) output << "newdef " << name() << ":Phase " << ix 
		    << " " << phase_[ix] << "\n";
    else     output << "insert " << name() << ":Phase " << ix 
		    << " " << phase_[ix] << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./EtaOmegaCurrent.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaOmegaCurrent class.
//

#include "EtaOmegaCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

EtaOmegaCurrent::EtaOmegaCurrent() {
  addDecayMode(3,-3);
  setInitialModes(3);
  // Masses for the resonances
  resMasses_ = {1.425*GeV,1.67*GeV};
  // widths for the resonances
  resWidths_ = {215*MeV  , 113*MeV};
  // amplitudes
  amp_   = {0.0862/GeV,0.0648/GeV};
  // phases
  phase_ = {0.,180.};
}

IBPtr EtaOmegaCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr EtaOmegaCurrent::fullclone() const {
  return new_ptr(*this);
}

void EtaOmegaCurrent::doinit() {
  WeakCurrent::doinit();
  assert(phase_.size()==amp_.size());
  couplings_.clear();
  Complex ii(0.,1.);
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    double phi = phase_[ix]/180.*Constants::pi;
    couplings_.push_back(amp_[ix]*(cos(phi)+ii*sin(phi)));
  }
}

void EtaOmegaCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(resMasses_,GeV) << ounit(resWidths_,GeV)
     << ounit(amp_,1./GeV) << phase_ << ounit(couplings_,1./GeV);
}

void EtaOmegaCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(resMasses_,GeV) >> iunit(resWidths_,GeV)
     >> iunit(amp_,1./GeV) >> phase_ >> iunit(couplings_,1./GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EtaOmegaCurrent,WeakCurrent>
describeHerwigEtaOmegaCurrent("Herwig::EtaOmegaCurrent",
			    "HwWeakCurrents.so");

void EtaOmegaCurrent::Init() {

  static ClassDocumentation<EtaOmegaCurrent> documentation
    ("The EtaOmegaCurrent class implements a current based"
     " on the model of SND for eta + omega "
     "The current based on the model of \\cite{Achasov:2016qvd}"
     " for eta and omega was used.",
     "\\bibitem{Achasov:2016qvd}\n"
     "M.~N.~Achasov {\\it et al.},\n"
     "%``Measurement of the $e^+e^- \\to \\omega\\eta$ cross section below $\\sqrt{s}=2$ GeV,''\n"
     "Phys.\\ Rev.\\ D {\\bf 94} (2016) no.9,  092002\n"
     "doi:10.1103/PhysRevD.94.092002\n"
     "[arXiv:1607.00371 [hep-ex]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.94.092002;%%\n"
     "%18 citations counted in INSPIRE as of 12 Oct 2018\n");

  static ParVector<EtaOmegaCurrent,Energy> interfaceResonanceMasses
    ("ResonanceMasses",
     "The masses of the resonances for the form factor",
     &EtaOmegaCurrent::resMasses_, GeV, 1, 1680*MeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<EtaOmegaCurrent,Energy> interfaceResonanceWidths
    ("ResonanceWidths",
     "The widths of the resonances for the form factor",
     &EtaOmegaCurrent::resWidths_, GeV, 1, 150*MeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<EtaOmegaCurrent,InvEnergy> interfaceAmplitude
    ("Amplitude",
     "The amplitudes of the couplings",
     &EtaOmegaCurrent::amp_, 1./GeV, 2, 0.0648/GeV, 0.0/GeV, 10/GeV,
     false, false, Interface::limited);

  static ParVector<EtaOmegaCurrent,double> interfacePhase
    ("Phase",
     "The phases of the couplings in degrees",
     &EtaOmegaCurrent::phase_, 1, 0., 0.0, 360.0,
     false, false, Interface::limited);
}

// complete the construction of the decay mode for integration
bool EtaOmegaCurrent::createMode(int icharge, tcPDPtr resonance,
			       FlavourInfo flavour,
			       unsigned int, PhaseSpaceModePtr mode,
			       unsigned int iloc,int ires,
			       PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(icharge!=0) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IZero) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    if(flavour.I3!=IsoSpin::I3Zero) return false;
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero ) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       ) return false;
  // check that the mode is are kinematical allowed
  Energy min = getParticleData(ParticleID::eta)->mass()+
               getParticleData(ParticleID::omega)->massMin();
  if(min>upp) return false;
  // resonances for the intermediate channels
  tPDVector res = {getParticleData(100223),getParticleData(30223)};
  // set up the integration channels;
  for(unsigned int ix=0;ix<res.size();++ix) {
    if(resonance && resonance!=res[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],
 		      ires+1,iloc+1,ires+1,iloc+2));
  }
  // reset the masses and widths of the resonances if needed
  for(unsigned int ix=0;ix<res.size();++ix) {
    mode->resetIntermediate(res[ix],resMasses_[ix],resWidths_[ix]);
  }
  return true;
}

// the particles produced by the current
tPDVector EtaOmegaCurrent::particles(int icharge, unsigned int imode,int,int) {
  assert(icharge==0 && imode<=1);
  return {getParticleData(ParticleID::eta),getParticleData(ParticleID::omega)};
}

void EtaOmegaCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    temp[ix] = HelicityFunctions::polarizationVector(-decay[1]->momentum(),
						     ix,Helicity::outgoing);
  }
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  VectorWaveFunction::constructSpinInfo(temp,decay[1],
					outgoing,true,true);
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
EtaOmegaCurrent::current(tcPDPtr resonance,
		       FlavourInfo flavour,
		       const int, const int ichan, Energy & scale, 
		       const tPDVector & ,
		       const vector<Lorentz5Momentum> & momenta,
		       DecayIntegrator::MEOption) const {
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IZero) return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    if(flavour.I3!=IsoSpin::I3Zero) return vector<LorentzPolarizationVectorE>();
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return vector<LorentzPolarizationVectorE>();
  useMe();
  // polarization vectors of the photon
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // total momentum of the system
  Lorentz5Momentum q(momenta[0]+momenta[1]);
  // overall hadronic mass
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2(q.m2());
  unsigned int imin = 0;
  unsigned int imax = couplings_.size();
  if(ichan>0) {
    imin = ichan;
    imax = imin+1;
  }
  if(resonance) {
    switch(abs(resonance->id())) {
    case 100223:
      imin=0;
      break;
    default:
      assert(false);
    }
    imax=imin+1;
  }
  // compute the form factor
  complex<InvEnergy> formFactor(ZERO);
  // loop over the resonances
  for(unsigned int ix=imin;ix<imax;++ix) {
    Energy2 mR2(sqr(resMasses_[ix]));
    // compute the width
    Energy width = resWidths_[ix];
    formFactor += couplings_[ix]*mR2/(mR2-q2-Complex(0.,1.)*q.mass()*width);
  }
  // calculate the current
  vector<LorentzPolarizationVectorE> ret(3);
  for(unsigned int ix=0;ix<3;++ix) {
    ret[ix] += formFactor*Helicity::epsilon(q,temp[ix],momenta[1]);
  }
  return ret;
}

bool EtaOmegaCurrent::accept(vector<int> id) {
  if(id.size()!=2) return false;
  unsigned int neta(0),nomega(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::eta)   ++neta;
    else if(id[ix]==ParticleID::omega) ++nomega;
  }
  return nomega == 1 && neta==1;
}

unsigned int EtaOmegaCurrent::decayMode(vector<int>) {
  return 0;
}

// output the information for the database
void EtaOmegaCurrent::dataBaseOutput(ofstream & output,bool header,
				   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::EtaOmegaCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<resMasses_.size();++ix) {
    if(ix<1) output << "newdef " << name() << ":ResonanceMasses " << ix 
		    << " " << resMasses_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":ResonanceMasses " << ix 
		    << " " << resMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<resWidths_.size();++ix) {
    if(ix<1) output << "newdef " << name() << ":ResonanceWidths " << ix 
		    << " " << resWidths_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":ResonanceWidths " << ix 
		    << " " << resWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    if(ix<1) output << "newdef " << name() << ":Amplitude " << ix 
		    << " " << amp_[ix]*GeV << "\n";
    else     output << "insert " << name() << ":Amplitude " << ix 
		    << " " << amp_[ix]*GeV << "\n";
  }
  for(unsigned int ix=0;ix<phase_.size();++ix) {
    if(ix<1) output << "newdef " << name() << ":Phase " << ix 
		    << " " << phase_[ix] << "\n";
    else     output << "insert " << name() << ":Phase " << ix 
		    << " " << phase_[ix] << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./OmegaPiPiCurrent.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OmegaPiPiCurrent class.
//

#include "OmegaPiPiCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

OmegaPiPiCurrent::OmegaPiPiCurrent() {
  mRes_ = 1.69*GeV;
  wRes_ = 0.285*GeV;
  gRes_ = 1.63*GeV;
  
  mSigma_ = 0.6*GeV;
  wSigma_ = 1.0*GeV;
  mf0_ = 0.98*GeV;
  gSigma_ = 1.0 *GeV2;
  gf0_    = 0.883*GeV2;
  gPiPi_ = 0.165*GeV2;
  gKK_   = 0.695*GeV2;
  addDecayMode(1,-1);
  addDecayMode(1,-1);
  setInitialModes(2);
}

IBPtr OmegaPiPiCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr OmegaPiPiCurrent::fullclone() const {
  return new_ptr(*this);
}

void OmegaPiPiCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(mRes_,GeV) << ounit(wRes_,GeV) << ounit(gRes_,GeV)
     << ounit(mSigma_,GeV) << ounit(wSigma_,GeV) << ounit(mf0_,GeV)
     << ounit(gPiPi_,GeV2) << ounit(gKK_,GeV2) << ounit(gSigma_,GeV2) << ounit(gf0_,GeV2);
}

void OmegaPiPiCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(mRes_,GeV) >> iunit(wRes_,GeV) >> iunit(gRes_,GeV)
     >> iunit(mSigma_,GeV) >> iunit(wSigma_,GeV) >> iunit(mf0_,GeV)
     >> iunit(gPiPi_,GeV2) >> iunit(gKK_,GeV2) >> iunit(gSigma_,GeV2) >> iunit(gf0_,GeV2);
}

void OmegaPiPiCurrent::doinit() {
  WeakCurrent::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<OmegaPiPiCurrent,WeakCurrent>
describeHerwigOmegaPiPiCurrent("Herwig::OmegaPiPiCurrent", "HwWeakCurrents.so");

void OmegaPiPiCurrent::Init() {

  static ClassDocumentation<OmegaPiPiCurrent> documentation
    ("The OmegaPiPiCurrent class provides the current for I=0 omega pi pi");
  
  static Parameter<OmegaPiPiCurrent,Energy> interfacemRes
    ("mRes",
     "The mass of the s-channel resonance",
     &OmegaPiPiCurrent::mRes_, GeV, 1.62*GeV, 0.*GeV, 10.*GeV,
     false, false, Interface::limited);
  
  static Parameter<OmegaPiPiCurrent,Energy> interfacewRes
    ("wRes",
     "The width of the s-channel resonance",
     &OmegaPiPiCurrent::wRes_, GeV, 0.288*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<OmegaPiPiCurrent,Energy> interfacegRes
    ("gRes",
     "The coupling of the s-channel resonance",
     &OmegaPiPiCurrent::gRes_, GeV, 2.83*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<OmegaPiPiCurrent,Energy> interfacemSigma
    ("mSigma",
     "The mass of the Sigma",
     &OmegaPiPiCurrent::mSigma_, GeV, 0.6*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<OmegaPiPiCurrent,Energy> interfacewSigma
    ("wSigma",
     "The width of the Sigma",
     &OmegaPiPiCurrent::wSigma_, GeV, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<OmegaPiPiCurrent,Energy2> interfacegSigma
    ("gSigma",
     "The coupling of the Sigma resonance",
     &OmegaPiPiCurrent::gSigma_, GeV2, 1.0*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);
  
  static Parameter<OmegaPiPiCurrent,Energy> interfacemf0
    ("mf0",
     "The mass of the f_0(980)",
     &OmegaPiPiCurrent::mf0_, GeV, 0.98*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<OmegaPiPiCurrent,Energy2> interfacegf0
    ("gf0",
     "The coupling of the f0(980) resonance",
     &OmegaPiPiCurrent::gf0_, GeV2, 1.0*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);
  
  static Parameter<OmegaPiPiCurrent,Energy2> interfacegPiPi
    ("gPiPi",
     "The coupling of the f_0(980) to pipi",
     &OmegaPiPiCurrent::gPiPi_, GeV2, 0.165*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);
  
  static Parameter<OmegaPiPiCurrent,Energy2> interfacegKK
    ("gKK",
     "The coupling of the f_0(980) to KK",
     &OmegaPiPiCurrent::gKK_, GeV2, 0.695*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);
}



// complete the construction of the decay mode for integration
bool OmegaPiPiCurrent::createMode(int icharge, tcPDPtr resonance,
			       FlavourInfo flavour,
			       unsigned int, PhaseSpaceModePtr mode,
			       unsigned int iloc,int ires,
			       PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(icharge!=0) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IZero) return false;
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown && flavour.I3!=IsoSpin::I3Zero) return false;
  // and other flavour
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
  // check that the mode is are kinematical allowed
  Energy min = getParticleData(ParticleID::omega)->massMin()+
    2.*getParticleData(ParticleID::pi0)->mass();
  if(min>upp) return false;
  // resonances for the intermediate channels
  tPDVector res = {getParticleData(30223)};
  tPDVector res2 = {getParticleData(9000221),getParticleData(9010221)};
  // set up the integration channels;
  for(unsigned int ix=0;ix<res.size();++ix) {
    if(resonance && resonance!=res[ix]) continue;
    for(unsigned int iy=0;iy<res2.size();++iy) {
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],
   		      ires+1,iloc+1,ires+1,res2[iy],ires+2,iloc+2,ires+2,iloc+3));
    }
  }
  return true;
}

// the particles produced by the current
tPDVector OmegaPiPiCurrent::particles(int icharge, unsigned int imode,int,int) {
  assert(icharge==0 && imode<=1);
  if(imode==0) 
    return {getParticleData(ParticleID::omega),
	    getParticleData(ParticleID::piplus),
	    getParticleData(ParticleID::piminus)};
  else if(imode==1) 
    return {getParticleData(ParticleID::omega),
	getParticleData(ParticleID::pi0),
	getParticleData(ParticleID::pi0)};
  else
    assert(false);
}

void OmegaPiPiCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    temp[ix] = HelicityFunctions::polarizationVector(-decay[0]->momentum(),
  						     ix,Helicity::outgoing);
  }
  VectorWaveFunction::constructSpinInfo(temp,decay[0],
  					outgoing,true,true);
  for(unsigned int ix=1;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
OmegaPiPiCurrent::current(tcPDPtr resonance,
			  FlavourInfo flavour,
			  const int, const int ichan, Energy & scale, 
			  const tPDVector & ,
			  const vector<Lorentz5Momentum> & momenta,
			  DecayIntegrator::MEOption) const {
  // no isospin/flavour here
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IZero) return vector<LorentzPolarizationVectorE>();
  if(flavour.I3!=IsoSpin::I3Unknown && flavour.I3!=IsoSpin::I3Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(resonance and resonance->id()!=30223) return vector<LorentzPolarizationVectorE>();
  useMe();
  // polarization vectors of the omega
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[0],ix,Helicity::outgoing);
  }
  // total momentum of the system
  Lorentz5Momentum q(momenta[0]+momenta[1]+momenta[2]);
  // overall hadronic mass
  q.rescaleMass();
  scale=q.mass();
  Complex ii(0.,1.);
  Energy2 q2(q.m2());
  // resonance factor for s channel resonance
  Energy2 mR2=sqr(mRes_);
  complex<Energy> pre= mR2*gRes_/(q2-mR2 - ii*scale*wRes_);
  
  //cerr << "pre factor " << scale/GeV << " " << q2/GeV2 << " " << pre/GeV << "\n";
  // for(auto p : momenta) cerr << p/GeV << " " << p.m()/GeV << "\n";

  
  
  // virtual mass energy for intermediate f0 channel
  Energy2 s1 = (momenta[1]+momenta[2]).m2();
  Energy sqrs1 = sqrt(s1);
  
  //sigma meson
  Energy2 mSigma2 = sqr(mSigma_);
  Complex Sigma_form = gSigma_/(mSigma2 -s1 - ii*sqrs1*wSigma_);

  // compute the form factor
  Energy2 mf02 = sqr(mf0_);
  Energy   mPi = getParticleData(211)->mass();
  Energy2 mPi2 = sqr(mPi);

  complex<Energy2> mGamma = gPiPi_*sqrt(max(0.,1.-4.*mPi2/s1));
  // cerr << "testing pi " << mGamma/GeV2 << " " << gPiPi_/GeV2 << " " << sqrt(max(0.,1.-4.*mPi2/s1))<< "\n";
  Energy2 mKp2 = sqr(getParticleData(321)->mass());
  double val = 1.-4.*mKp2/s1;
  if(val>=0.) mGamma += 0.5*gKK_*   sqrt( val);
  else        mGamma += 0.5*gKK_*ii*sqrt(-val);
  Energy2 mK02 = sqr(getParticleData(311)->mass());
  val = 1.-4.*mK02/s1;
  if(val>=0.) mGamma += 0.5*gKK_*   sqrt( val);
  else        mGamma += 0.5*gKK_*ii*sqrt(-val);
  
  Complex f0_form = gf0_/(mf02-s1-Complex(0.,1.)*mGamma);
  // cerr << "f0 pieces " << mGamma/GeV2 << "\n";
  // cerr << "testing form factor " << s1/GeV2 << " " << f0_form << " " << Sigma_form << "\n";
  complex<Energy> formFactor(ZERO);
  if(ichan<0) 
    formFactor = (Sigma_form+f0_form)*pre;
  else if(ichan==0)
    formFactor = Sigma_form*pre;
  else if(ichan==1)
    formFactor = f0_form*pre;
  // calculate the current
  vector<LorentzPolarizationVectorE> ret(3);
  for(unsigned int ix=0;ix<3;++ix) {
    ret[ix] = formFactor*temp[ix];
  }
  return ret;
}

bool OmegaPiPiCurrent::accept(vector<int> id) {
  if(id.size()!=3) return false;
   unsigned int nomega(0),npip(0),npim(0),npi0(0);
   for(unsigned int ix=0;ix<id.size();++ix) {
     if(abs(id[ix])==ParticleID::piminus) ++npim;
     if(abs(id[ix])==ParticleID::pi0    ) ++npi0;
     if(abs(id[ix])==ParticleID::piplus ) ++npip;
     else if(id[ix]==ParticleID::omega)  ++nomega;
   }
   return nomega == 1 && (npi0==2 || (npip==1&&npim==1));
}

unsigned int OmegaPiPiCurrent::decayMode(vector<int> id) {
   unsigned int npi0(0);
   for(unsigned int ix=0;ix<id.size();++ix) {
     if(abs(id[ix])==ParticleID::pi0    ) ++npi0;
   }
   return npi0==2 ? 1 : 0;
}

// output the information for the database
void OmegaPiPiCurrent::dataBaseOutput(ofstream & output,bool header,
				   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::OmegaPiPiCurrent " << name() 
  		    << " HwWeakCurrents.so\n";
  output << "newdef " << name() << ":mRes "   << " " << mRes_/GeV   << "\n";
  output << "newdef " << name() << ":wRes "   << " " << wRes_/GeV   << "\n";
  output << "newdef " << name() << ":mSigma " << " " << mSigma_/GeV << "\n";
  output << "newdef " << name() << ":wSigma " << " " << wSigma_/GeV << "\n";
  output << "newdef " << name() << ":mf0 "    << " " << mf0_/GeV    << "\n";
  output << "newdef " << name() << ":gRes "   << " " << gRes_/GeV   << "\n";
  output << "newdef " << name() << ":gSigma " << " " << gSigma_/GeV2<< "\n";
  output << "newdef " << name() << ":gf0 "    << " " << gf0_/GeV2   << "\n";
  output << "newdef " << name() << ":gPiPi "  << " " << gPiPi_/GeV2 << "\n";
  output << "newdef " << name() << ":gKK "    << " " << gKK_/GeV2   << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
   		    << fullName() << "\";" << endl;
}
