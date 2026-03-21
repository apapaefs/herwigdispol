// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PolarizedPartonExtractor class.
//

#include "PolarizedPartonExtractor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;
IBPtr PolarizedPartonExtractor::clone() const {
  return new_ptr(*this);
}

IBPtr PolarizedPartonExtractor::fullclone() const {
  return new_ptr(*this);
}

void PolarizedPartonExtractor::persistentOutput(PersistentOStream & os) const {
  os << theFirstLongPol << theSecondLongPol << theFirstLongDiffPDF << theSecondLongDiffPDF
     << theFirstTransXPol << theSecondTransXPol << theFirstTransYPol << theSecondTransYPol
     << theFirstTransDiffPDF << theSecondTransDiffPDF
     << isPolarizedFirst << isPolarizedSecond;
}

void PolarizedPartonExtractor::persistentInput(PersistentIStream & is, int) {
  is >> theFirstLongPol >> theSecondLongPol >> theFirstLongDiffPDF >> theSecondLongDiffPDF
     >> theFirstTransXPol >> theSecondTransXPol >> theFirstTransYPol >> theSecondTransYPol
     >> theFirstTransDiffPDF >> theSecondTransDiffPDF
     >> isPolarizedFirst >> isPolarizedSecond;
}

void PolarizedPartonExtractor::doinit() {
  PartonExtractor::doinit();
  isPolarizedFirst  = theFirstLongPol  != 0. || theFirstTransXPol  != 0. || theFirstTransYPol  != 0.;
  isPolarizedSecond = theSecondLongPol != 0. || theSecondTransXPol != 0. || theSecondTransYPol != 0.;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PolarizedPartonExtractor,PartonExtractor>
describeThePEGPolarizedPartonExtractor("ThePEG::PolarizedPartonExtractor", "libThePEG.so");

void PolarizedPartonExtractor::Init() {

  static ClassDocumentation<PolarizedPartonExtractor> documentation
    ("The PolarizedPartonExtractor class inherits from the PartonExtractor"
     " and provides the functionality for polorized beams");

  static Parameter<PolarizedPartonExtractor,double> interfaceFirstLongPolarization
    ("FirstLongitudinalPolarization",
     "The longitudinal polarization of the first beam",
     &PolarizedPartonExtractor::theFirstLongPol, 0.0, -1., 1.,
     false, false, Interface::limited);

  static Parameter<PolarizedPartonExtractor,double> interfaceSecondLongPolarization
    ("SecondLongitudinalPolarization",
     "The longitudinal polarization of the second beam",
     &PolarizedPartonExtractor::theSecondLongPol, 0.0, -1., 1.,
     false, false, Interface::limited);

  static Reference<PolarizedPartonExtractor,PDFBase> interfaceFirstLongDifferencePDF
    ("FirstLongitudinalDifferencePDF",
     "The PDF providing the differences in longitudinal helicities for the first particle",
     &PolarizedPartonExtractor::theFirstLongDiffPDF, false, false, true, true, false);

  static Reference<PolarizedPartonExtractor,PDFBase> interfaceSecondLongDifferencePDF
    ("SecondLongitudinalDifferencePDF",
     "The PDF providing the differences in longitudinal helicities for the second particle",
     &PolarizedPartonExtractor::theSecondLongDiffPDF, false, false, true, true, false);

  static Parameter<PolarizedPartonExtractor,double> interfaceFirstTransXPolarization
    ("FirstTransverseXPolarization",
     "The transverse x polarization of the first beam",
     &PolarizedPartonExtractor::theFirstTransXPol, 0.0, -1., 1.,
     false, false, Interface::limited);

  static Parameter<PolarizedPartonExtractor,double> interfaceSecondTransXPolarization
    ("SecondTransverseXPolarization",
     "The transverse x polarization of the second beam",
     &PolarizedPartonExtractor::theSecondTransXPol, 0.0, -1., 1.,
     false, false, Interface::limited);

  static Parameter<PolarizedPartonExtractor,double> interfaceFirstTransYPolarization
    ("FirstTransverseYPolarization",
     "The transverse y polarization of the first beam",
     &PolarizedPartonExtractor::theFirstTransYPol, 0.0, -1., 1.,
     false, false, Interface::limited);

  static Parameter<PolarizedPartonExtractor,double> interfaceSecondTransYPolarization
    ("SecondTransverseYPolarization",
     "The transverse y polarization of the second beam",
     &PolarizedPartonExtractor::theSecondTransYPol, 0.0, -1., 1.,
     false, false, Interface::limited);

  static Reference<PolarizedPartonExtractor,PDFBase> interfaceFirstTransDifferencePDF
    ("FirstTransverseDifferencePDF",
     "The PDF providing the differences in transverse helicities for the first particle",
     &PolarizedPartonExtractor::theFirstTransDiffPDF, false, false, true, true, false);

  static Reference<PolarizedPartonExtractor,PDFBase> interfaceSecondTransDifferencePDF
    ("SecondTransverseDifferencePDF",
     "The PDF providing the differences in transverse helicities for the second particle",
     &PolarizedPartonExtractor::theSecondTransDiffPDF, false, false, true, true, false);

}

pair<RhoDMatrix,RhoDMatrix> PolarizedPartonExtractor::getRhoMatrices(const PBIPair & pbins) const {
  // spin density matrices for output
  pair<RhoDMatrix,RhoDMatrix> output = make_pair(RhoDMatrix(pbins.first ->partonData()->iSpin()), 
                                                 RhoDMatrix(pbins.second->partonData()->iSpin()));
  if(pbins.first ->partonData()->id()==21 || pbins.first ->partonData()->id()==22) {
    output.first(0,0) = 0.5;
    output.first(1,1) = 0.0;
    output.first(2,2) = 0.5;
  }
  if(pbins.second ->partonData()->id()==21 || pbins.second ->partonData()->id()==22) {
    output.second(0,0) = 0.5;
    output.second(1,1) = 0.0;
    output.second(2,2) = 0.5;
  }
  // first parton
  if(isPolarizedFirst) {
    unsigned int imax = output.first.iSpin()-1;
    double  pL(theFirstLongPol);
    Complex pT(theFirstTransXPol,theFirstTransYPol);
    // case with pdf
    if (pbins.first ->particleData()  &&  pbins.first ->particleData() != pbins.first ->partonData() ) {
      double sum  = pbins.first->pdf()    ->xfl(pbins.first->particleData(), pbins.first->partonData(), pbins.first->scale(),pbins.first->li(), pbins.first->incoming()->scale());
      // longitudinal polarization
      if (  theFirstLongPol != 0. && theFirstLongDiffPDF) {
        double diff = theFirstLongDiffPDF->xfl(pbins.first->particleData(), pbins.first->partonData(), pbins.first->scale(),pbins.first->li(), pbins.first->incoming()->scale());
        pL *= diff/sum;
      }
      // transverse polarization
      if( (theFirstTransXPol !=0. || theFirstTransYPol !=0.) && theFirstTransDiffPDF) {
        double diff = theFirstTransDiffPDF->xfl(pbins.first->particleData(), pbins.first->partonData(), pbins.first->scale(),pbins.first->li(), pbins.first->incoming()->scale());
        pT *= diff/sum;
      }
    }
    // setup the spin density matrix
    output.first(0   ,0   ) = 0.5*(1.-pL);
    output.first(imax,imax) = 0.5*(1.+pL);
    output.first(0   ,imax) = 0.5*pT;
    output.first(imax,0   ) = 0.5*conj(pT);
  }
  // second parton
  if(isPolarizedSecond) {
    unsigned int imax = output.second.iSpin()-1;
    double pL = theSecondLongPol;
    Complex pT(theSecondTransXPol,theSecondTransYPol);
    // case with pdf
    if (pbins.second ->particleData()  &&  pbins.second ->particleData() != pbins.second ->partonData() ) {
      double sum  = pbins.second->pdf()    ->xfl(pbins.second->particleData(), pbins.second->partonData(), pbins.second->scale(),pbins.second->li(), pbins.second->incoming()->scale());
      // longitudinal polarization
      if (  theSecondLongPol != 0. && theSecondLongDiffPDF) {
        double diff = theSecondLongDiffPDF->xfl(pbins.second->particleData(), pbins.second->partonData(), pbins.second->scale(),pbins.second->li(), pbins.second->incoming()->scale());
        pL *= diff/sum;
      }
      // transverse polarization
      if( (theSecondTransXPol !=0. || theSecondTransYPol !=0.) && theSecondTransDiffPDF) {
        double diff = theSecondTransDiffPDF->xfl(pbins.second->particleData(), pbins.second->partonData(), pbins.second->scale(),pbins.second->li(), pbins.second->incoming()->scale());
        pT *= diff/sum;
      }
    }
    // setup the spin density matrix
    output.second(0   ,0   ) = 0.5*(1.-pL);
    output.second(imax,imax) = 0.5*(1.+pL);
    output.second(0   ,imax) = 0.5*pT;
    output.second(imax,0   ) = 0.5*conj(pT);
  }
  // return the result
  return output;
}
