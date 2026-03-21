// -*- C++ -*-
#ifndef ThePEG_PolarizedPartonExtractor_H
#define ThePEG_PolarizedPartonExtractor_H
//
// This is the declaration of the PolarizedPartonExtractor class.
//

#include "PartonExtractor.h"
#include "ThePEG/Vectors/ThreeVector.h"

namespace ThePEG {

/**
 * Here is the documentation of the PolarizedPartonExtractor class.
 *
 * @see \ref PolarizedPartonExtractorInterfaces "The interfaces"
 * defined for PolarizedPartonExtractor.
 */
class PolarizedPartonExtractor: public PartonExtractor {

public:

  /**
   *    Spin density matrices for the incoming particles
   */
  pair<RhoDMatrix,RhoDMatrix> getRhoMatrices(const PBIPair & pbins) const;

public:

  /**
   *  Whether or not the beam is polarized
   */
  pair<bool,bool> isPolarized() const {
    return make_pair(isPolarizedFirst,isPolarizedSecond);
  }
  
  /**
   *  Polarization of the beams
   */
  ThreeVector<double> polarization(bool first=true) const {
    return ThreeVector<double>(first ? theFirstTransXPol :  theSecondTransXPol,
                               first ? theFirstTransYPol :  theSecondTransYPol,
                               first ? theFirstLongPol   :  theSecondLongPol  );
  }

  /**
   *  Longitudinal difference PDF
   */
  pair<PDFPtr,PDFPtr> longitudinalDifferencePDF() const {
    return make_pair(theFirstLongDiffPDF,theSecondLongDiffPDF);
  }

  /**
   *  Transverse difference PDF
   */
  pair<PDFPtr,PDFPtr> transverseDifferencePDF() const {
    return make_pair(theFirstTransDiffPDF,theSecondTransDiffPDF);
  }
  
public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PolarizedPartonExtractor & operator=(const PolarizedPartonExtractor &) = delete;

private:

  /**
   *  Longitudinal polarization of the first beam
   */
  double theFirstLongPol=0.;

  /**
   *  Longitudinal polarization of the second beam
   */
  double theSecondLongPol=0.;

  /**
   *  PDFBase object for first longitusinal difference PDF
   */
  PDFPtr theFirstLongDiffPDF;

  /**
   *  PDFBase object for second longitusinal difference PDF
   */
  PDFPtr theSecondLongDiffPDF;

  /**
   *  Transverse x polarization of the first beam
   */
  double theFirstTransXPol=0.;

  /**
   *  Transverse x polarization of the second beam
   */
  double theSecondTransXPol=0.;

  /**
   *  Transverse y polarization of the first beam
   */
  double theFirstTransYPol=0.;

  /**
   *  Transverse y polarization of the second beam
   */
  double theSecondTransYPol=0.;

  /**
   *  PDFBase object for first transverse difference PDF
   */
  PDFPtr theFirstTransDiffPDF;

  /**
   *  PDFBase object for second transverse difference PDF
   */
  PDFPtr theSecondTransDiffPDF;

  /**
   *   Whether need polarization for first particle
   */
  bool isPolarizedFirst;

  /**
   *   Whether need polarization for second particle
   */
  bool isPolarizedSecond;

};

}

#endif /* ThePEG_PolarizedPartonExtractor_H */
