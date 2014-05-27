/***************************************************************************
 *   Copyright (C) 2008 by Antonio Jesús Adsuar Gómez                      *
 *   adsuar@lsi.upc.edu                                                    *
 *                                                                         *
 *   This file is part of CORTEX                                           *
 *   CORTEX is free software; you can redistribute it and/or modify        *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   CORTEX is distributed in the hope that it will be useful,             *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
 
/****************************************************************************
 ** CORTEX
 ** Tool designed by:
 **    Antonio Jesús Adsuar Gómez
 **    adsuar@lsi.upc.edu
 ** Sevilla 2008
 **
 **
 ** File: fsevaluation.h
 ** Description: Definition of the Evaluation Class, that estimates the
 **              goodness of a feature selection
 ***************************************************************************/
 
#ifndef _EVALUATION_H_
#define _EVALUATION_H_

#include "../fsvector/fsvector.h"
#include "../fsmatrix/fsmatrix.h"
#include "../fsvotes/fsvotes.h"
#include "../fsdata/fsdata.h"
#include "../fsdefinitions/fstypes.h"
#include "../fsdistance/fsdistance.h"

using namespace std;

class fsEvaluation
{
   private:
      fsReal (fsEvaluation::*eval)(fsVector *, fsUInt);
      fsReal NearestNeighbours(fsVector*, fsUInt);
      fsReal NaiveBayes(fsVector*, fsUInt);
      fsReal Scorer(fsVector *, fsUInt);
      
      fsVoid initializeNaiveBayes();
      
      fsEvaluationAlgorithm fsEA;                 // Evaluation Algorithm
      fsDistanceAlgorithm   fsDA;                 // Distance Algorithm
      fsUInt kFoldCrossValidation;                // K-Fold Cross-Validation
      fsUInt kFoldSteps;                          // Size of K-Fold Blocks
      fsBool randomized;                          // Data Splitted Randomly
      fsUInt firstBlock;                          // First block of splitted
                                                  // data
      fsDistance *fsDist;                         // Distance Algorithm

      fsData *study;
      fsData *studyTest;

      vector< vector< vector<fsReal> > > matrixNaiveBayes;
      vector< vector< vector<fsReal> > > matrixNaiveBayesTest;
      
      fsVotes items;
   public:
      fsEvaluation(fsEvaluationAlgorithm fsea=NN,
                   fsDistanceAlgorithm fsda=SQUAREEUCLIDEANDISTANCE,
                   fsUInt kFCV=2,
                   fsBool r=false,
                   fsChar *fsFileData=NULL);
      ~fsEvaluation();
      
      fsReal J(fsVector *);
      fsUInt numberOfFeatures();
      fsUInt numberOfObservations();
};

#endif // _EVALUATION_H_
