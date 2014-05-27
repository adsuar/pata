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
 ** File: fsaprs.cc
 ** Description: Source code for the fsaprs class
 ***************************************************************************/

#include "fsaprs.h"

/****************************************************************************
 ** Name: fsAprs
 ** Class: fsAprs
 ** Created:           04/05/2008
 ** Last Modification: 25/05/2008
 ** Args:
 **    Input:
 **       - fsEvaluationAlgorithm fsea: Evaluation Algorithm for Data
 **       - fsDistanceAlgorithm fsda: Distance Algorithm to calculate de
 **                                   distance between two vectors
 **       - fsUInt kFCV: K-Fold Cross Validation value
 **       - fsBool r: Boolean that indicates whether Data is splitted randomly
 **                   or not
 **       - fsChar *fsFileData:  Data file
 **       - fsAlgorithm fsssalg: algorithm to treat the subpartitions if
 **                              necessary at sweep search algorithm
 **       - fsSubPartitionsTreatment fsssspt: for each partition that's not 
 **                                           being optimized at a certain
 **                                           iteration, this variable indicates
 **                                           how it does stay
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsAprs class
 **
 ***************************************************************************/

fsAprs::fsAprs(
             fsAlgorithm fsalg,
             fsEvaluationAlgorithm fsea,
             fsDistanceAlgorithm fsda,
             fsUInt kFCV,
             fsBool r,
             fsChar *fsFileData,
             fsAlgorithm fsssalg,
             fsSubPartitionsTreatment fsssspt,
             fsUInt k)
{
   // Algorithm must has been defined
   assert(fsalg != NOALGORITHM);

   switch(fsalg)
   {
      case GRASEA:
         algorithm = &fsAprs::granularitySearch;
         break;
      case SWESEA:
         algorithm = &fsAprs::sweepSearch;
         break;
      case CODSEA:
         algorithm = &fsAprs::codichotomousSearch;
         break;
      default:
         algorithm = &fsAprs::granularitySearch;
         break;
   }

   // We instantiate the evaluation algorithm
   evaluator = new fsEvaluation(fsea,fsda,kFCV,r,fsFileData);
   
   // Save the values of Sweep Search configuration
   sweep.algorithm = fsssalg;
   sweep.treatment = fsssspt;
   sweep.k         = k;
}

/****************************************************************************
 ** Name: ~fsAprs
 ** Class: fsAprs
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsAprs class
 **
 ***************************************************************************/

fsAprs::~fsAprs()
{
}

/****************************************************************************
 ** Name: granularitySearch
 ** Class: fsAprs
 ** Created:           05/05/2008
 ** Last Modification: 24/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that executes the Granularity Search Algorithm. This
 **              implementation will be made taking profit from the Hill
 **              Climbing algorithm.
 **              The GS algorithm splits X (the set of features) into n/k0
 **              consecutive segments of k0 size. We initialize each
 **              segment randomly.
 **              Afterwards, we assign to each segment a binary variable
 **              named Zj, that represents the fact that we take segment j
 **              or its negation (that consists in negate each one of the
 **              variables of the segment) and let Z the set of Zj. We look
 **              for optimize J in Z.
 **
 ***************************************************************************/

fsVoid fsAprs::granularitySearch()
{
   // We get the number of features
   fsUInt fsAprsFeatures = (*evaluator).numberOfFeatures();
   
   // We declare the vector that will save the features
   // X = {X1,...,Xn}, with Xi belonging to {0,1}
   fsVector X(fsAprsFeatures);

   // We initialize the features vector randomly at the {0,1} dominion,
   // being 0 not taken, and 1 taken feature
   X.randomInitialization(0,1);

   // Estimation value
   fsReal estimation;

   // Best estimation value
   fsReal bestEstimation = -1;
   
   // We create a new vector, Z. Its function is to save the
   // information of the segments for each iteration of the main loop
   // where X data is splitted in ki different segments
   // Each position of Z indicates if the correspondent partition
   // of features is taken as it is (1) or its negation (0)
   fsVector Z;
   
   // We create a new vector that saves the best Z configuration
   fsVector Zbest;
   
   // The granularity value is defined at k
   // Initially there's only one set
   // The use of k will be in the way k=n div (i+1), being n the number of
   // features and i the present iteration
   fsUInt k;

   // This is the main loop, the will increase i, so k will change its
   // value changing the size of the partitions of features.
   // This loop will stop when each partition has only 1 feature
   for(fsUInt i=1;i<fsAprsFeatures;i++)
   {
      // By now, only exact divisions with unoverlapped partitions
      if((fsAprsFeatures%(i+1))==0)
      {
         // We get the size of each partition
         k = fsAprsFeatures/(i+1);
         
         // We construct the vector Z, that will save the 
         Z = fsVector(i+1);

         // Initially, the vector stays as it is right now
         Z.initialization(1.0);
         
         // Zbest it's Z
         Zbest = Z;

         // We get a copy of features vector
         fsVector Xaux;

         // Initially, the best Estimation is when Z has all it's fields equal
         // to 1
         bestEstimation = (*evaluator).J(&Xaux);

         // For each position of Z, we select if we need its current
         // value or the negated one
         for(fsUInt j=0;j<(i+1);j++)
         {
            Z[j] = 0;

            applyBinarySelector(&Z,&Zbest,&Xaux);

            estimation = (*evaluator).J(&Xaux);

            // If we improve the evaluation of X, we save the estimation
            // If we don't, we turn back the changing of Z[j]
            if(estimation > bestEstimation)
            {
               bestEstimation = estimation;
               Zbest = Z;
            }
            else
            {
               // Undo the Xaux changes
               applyBinarySelector(&Z,&Zbest,&Xaux);
               
               // Undo the Z change
               Z[j] = 1;
            }
         }
         
         // At the end of this loop, we've got the best configuration of Z, and
         // so to Xaux, so we assign this one to X
         X = Xaux;
      }
   }
}

/****************************************************************************
 ** Name: sweepSearch
 ** Class: fsAprs
 ** Created:           05/05/2008
 ** Last Modification: 26/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that executes the Sweep Search Algorithm
 **
 ***************************************************************************/

fsVoid fsAprs::sweepSearch()
{
   // We get the number of features
   fsUInt fsAprsFeatures = (*evaluator).numberOfFeatures();
   
   // We must assure the partitions won't be overlapped
   assert((fsAprsFeatures%(sweep.k)) == 0);
   
   // We initialize features set depending on the configuration given
   switch(sweep.treatment)
   {
      case TOONE:    X.initialization(1);         break;
      case TOZERO:   X.initialization(0);         break;
      case RANDOM:   X.randomInitialization(0,1); break;
      case LASTDONE: X.initialization(0);         break;
   }
   
   // We copy the features set solution to an auxiliary one that it's the one
   // we'll rely on to make the correspondent computations
   fsVector Xaux(X);
   
   // We get the number of divisions
   fsUInt numberOfDivisions = (fsAprsFeatures/(sweep.k));
   
   // The Algorithm starts
   for(fsUInt i=0; i<numberOfDivisions; i++)
   {
      // Calling for the method that'll solve the current partition
      switch(sweep.algorithm)
      {
         case SFS:        subSFS(i,sweep.k,&Xaux); break;
         case SBS:        subSBS(i,sweep.k,&Xaux); break;
         case EXHAUSTIVE:
         default:         subExhaustive(i,sweep.k,&Xaux); break;
      }
      
      // Updating the solution, and Xaux in function of the values that have to
      // have the latter's partitionsat the begining of a new iteration
      switch(sweep.treatment)
      {
         case TOONE:
         case TOZERO:
            for(fsUInt j=0;j<sweep.k;j++)
            {
               X[(i*sweep.k)+j] = Xaux[(i*sweep.k)+j];
               Xaux[(i*sweep.k)+j] = (sweep.treatment == TOONE)?1:0;
            }
            
            break;
         case RANDOM:
         case LASTDONE:
             for(fsUInt j=0;j<sweep.k;j++)
            {
               X[(i*sweep.k)+j] = Xaux[(i*sweep.k)+j];
            }
            
            break;
     }
   
      
   }
}

/****************************************************************************
 ** Name: codichotomousSearch
 ** Class: fsAprs
 ** Created:           24/05/2008
 ** Last Modification: 24/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that executes the Co-dichotomous Search Algorithm
 **
 ***************************************************************************/

fsVoid fsAprs::codichotomousSearch()
{
   ;
}


/****************************************************************************
 ** Name: process
 ** Class: fsAprs
 ** Created:           05/05/2008
 ** Last Modification: 24/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that executes the selected feature selection
 **              algorithm in the Adaptive Partitioned Random Search
 **              way
 **
 ***************************************************************************/

fsVoid fsAprs::process()
{
   (this->*algorithm)();
}

/****************************************************************************
 ** Name: applyBinarySelector
 ** Class: fsAprs
 ** Created:           23/05/2008
 ** Last Modification: 24/05/2008
 ** Args:
 **    Input:
 **       - fsVector *Z: Vector that indicates which partitions are taken
 **                      as they are, and which ones are taken by their
 **                      negation
 **       - fsVector *Zbest: Vector that indicates how the partitions are
 **                          taken for the best solutions so far
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that applies a selector of partitions to a vector
 **              of features, considering the best solution. Thereby, the
 **              construction of Xaux from Z will be optimized 'cause
 **              will be change only those partitions where Z and Zbest are
 **              different.
 **              This algorithms optimizes the problem of apply a selector
 **              when size of each partition is huge.
 **
 ***************************************************************************/

fsVoid fsAprs::applyBinarySelector(fsVector *Z,fsVector *Zbest,fsVector *Xaux)
{
   // Get the size of each partition
   fsUInt sizeOfPartitions = (*Xaux).sizeOf()/(*Z).sizeOf();
         
   for(fsUInt i=0;i<(*Z).sizeOf();i++)
   {
      // If two selector vectors are different at this position, we've to
      // change the values of this partition
      if((*Z)[i] != (*Zbest)[i])
      {
         for(fsUInt j=0,offset=i*sizeOfPartitions;j<sizeOfPartitions;j++,offset++)
         {
            // We invert this position
            (*Xaux)[offset] = 1 - (*Xaux)[offset];
         }
      }
   }
   
   return;
} 

/****************************************************************************
 ** Name: subSFS
 ** Class: fsAprs
 ** Created:           25/05/2008
 ** Last Modification: 25/05/2008
 ** Args:
 **    Input:
 **       - fsUInt partition:       Partition to be studied
 **       - fsUInt sizeOfPartition: Size of the partition
 **       - fsVector *Xaux:         Vector of features to be modified
 **    Output:
 **       - fsVector *Xaux:         Vector of features to be modified
 ** Return:
 **    - None
 ** 
 ** Description: Method that executes the Sequential Forward Search
 **              algorithm for a reduced subset of features
 **
 ***************************************************************************/

fsVoid fsAprs::subSFS(fsUInt partition,fsUInt sizeOfPartition,fsVector *Xaux)
{
   //fsUInt fsAprFeatures = (*evaluator).numberOfFeatures();
   fsUInt startPosition = partition*sizeOfPartition;
   
   fsVector subX(sizeOfPartition);
   subX.initialization(0);
   fsReal estimationSubX = -1;
   
   for(fsUInt i=0; i<sizeOfPartition; i++)
      (*Xaux)[i+startPosition] = 0;
   
   for(fsUInt i=0;i<sizeOfPartition; i++)
   {
      fsReal estimation = -1;
      fsUInt posBest;
      fsReal estimationBest = -1;
      
      for(fsUInt j=0; j<sizeOfPartition; j++)
      {
         if((*Xaux)[j+startPosition] == 0)
         {
            (*Xaux)[j+startPosition] = 1;
            estimation = (*evaluator).J(Xaux);
            (*Xaux)[j+startPosition] = 0;
            
            if(estimationBest < estimation)
            {
               estimationBest = estimation;
               posBest = j;
            }
         }
      }
      
      (*Xaux)[posBest+startPosition] = 1;
      
      // Because of the > comparison, we assure that we'll get the best
      // solution with minimum features
      if(estimationBest > estimationSubX)
      {
         subX[posBest] = 1;
         estimationSubX = estimationBest;
      }
   }

   for(fsUInt i=0;i<startPosition;i++)
      (*Xaux)[i+startPosition] = subX[i];
   
   return;
}

/****************************************************************************
 ** Name: subSBS
 ** Class: fsAprs
 ** Created:           25/05/2008
 ** Last Modification: 25/05/2008
 ** Args:
 **    Input:
 **       - fsUInt partition:       Partition to be studied
 **       - fsUInt sizeOfPartition: Size of the partition
 **       - fsVector *Xaux:         Vector of features to be modified
 **    Output:
 **       - fsVector *Xaux:         Vector of features to be modified
 ** Return:
 **    - None
 ** 
 ** Description: Method that executes the Sequential Backward Search
 **              algorithm for a reduced subset of features
 **
 ***************************************************************************/

fsVoid fsAprs::subSBS(fsUInt partition,fsUInt sizeOfPartition,fsVector *Xaux)
{
   //fsUInt fsAprFeatures = (*evaluator).numberOfFeatures();
   fsUInt startPosition = partition*sizeOfPartition;
   
   fsVector subX(sizeOfPartition);
   subX.initialization(1);
   fsReal estimationSubX = -1;
   
   for(fsUInt i=0; i<sizeOfPartition; i++)
      (*Xaux)[i+startPosition] = 1;
   
   for(fsUInt i=0;i<sizeOfPartition; i++)
   {
      fsReal estimation = -1;
      fsUInt posBest;
      fsReal estimationBest = -1;
      
      for(fsUInt j=0; j<sizeOfPartition; j++)
      {
         if((*Xaux)[j+startPosition] == 0)
         {
            (*Xaux)[j+startPosition] = 0;
            estimation = (*evaluator).J(Xaux);
            (*Xaux)[j+startPosition] = 1;
            
            if(estimationBest < estimation)
            {
               estimationBest = estimation;
               posBest = j;
            }
         }
      }
      
      (*Xaux)[posBest+startPosition] = 0;
      
      // Because of the >= comparison, we assure that we'll get the best
      // solution with minimum features
      if(estimationBest >= estimationSubX)
      {
         subX[posBest] = 0;
         estimationSubX = estimationBest;
      }
   }

   for(fsUInt i=0;i<startPosition;i++)
      (*Xaux)[i+startPosition] = subX[i];
   
   return;
}

/****************************************************************************
 ** Name: subExhaustive
 ** Class: fsAprs
 ** Created:           25/05/2008
 ** Last Modification: 26/05/2008
 ** Args:
 **    Input:
 **       - fsUInt partition:       Partition to be studied
 **       - fsUInt sizeOfPartition: Size of the partition
 **       - fsVector *Xaux:         Vector of features to be modified
 **    Output:
 **       - fsVector *Xaux:         Vector of features to be modified
 ** Return:
 **    - None
 ** 
 ** Description: Method that executes the Exhaustive Search
 **              algorithm for a reduced subset of features
 **
 ***************************************************************************/

fsVoid fsAprs::subExhaustive(fsUInt partition,fsUInt sizeOfPartition,fsVector *Xaux)
{
   fsUInt startPosition = partition*sizeOfPartition;
   
   for(fsUInt i=0,position=startPosition;i<sizeOfPartition;i++,position++)
      (*Xaux)[position] = 0;
   
   fsVector Xbest(sizeOfPartition);
   fsReal estimation = -1,estimationBest = -1;
   
   fsUInt nOone = 0;
   
   estimationBest = (*evaluator).J(Xaux);
   Xbest.initialization(0);
   fsUInt k=0;
   
   while(nOone < sizeOfPartition)
   {
      (*Xaux)[startPosition]++;
      
      if((*Xaux)[startPosition] > 1)
      {
         while(k<sizeOfPartition && (*Xaux)[k+startPosition]>1)
         {
            (*Xaux)[k+startPosition] = 0;
            nOone--;
            (*Xaux)[++k+startPosition]++;
         }
      }

      nOone++;
      
      estimation = (*evaluator).J(Xaux);
      
      if(estimation > estimationBest)
      {
         estimationBest = estimation;
         
         for(fsUInt i=0,position=startPosition;i<sizeOfPartition;i++,position++)
         {
            Xbest[i] = (*Xaux)[position];
         }
      }
   }
   
   for(fsUInt i=0,position=startPosition;i<sizeOfPartition;i++,position++)
   {
      (*Xaux)[position] = Xbest[i];
   }
   
   return;
}
