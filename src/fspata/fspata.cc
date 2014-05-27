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
 ** Description: Source code for the fspata class
 ***************************************************************************/

#include "fspata.h"

/****************************************************************************
 ** Name: fsPata
 ** Class: fsPata
 ** Created:           05/05/2008
 ** Last Modification: 18/05/2008
 ** Args:
 **    Input:
 **       - fsEvaluationAlgorithm fsea: Evaluation Algorithm for Data
 **       - fsDistanceAlgorithm fsda: Distance Algorithm to calculate de
 **                                   distance between two vectors
 **       - fsUInt kFCV: K-Fold Cross Validation value
 **       - fsBool r: Boolean that indicates whether Data is splitted randomly
 **                   or not
 **       - fsChar *fsFileData:  Data file
 **       - fsInt sf:  Steps Forward
 **       - fsInt sb:  Steps Backward
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsPata class
 **
 ***************************************************************************/

fsPata::fsPata(
             fsAlgorithm fsalg,
             fsEvaluationAlgorithm fsea,
             fsDistanceAlgorithm fsda,
             fsUInt kFCV,
             fsBool r,
             fsChar *fsFileData,
             fsInt sf,
             fsInt sb)
{
   // Algorithm must has been defined
   assert(fsalg != NOALGORITHM);

   switch(fsalg)
   {
      case SFS:
         algorithm = &fsPata::sfs;
         break;
      case SBS:
         algorithm = &fsPata::sbs;
         break;
      case PATA:
         algorithm = &fsPata::pata;

         // Steps can't be both 0
         assert(!(sf==0 && sb==0));

         stepsForward = sf;
         stepsBackward = sb;
         break;
      case SFFS:
         algorithm = &fsPata::sffs;
         break;
      case SBFS:
         algorithm = &fsPata::sbfs;
         break;
      default:
         algorithm = &fsPata::sfs;
         break;
   }

   // We reserve the memory for different execution threads
   threads = new fsPD[kTHREADS];
   
   // Executed Js
   executedJs = 0;
   
   // Add of all J
   addOfAllJ = 0;

   // We open the file to save Js
   alg_plot_J.open("Js.txt",ios::out);
   
   // We initialize main_k, that's the value of the current k of the problem
   main_k = 0;
}

/****************************************************************************
 ** Name: ~fsPata
 ** Class: fsPata
 ** Created:           05/05/2008
 ** Last Modification: 18/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsPata class
 **
 ***************************************************************************/

fsPata::~fsPata()
{
   // We close the file
   alg_plot_J.close();
}

/****************************************************************************
 ** Name: sfs
 ** Class: fsPata
 ** Created:           05/05/2008
 ** Last Modification: 18/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that executes the Sequential Forward Search
 **              Algorithm
 **
 ***************************************************************************/

fsVoid fsPata::sfs()
{
   main_k = 0;
   
   // Take a look at Boundary_K
   
   while(main_k<evaluator.numberOfFeatures())
   {
      goForward(1,0);
   }
}

/****************************************************************************
 ** Name: sbs
 ** Class: fsPata
 ** Created:           05/05/2008
 ** Last Modification: 18/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that executes the Sequential Backward Search
 **              Algorithm
 **
 ***************************************************************************/

fsVoid fsPata::sbs()
{
   main_k = evaluator.numberOfFeatures();
   
   for(fsUInt i=0;i<evaluator.numberOfFeatures();i++)
   {
      // We select all the features
      threads[0].X[main_k][i] = 1;
   }
   
   // Take a look at Boundary_K
   
   while(main_k>1)
   {
      goBackward(1,0);
   }
}


/****************************************************************************
 ** Name: sffs
 ** Class: fsPata
 ** Created:           07/05/2008
 ** Last Modification: 07/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that executes the Sequential Forward Floating
 **              Search Algorithm
 **
 ***************************************************************************/

fsVoid fsPata::sffs()
{
}


/****************************************************************************
 ** Name: sbfs
 ** Class: fsPata
 ** Created:           07/05/2008
 ** Last Modification: 07/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that executes the Sequential Backward Floating
 **              Search Algorithm
 **
 ***************************************************************************/

fsVoid fsPata::sbfs()
{
}


/****************************************************************************
 ** Name: pata
 ** Class: fsPata
 ** Created:           07/05/2008
 ** Last Modification: 07/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that executes Generic Plus And Take Away Algorithms
 **
 ***************************************************************************/

fsVoid fsPata::pata()
{
}


/****************************************************************************
 ** Name: process
 ** Class: fsPata
 ** Created:           05/05/2008
 ** Last Modification: 05/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that executes the selected feature selection
 **              algorithm in the Plus And Take Away Search way
 **
 ***************************************************************************/

fsVoid fsPata::process()
{
   (this->*algorithm)();
}

/****************************************************************************
 ** Name: bestFeatureToAdd
 ** Class: fsPata
 ** Created:           08/05/2008
 ** Last Modification: 18/05/2008
 ** Args:
 **    Input:
 **       - fsUInt thrd: Thread in which make the calculation
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that gets the best feature of which are not in the
 **              solution
 **
 ***************************************************************************/

fsVoid fsPata::bestFeatureToAdd(fsUInt thrd)
{
   fsReal estimation;
   
   fsUInt k = threads[thrd].k;
   
   threads[thrd].add_est = -1.0;

   fsUInt features = X.sizeOf();

   for(fsUInt i=0;i<features;i++)
   {
      // If the i feature isn't selected yet, we see what happens if it's added
      if(X[i] == 0)
      {
         // We add the i feature
         X[i] = 1;
         // We estimate the value of the new solution
         estimation = evaluator.J(&X);

         // We update the variables of the system
         J_best(estimation,thrd);
         J_post(estimation,k+1);
         
         // We substract the i feature
         X[i] = 0;
         
         if(estimation>threads[thrd].add_est)
         {
            threads[thrd].add_est = estimation;
            threads[thrd].add = i;
         }
      }
   }
}

/****************************************************************************
 ** Name: bestFeatureToSubstract
 ** Class: fsPata
 ** Created:           08/05/2008
 ** Last Modification: 18/05/2008
 ** Args:
 **    Input:
 **       - fsUInt thrd: Thread in which make the calculation
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that gets the worst feature of which are in the
 **              solution
 **
 ***************************************************************************/

fsVoid fsPata::bestFeatureToSubstract(fsUInt thrd)
{
   fsReal estimation;
   
   threads[thrd].sub_est = -1.0;

   fsUInt k = threads[thrd].k;
   
   fsUInt features = X.sizeOf();

   for(fsUInt i=0;i<features;i++)
   {
      // If the i feature isn't selected yet, we see what happens if it's added
      if(X[i] == 1)
      {
         // We add the i feature
         X[i] = 0;
         // We estimate the value of the new solution
         estimation = evaluator.J(&X);
         
         // We update the variables of the system
         J_best(estimation,thrd);
         J_post(estimation,k-1);
         
         // We substract the i feature
         X[i] = 1;
         
         if(estimation>threads[thrd].sub_est)
         {
            threads[thrd].sub_est = estimation;
            threads[thrd].sub = i;
         }
      }
   }
}

/****************************************************************************
 ** Name: goForward
 ** Class: fsPata
 ** Created:           14/05/2008
 ** Last Modification: 18/05/2008
 ** Args:
 **    Input:
 **       - fsUInt steps: steps to go forward
 **       - fsUInt thrd: number of thread
 **    Output:
 **       - None
 ** Return:
 **    - fsUInt: Effective steps forward
 ** 
 ** Description: Method makes a certain number of steps forward
 **
 ***************************************************************************/

fsUInt fsPata::goForward(fsUInt steps,fsUInt thrd)
{
   fsUInt effectiveSteps = 0;
    
   cout << "=====================================================" << endl;    
   cout << "Execution of " << steps << " forward step(s)" << endl;
   cout << "=====================================================" << endl;    

   // We'll do a maximum of "steps" loops
   while(effectiveSteps<steps)
   {
      // We get the best feature to be added at 'add', saving it's estimation at
      // 'add_est'
      (threads[thrd].add_est) = -1;
      
      if(Boundary_K(threads[thrd].k+1))
      {
         //cdebug_2("El valor de k en la salida es ",threads[opc].k);
         return effectiveSteps;
      }
       
      bestFeatureToSubstract(thrd);

      // We substract the feature selected by bestFeatureToSubstract
      threads[thrd].vJ[(threads[thrd].k)+1]=(threads[thrd].add_est);
      
      // The result set with (k-1) features, would be made up with the
      // current features minus the worst feature, given by
      // bestFeatureToSubstract
      (threads[thrd].X[(threads[thrd].k)+1])[threads[thrd].add] = 1; 

      //cout << threads[opc].k << endl;
      threads[thrd].k = threads[thrd].k + 1;
      //cout << threads[opc].k << endl;
       
      threads[thrd].last = threads[thrd].sub;
      
      // Incrementamos los pasos efectivos
      //next_step_exec();
       
      incK();
       
      effectiveSteps++;
       
      // Update the best J once the step is over
      J_best(threads[thrd].vJ[threads[thrd].k],thrd);
       
      //getOut(-effectiveSteps);
   }
    
   return effectiveSteps;
}

/****************************************************************************
 ** Name: goBackward
 ** Class: fsPata
 ** Created:           14/05/2008
 ** Last Modification: 18/05/2008
 ** Args:
 **    Input:
 **       - fsUInt steps: steps to go backward
 **       - fsUInt thrd: number of thread
 **    Output:
 **       - None
 ** Return:
 **    - fsUInt: Effective steps backward
 ** 
 ** Description: Method makes a certain number of steps backward
 **
 ***************************************************************************/

fsUInt fsPata::goBackward(fsUInt steps,fsUInt thrd)
{
   fsUInt effectiveSteps = 0;
    
   cout << "=====================================================" << endl;    
   cout << "Execution of " << steps << " backward step(s)" << endl;
   cout << "=====================================================" << endl;    

   // We'll do a maximum of "steps" loops
   while(effectiveSteps<steps)
   {
      // We get the best feature to be substracted at 'sub', saving it's
      // estimation at 'sub_est'
      (threads[thrd].sub_est) = -1;
      
      if(Boundary_K(threads[thrd].k-1))
      {
         //cdebug_2("El valor de k en la salida es ",threads[opc].k);
         return -effectiveSteps;
      }
       
      bestFeatureToSubstract(thrd);

      // We substract the feature selected by bestFeatureToSubstract
      threads[thrd].vJ[(threads[thrd].k)-1]=(threads[thrd].sub_est);
      
      // The result set with (k-1) features, would be made up with the
      // current features minus the worst feature, given by
      // bestFeatureToSubstract
      (threads[thrd].X[(threads[thrd].k)-1])[threads[thrd].sub] = 0; 

      //cout << threads[opc].k << endl;
      threads[thrd].k = threads[thrd].k - 1;
      //cout << threads[opc].k << endl;
       
      threads[thrd].last = threads[thrd].sub;
      
      // Incrementamos los pasos efectivos
      //next_step_exec();
       
      decK();
       
      effectiveSteps++;
       
      // Update the best J once the step is over
      J_best(threads[thrd].vJ[threads[thrd].k],thrd);
       
      //getOut(-effectiveSteps);
   }
    
   return -effectiveSteps;
}
