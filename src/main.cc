/***************************************************************************
 *   Copyright (C) 2008 by Antonio Jes�s Adsuar G�mez                      *
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

/***************************************************************************
 ** CORTEX
 ** Tool designed by:
 **    Antonio Jes�s Adsuar G�mez
 **    adsuar@lsi.upc.edu
 ** Sevilla 2008
 **
 **
 ** File: main.cc
 ** Description: Execution manager that allows the launch of different
 **              algorithms of feature selection.
 ***************************************************************************/
#include "main.h"

int main(int argc, char**argv)
{
	fsInt i=1;
	fsUInt k=10;
	fsVector v(k),y;
	
	v[i] = 1992;
	v[3] = 7;
	

	cout << v << endl;
	cout << (y=v) << endl;
	cout << (y==v) << endl;
	
	y[8] = 4;
	
	cout << v << endl;
	cout << y << endl;

	cout << (y==v) << endl;
	
	fsVector X(10);
	
	for(fsInt i=0;i<10;i++)
	   X[i] = 1;
	
	fsDistance d1(SQUAREEUCLIDEANDISTANCE);
	fsDistance d2(EUCLIDEANDISTANCE);

    cout << X << endl;
    
    cout << "d1: " << d1.get(&v,&y,&X) << " ";
    
    cout << "d2: " << d2.get(&v,&y,&X) << endl;
    
    cout << sizeof(fsReal) << endl;
    cout << sizeof(fsInt) << endl;
    cout << sizeof(fsUInt) << endl;
    
    fsReal a = 3e38;
    
    cout << a << endl << endl << endl;

	system("sleep 2");

	return EXIT_SUCCESS;
}

