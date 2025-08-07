#ifndef _MelderThread_h_
#define _MelderThread_h_
/* MelderThread.h
 *
 * Copyright (C) 2014-2018,2020,2025 Paul Boersma
 *
 * This code is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This code is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this work. If not, see <http://www.gnu.org/licenses/>.
 */

#include "melder.h"
#include <vector>
#include <thread>

integer MelderThread_getNumberOfProcessors ();

integer MelderThread_computeNumberOfThreads (
	integer numberOfElements,
	integer thresholdNumberOfElementsPerThread,
	bool useRandom
);

#define MelderThread_DEFINE(threadNumber, firstElement, lastElement)  \
{  \
	std::atomic <bool> _thrown_ = false;  \
	/* We hand all local variables over to the thread lambda by reference, */  \
	/* because many variables (namely those that have to change, such as `_thrown_`, */  \
	/* and those whose copy constructor has been deleted, such as our autoPitch `thee` and our autoVECs `window` and `windowR`) */  \
	/* cannot be copied into the lambda. */  \
	auto _lambda_ = [&] (integer threadNumber, integer firstElement, integer lastElement) {  \
		try {

#define MelderThread_TRY  \
	try {

#define MelderThread_CATCH(errorFlag)  \
	} catch (MelderError) {  \
		errorFlag = true;  \
		return;  \
	}

void MelderThread_run (
	std::atomic <bool> *p_errorFlag,
	integer numberOfElements,
	integer thresholdNumberOfElementsPerThread,
	bool useRandom,
	std::function <void (integer threadNumber, integer firstElement, integer lastElement)> const& threadFunction
);

/* End of file MelderThread.h */
#endif
