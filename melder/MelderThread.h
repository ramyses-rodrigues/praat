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
#include <atomic>

integer MelderThread_getNumberOfProcessors ();

integer MelderThread_computeNumberOfThreads (
	integer numberOfElements,
	integer thresholdNumberOfElementsPerThread,
	bool useRandom
);

/*
	The following is a replacement for std::this_thread::get_id(),
	which is not guaranteeed to return a unique ID.

	A case where this matters is if a throwing thread has finished,
	and then a new thread is created with the same idea and quickly throws as well.
	Admittedly, this is rare, but it is not impossible and it is easy to prevent it.
*/
integer Melder_thisThread_getUniqueID ();

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
