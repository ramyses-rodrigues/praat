/* MelderThread.cpp
 *
 * Copyright (C) 2025 Paul Boersma
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

#include "MelderThread.h"

integer MelderThread_getNumberOfProcessors () {
	//return 1;   // un-comment-out to force single-threading
	return Melder_clippedLeft (1_integer, uinteger_to_integer_a (std::thread::hardware_concurrency ()));
}

/* global */ std::thread::id theMelder_error_threadId;

integer MelderThread_computeNumberOfThreads (
	const integer numberOfElements,
	const integer thresholdNumberOfElementsPerThread,
	const bool useRandom)
{
	/* mutable clip */ integer numberOfThreads =
		#if defined (macintosh)
			Melder_iroundDown ((double) numberOfElements / thresholdNumberOfElementsPerThread);
				// round down, assuming that the first spawned thread is the costliest
		#elif defined (_WIN32)
			Melder_iroundDown ((double) numberOfElements / 2.0 / thresholdNumberOfElementsPerThread);
				// round down, assuming that the first spawned thread is the costliest
		#elif defined (linux)
			Melder_iround ((double) numberOfElements / 1.5 / thresholdNumberOfElementsPerThread);
				// round to nearest, assuming that all spawned threads are equally costly
		#else
			#error Undefined platform for MelderThread_computeNumberOfThreads().
		#endif
	Melder_clipRight (& numberOfThreads, MelderThread_getNumberOfProcessors ());
	if (useRandom)
		Melder_clipRight (& numberOfThreads, NUMrandom_maximumNumberOfParallelThreads);
	Melder_clipLeft (1_integer, & numberOfThreads);
	return numberOfThreads;
}

void MelderThread_run (
	std::atomic <bool> *p_errorFlag,
	const integer numberOfElements,
	const integer thresholdNumberOfElementsPerThread,
	const bool useRandom,
	std::function <void (integer threadNumber, integer firstElement, integer lastElement)> const& threadFunction
) {
	const integer numberOfThreads = MelderThread_computeNumberOfThreads (numberOfElements, thresholdNumberOfElementsPerThread, useRandom);
	if (numberOfThreads == 1) {
		threadFunction (0, 1, numberOfElements);
	} else {
		const integer numberOfExtraThreads = numberOfThreads - 1;   // at least 1
		std::vector <std::thread> spawns;
		try {
			spawns. resize (uinteger (numberOfExtraThreads));   // default-construct a number of empty (non-joinable) threads
		} catch (...) {
			Melder_throw (U"Out of memory creating a thread vector. Contact the author if this happens more often.");
		}
		const integer base = numberOfElements / numberOfThreads;
		const integer remainder = numberOfElements % numberOfThreads;
		integer firstElement = 1;
		try {
			for (integer ispawn1 = 1; ispawn1 <= numberOfExtraThreads; ispawn1 ++) {   // ispawn1 is base-1
				const integer lastElement = firstElement + base - 1 + ( ispawn1 <= remainder );
				spawns [uinteger (ispawn1 - 1)] = std::thread (threadFunction, ispawn1, firstElement, lastElement);
				firstElement = lastElement + 1;
			}
		} catch (...) {
			*p_errorFlag = true;   // try to stop any threads that were already spawned
			for (size_t ispawn0 = 0; ispawn0 < spawns.size(); ispawn0 ++)   // ispawn0 is base-0
				if (spawns [ispawn0]. joinable ())   // any extra thread already spawned
					spawns [ispawn0]. join ();   // wait for the spawned thread to finish, hopefully soon
			Melder_throw (U"Couldn't start a thread. Contact the author.");
		}
		Melder_assert (firstElement + base - 1 == numberOfElements);
		threadFunction (0, firstElement, numberOfElements);
		for (size_t ispawn0 = 0; ispawn0 < spawns.size(); ispawn0 ++)   // ispawn0 is base-0
			spawns [ispawn0]. join ();
	}
	if (*p_errorFlag) {
		theMelder_error_threadId = std::this_thread::get_id ();   // TODO: make this truly unique
		throw MelderError();   // turn the error flag back into a MelderError
	}
}

/* End of file MelderThread.cpp */
