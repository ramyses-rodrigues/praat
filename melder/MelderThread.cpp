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

#include "melder.h"
#include <vector>
#include <thread>

integer MelderThread_getNumberOfProcessors () {
	//return 1;   // un-comment-out to force single-threading
	return Melder_clippedLeft (1_integer, uinteger_to_integer_a (std::thread::hardware_concurrency ()));
}

/* global */ std::atomic <integer> theMelder_error_threadId;

integer Melder_thisThread_getUniqueID () {
	static std::atomic <integer> uniqueID = 0;
	static thread_local integer thisThread_uniqueID = uniqueID ++;
	//TRACE
	trace (thisThread_uniqueID);
	return thisThread_uniqueID;
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
		theMelder_error_threadId = Melder_thisThread_getUniqueID ();
		throw MelderError();   // turn the error flag back into a MelderError
	}
}

integer MelderThread_computeNumberOfThreads (
	const integer numberOfElements,
	const integer thresholdNumberOfElementsPerThread,
	const bool useRandom)
{
	if (! MelderThread_getUseMultithreading ())
		return 1;
	integer minimumNumberOfElementsPerThread = MelderThread_getMinimumNumberOfElementsPerThread ();
	if (minimumNumberOfElementsPerThread <= 0)
		minimumNumberOfElementsPerThread = thresholdNumberOfElementsPerThread;
	/* mutable clip */ integer numberOfThreads =
		#if defined (macintosh)
			Melder_iroundDown ((double) numberOfElements / minimumNumberOfElementsPerThread);
				// round down, assuming that the first spawned thread is the costliest
		#elif defined (_WIN32)
			Melder_iroundDown ((double) numberOfElements / 2.0 / minimumNumberOfElementsPerThread);
				// round down, assuming that the first spawned thread is the costliest
		#elif defined (linux)
			Melder_iround ((double) numberOfElements / 1.5 / minimumNumberOfElementsPerThread);
				// round to nearest, assuming that all spawned threads are equally costly
		#else
			#error Undefined platform for MelderThread_computeNumberOfThreads().
		#endif
	Melder_clipRight (& numberOfThreads, MelderThread_getMaximumNumberOfConcurrentThreads ());
	if (useRandom)
		Melder_clipRight (& numberOfThreads, NUMrandom_maximumNumberOfParallelThreads);
	Melder_clipLeft (1_integer, & numberOfThreads);
	return numberOfThreads;
}

void MelderThread_getInfo (integer numberOfElements, integer *p_numberOfThreads, integer *p_numberOfElementsPerThread) {
	const integer maximumNumberOfConcurrentThreads = MelderThread_getMaximumNumberOfConcurrentThreads ();
	integer minimumNumberOfElementsPerThread = MelderThread_getMinimumNumberOfElementsPerThread ();
	if (minimumNumberOfElementsPerThread <= 0)
		minimumNumberOfElementsPerThread = 40;   // BUG: hard-coded here, but should be factory-tuned
	const integer maximumNumberOfElementsPerThread =
			Melder_clippedLeft (minimumNumberOfElementsPerThread, MelderThread_getMaximumNumberOfElementsPerThread ());
	Melder_assert (maximumNumberOfConcurrentThreads > 0);
	Melder_assert (minimumNumberOfElementsPerThread > 0);
	Melder_assert (maximumNumberOfElementsPerThread >= minimumNumberOfElementsPerThread);
	if (MelderThread_getUseMultithreading ()) {
		*p_numberOfElementsPerThread = Melder_iroundUp ((double) numberOfElements / maximumNumberOfConcurrentThreads);
		Melder_clip (minimumNumberOfElementsPerThread, p_numberOfElementsPerThread, maximumNumberOfElementsPerThread);
		*p_numberOfThreads = Melder_iroundUp ((double) numberOfElements / *p_numberOfElementsPerThread);
		Melder_clipLeft (1_integer, p_numberOfThreads);
		//TRACE
		trace (numberOfElements, U" ", maximumNumberOfConcurrentThreads, U" ", *p_numberOfElementsPerThread, U" ", *p_numberOfThreads);
	} else {
		*p_numberOfThreads = 1;
		*p_numberOfElementsPerThread = numberOfElements;
	}
	Melder_assert (*p_numberOfThreads > 0);
	Melder_assert (*p_numberOfElementsPerThread > 0 || numberOfElements == 0);   // note edge case
}

/*
	Preferences.
*/

static struct {
	bool useMultithreading = true;
	integer maximumNumberOfConcurrentThreads = 0;   // "0" signals automatic
	integer minimumNumberOfElementsPerThread = 0;   // "0" signals the factory-tuned value
	integer maximumNumberOfElementsPerThread = 0;   // "0" signals no limit
	bool traceThreads = false;
} preferences;

void MelderThread_debugMultithreading (bool useMultithreading, integer maximumNumberOfConcurrentThreads,
	integer minimumNumberOfElementsPerThread, integer maximumNumberOfElementsPerThread, bool traceThreads)
{
	preferences. useMultithreading = useMultithreading;
	preferences. maximumNumberOfConcurrentThreads = maximumNumberOfConcurrentThreads;
	preferences. minimumNumberOfElementsPerThread = minimumNumberOfElementsPerThread;
	preferences. maximumNumberOfElementsPerThread = maximumNumberOfElementsPerThread;
	preferences. traceThreads = traceThreads;
}

bool MelderThread_getUseMultithreading () {
	return preferences. useMultithreading;
}

integer MelderThread_getMaximumNumberOfConcurrentThreads () {
	if (! preferences. useMultithreading)
		return 1;
	if (preferences. maximumNumberOfConcurrentThreads <= 0)
		return MelderThread_getNumberOfProcessors ();
	return preferences. maximumNumberOfConcurrentThreads;
}

integer MelderThread_getMinimumNumberOfElementsPerThread () {
	if (! preferences. useMultithreading)
		return 1;
	if (preferences. minimumNumberOfElementsPerThread <= 0)
		return 0;   // signals factory tuning
	return preferences. minimumNumberOfElementsPerThread;
}

integer MelderThread_getMaximumNumberOfElementsPerThread () {
	if (! preferences. useMultithreading)
		return INTEGER_MAX;
	if (preferences. maximumNumberOfElementsPerThread <= 0)
		return INTEGER_MAX;   // signals no limit
	return preferences. maximumNumberOfElementsPerThread;
}

bool MelderThread_getTraceThreads () {
	return preferences. traceThreads;
}

/* End of file MelderThread.cpp */
