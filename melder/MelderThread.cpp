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

/* global */ std::thread::id theMelder_error_threadId;

integer MelderThread_computeNumberOfThreads (
	const integer numberOfElements,
	const integer thresholdNumberOfElementsPerThread,
	const bool useRandom)
{
	/* mutable clip */ integer numberOfThreads = Melder_iroundDown ((double) numberOfElements / thresholdNumberOfElementsPerThread);
			// round down, assuming that the first spawned thread is the costiest
	Melder_clipRight (& numberOfThreads, MelderThread_getNumberOfProcessors ());
	if (useRandom)
		Melder_clipRight (& numberOfThreads, NUMrandom_maximumNumberOfParallelThreads);
	Melder_clipLeft (1_integer, & numberOfThreads);
	return numberOfThreads;
}


/* End of file MelderThread.cpp */
