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

inline integer MelderThread_getNumberOfProcessors () {
	return Melder_clippedLeft (1_integer, uinteger_to_integer_a (std::thread::hardware_concurrency ()));
}

integer MelderThread_computeNumberOfThreads (
	integer numberOfElements,
	integer thresholdNumberOfElementsPerThread,
	bool useRandom
);

#define MelderThread_BEGIN(numberOfElements, thresholdNumberOfElementsPerThread, useRandom, firstElement, lastElement)  \
{  \
	integer _numberOfElements_ = numberOfElements;   /* Duplicate, because it will be needed in another macro. */  \
	const integer _numberOfThreads_ = MelderThread_computeNumberOfThreads (_numberOfElements_, thresholdNumberOfElementsPerThread, useRandom);  \
	std::atomic <bool> _thrown_ = false;  \
	/* We hand all local variables over to the thread lambda by reference, */  \
	/* because many variables (namely those that have to change, such as `_thrown_`, */  \
	/* and those whose copy constructor has been deleted, such as our autoPitch `thee` and our autoVECs `window` and `windowR`) */  \
	/* cannot be copied into the lambda. */  \
	auto _lambda_ = [&] (integer _ithread_, integer firstElement, integer lastElement) {  \
		try {

#define MelderThread_END  \
		} catch (MelderError) {  \
			_thrown_ = true;  \
			return;  \
		}  \
	};  \
	if (_numberOfThreads_ == 1) {  \
		_lambda_ (MelderThread_MASTER, 1, _numberOfElements_);  \
	} else {  \
		const integer _numberOfExtraThreads = _numberOfThreads_ - 1;   /* At least 1. */  \
		std::vector <std::thread> _spawns { uinteger (_numberOfExtraThreads) };   /* Safe cast. */  \
		const integer _base = _numberOfElements_ / _numberOfThreads_;  \
		const integer _remainder = _numberOfElements_ % _numberOfThreads_;  \
		integer _firstElement = 1;  \
		try {  \
			for (integer _ispawn1 = 1; _ispawn1 <= _numberOfExtraThreads; _ispawn1 ++) {   /* _ispawn1 is base-1 */  \
				const integer _lastElement = _firstElement + _base - 1 + ( _ispawn1 <= _remainder );  \
				_spawns [uinteger (_ispawn1 - 1)] = std::thread (_lambda_, _ispawn1, _firstElement, _lastElement);  \
				_firstElement = _lastElement + 1;  \
			}  \
		} catch (...) {  \
			_thrown_ = true;   /* Try to stop any threads that were already spawned. */  \
			for (size_t _ispawn0 = 0; _ispawn0 < _spawns.size(); _ispawn0 ++)   /* _ispawn0 is base-0 */  \
				if (_spawns [_ispawn0]. joinable ())   /* Any extra thread already spawned. */  \
					_spawns [_ispawn0]. join ();   /* Wait for the spawned thread to finish, hopefully soon. */  \
			Melder_throw (U"Couldn't start a thread. Contact the author.");  \
		}  \
		Melder_assert (_firstElement + _base - 1 == _numberOfElements_);  \
		_lambda_ (MelderThread_MASTER, _firstElement, _numberOfElements_);  \
		for (size_t _ispawn0 = 0; _ispawn0 < _spawns.size(); _ispawn0 ++)   /* _ispawn0 is base-0 */  \
			_spawns [_ispawn0]. join ();  \
	}  \
	if (_thrown_) {  \
		theMelder_error_threadId = std::this_thread::get_id ();  \
		throw MelderError();  \
	}  \
}

/*
	In sound analysis, the following is typically called at the beginning of each frame
	(between MelderThread_BEGIN and MelderThread_END).
*/
#define MelderThread_OPPORTUNITY_TO_BAIL_OUT  \
	if (_thrown_)  \
		return;   // stop the thread

/*
	The following two definitions are typically used to do something in the master thread,
	such as a visualization or an interaction with the GUI
	(between MelderThread_BEGIN and MelderThread_END):

		if (MelderThread_CURRENT == MelderThread_MASTER)
			Melder_progress (...);

	or when debugging (between MelderThread_BEGIN and MelderThread_END):

		trace ("Testing thread ", MelderThread_CURRENT, U".");

	The value of MelderThread_CURRENT is either 0 (the master thread), or a number between 1 and the number of spawned threads.
*/
#define MelderThread_CURRENT  _ithread_
#define MelderThread_MASTER  0

/* End of file MelderThread.h */
#endif
