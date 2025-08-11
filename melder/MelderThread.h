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

void MelderThread_run (
	std::atomic <bool> *p_errorFlag,
	integer numberOfElements,
	integer thresholdNumberOfElementsPerThread,
	bool useRandom,
	std::function <void (integer threadNumber, integer firstElement, integer lastElement)> const& threadFunction
);

/*
	Here is how these functions help you create simple parallel procedures.
	Our example is pitch analysis in Praat (see fon/Sound_to_Pitch.cpp for details).

	You start with defining the function:

		autoPitch Sound_to_Pitch (Sound me, integer maxnCandidates, ...)
		{

	The function can throw:
			try {

	Then you create the result of the function:

				autoPitch thee = Pitch_create (...);

	Then you create the intermediate objects that are read-only for all threads
	(so that you have to create them only once):

				autoVEC window = ...;   // the window shape
				autoVEC windowR = ...;   // the autocorrelation of the window shape

	You can open a progress bar if you want:

				autoMelderProgress progress (U"Sound to Pitch...");

	Then you define the "thread function", which is a lambda that is called for every thread.
	The thread function is not allowed to throw exceptions, because exceptions cannot pass
	directly from a newly spawned thread into the master thread. Therefore,
	any exception generated with the thread function has to translate each MelderError
	into an error flag:

				std::atomic <bool> errorFlag = false;

	You then define the thread function as a lambda.
	We hand all local variables over to the thread function by reference,
	because many variables (namely those that have to change, such as `errorFlag`,
	and those whose copy constructor has been deleted,
	such as our autoPitch `thee` and our autoVECs `window` and `windowR`,
	cannot be copied into the lambda.

				auto Sound_to_Pitch_threadFunction = [&] (integer threadNumber, integer firstFrame, integer lastFrame) {

	We could also have said [=, & errorFlag, & thee, & window, & windowR],
	whereby most variables would have been handed over to the lambda by value,
	but simply handing over all variables by reference is more straightforward and not really slower.

	The thread function can throw, so you need a `try`-`catch` pair:
					try {

	Inside the thread function you first create the objects that are different for each thread, such as buffers:

						autoMAT frame = zero_MAT (my ny, ...);
						autoNUMFourierTable fftTable = NUMFourierTable_create (...);
						autoVEC ac = zero_VEC (...);
						autoVEC rbuffer = zero_VEC (...);
						double *r = & rbuffer [...];
						autoINTVEC imax = zero_INTVEC (maxnCandidates);
						autoVEC localMean = zero_VEC (my ny);

	You can now cycle over the frames that belong to the thread:

						for (integer iframe = firstFrame; iframe <= lastFrame; iframe ++) {

	Each frame provides an opportunity to bail out from the thread function (and therefore to abort the thread):

							if (errorFlag)
								return;   // abort this thread

	If not aborted, you can now do the pitch analysis for the present frame:
	
							Pitch_Frame pitchFrame = & thy frames [iframe];
							Sound_into_PitchFrame (me, pitchFrame, maxnCandidates,
								window.get(), windowR.get(),
								frame.get(), fftTable.get(), ac.get(),
								r, imax.get(), localMean.get()
							);

	You can update the progress bar if you are in the master thread, which has threadNumber 0.
	You can do that anywhere in the frame loop; here we do it at the end of it:

							if (threadNumber == 0) {
								const double estimatedFractionAnalysed = (iframe - firstFrame + 0.5) / (lastFrame - firstFrame + 1.0);
								Melder_progress (0.1 + 0.8 * estimatedFractionAnalysed,
									U"Sound to Pitch: analysed approximately ", Melder_iround (numberOfFrames * estimatedFractionAnalysed),
									U" out of ", numberOfFrames, U" frames"
								);
							}
						}

	The `catch` differs from the normal "catch (MelderError) { Melder_throw (...); }",
	in that it has to convert the MelderError into an error flag:

					} catch (MelderError) {
						errorFlag = true;
						return;
					}

	Finally, you close the lambda and call it:

				};
				MelderThread_run (& errorFlag, numberOfFrames, 5, false, Sound_to_Pitch_threadFunction);

	You also close the progress bar and return the result:

				Melder_progress (0.95, U"Sound to Pitch: path finder");
				return thee;

	You end with the `catch` of Sound_to_Pitch itself:

			} catch (MelderError) {
				Melder_throw (me, U": pitch analysis not performed.");
			}
		}
*/

/*
	The above procedure can be simplified a bit by using the following four macros.
	We surround the code by an extra pair of braces in order to allow multiple use within a function:
*/

#define MelderThread_PARALLELIZE(numberOfElements, thresholdNumberOfElementsPerThread, useRandom, _ithread)  \
	{/* start of scope of `_errorFlag_` and `_threadFunction_` */  \
		const integer _numberOfElements_ = numberOfElements;  \
		const integer _thresholdNumberOfElementsPerThread_ = thresholdNumberOfElementsPerThread;  \
		const bool _useRandom_ = useRandom;  \
		std::atomic <bool> _errorFlag_ = false;  \
		auto _threadFunction_ = [&] (integer _ithread, integer _firstElement_, integer _lastElement_) {  \
			try {

#define MelderThread_FOR(ielement)  \
				for (integer ielement = _firstElement_; ielement <= _lastElement_; ielement ++) {  \
					if (_errorFlag_)  \
						return;

#define MelderThread_ENDFOR  \
				}  \
			} catch (MelderError) {  \
				_errorFlag_ = true;  \
				return;  \
			}  \
		};  \
		MelderThread_run (& _errorFlag_, _numberOfElements_, _thresholdNumberOfElementsPerThread_, _useRandom_, _threadFunction_);  \
	}/* end of scope of `_errorFlag_` and `_threadFunction_` */

#define MelderThread_GET_ESTIMATED_FRACTION_ANALYSED(ielement)  \
	((ielement) - _firstElement_ + 0.5) / (_lastElement_ - _firstElement_ + 1.0)

/*
	The whole Sound_to_Pitch then becomes:

		autoPitch Sound_to_Pitch (Sound me, integer maxnCandidates, ...)
		{
			try {
				autoPitch thee = Pitch_create (...);
				autoVEC window = ...;   // the window shape
				autoVEC windowR = ...;   // the autocorrelation of the window shape

				autoMelderProgress progress (U"Sound to Pitch...");

				MelderThread_PARALLELIZE (numberOfFrames, 5, false, threadNumber)

				autoMAT frame = zero_MAT (my ny, ...);
				autoNUMFourierTable fftTable = NUMFourierTable_create (...);
				autoVEC ac = zero_VEC (...);
				autoVEC rbuffer = zero_VEC (...);
				double *r = & rbuffer [...];
				autoINTVEC imax = zero_INTVEC (maxnCandidates);
				autoVEC localMean = zero_VEC (my ny);

				MelderThread_FOR (iframe)

					Pitch_Frame pitchFrame = & thy frames [iframe];
					Sound_into_PitchFrame (me, pitchFrame, maxnCandidates,
						window.get(), windowR.get(),
						frame.get(), fftTable.get(), ac.get(),
						r, imax.get(), localMean.get()
					);

					if (threadNumber == 0) {   // are we in the master thread? then we can interact with the GUI
						const double estimatedFractionAnalysed = MelderThread_GET_ESTIMATED_FRACTION_ANALYSED (iframe);
						Melder_progress (0.1 + 0.8 * estimatedFractionAnalysed,
							U"Sound to Pitch: analysed approximately ", Melder_iround (numberOfFrames * estimatedFractionAnalysed),
							U" out of ", numberOfFrames, U" frames"
						);
					}

				MelderThread_ENDFOR

				Melder_progress (0.95, U"Sound to Pitch: path finder");
				return thee;
			} catch (MelderError) {
				Melder_throw (me, U": pitch analysis not performed.");
			}
		}
*/

/* End of file MelderThread.h */
#endif
