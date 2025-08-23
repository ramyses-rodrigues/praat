/* SampledIntoSampled.cpp
 *
 * Copyright (C) 2024,2025 David Weenink, 2025 Paul Boersma
 *
 * This code is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This code is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this work. If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>
#include <thread>

#include "SampledIntoSampled.h"
#include "Sound_and_LPC.h"
#include "Sound_extensions.h"
#include "NUM2.h"

#include "oo_DESTROY.h"
#include "SampledIntoSampled_def.h"
#include "oo_COPY.h"
#include "SampledIntoSampled_def.h"
#include "oo_EQUAL.h"
#include "SampledIntoSampled_def.h"
#include "oo_CAN_WRITE_AS_ENCODING.h"
#include "SampledIntoSampled_def.h"
#include "oo_WRITE_TEXT.h"
#include "SampledIntoSampled_def.h"
#include "oo_WRITE_BINARY.h"
#include "SampledIntoSampled_def.h"
#include "oo_READ_TEXT.h"
#include "SampledIntoSampled_def.h"
#include "oo_READ_BINARY.h"
#include "SampledIntoSampled_def.h"
#include "oo_DESCRIPTION.h"
#include "SampledIntoSampled_def.h"

Thing_implement (SampledIntoSampled, Daata, 0);

void SampledIntoSampled_init (mutableSampledIntoSampled me, constSampled input, mutableSampled output) {
	SampledIntoSampled_assertEqualDomains (input, output);
	my input = input;
	my output = output;
}

autoSampledIntoSampled SampledIntoSampled_create (constSampled input, mutableSampled output, autoSampledFrameIntoSampledFrame ws,
	autoSampledIntoSampledStatus status)
{
	try {
		autoSampledIntoSampled me = Thing_new (SampledIntoSampled);
		SampledIntoSampled_init (me.get(), input, output);
		my frameIntoFrame = ws.move();
		my frameIntoFrame -> allocateOutputFrames ();
		my status = status.move();
		const bool updateStatus = MelderThread_getTraceThreads ();
		SampledFrameIntoSampledFrame_initForStatusUpdates (my frameIntoFrame.get(), my status.get(), updateStatus);
		return me;
	} catch (MelderError) {
		Melder_throw (U"SampledIntoSampled not created.");
	}
}

integer SampledIntoSampled_analyseThreaded (mutableSampledIntoSampled me) {
	SampledFrameIntoSampledFrame frameIntoFrame = my frameIntoFrame.get();
	std::atomic<integer> globalFrameErrorCount (0);   // TODO: remove, because errors are notified, not counted
	try {
		std::atomic <bool> errorFlag = false;
		auto analyseFrames = [&] (int threadNumber, integer fromFrame, integer toFrame) {
			try {
				NUMrandom_setChannel (threadNumber);
				Melder_thisThread_setRange (fromFrame, toFrame);   // for progress bars in inner functions; otherwise not needed
				autoSampledFrameIntoSampledFrame frameIntoFrameCopy = Data_copy (frameIntoFrame);   // can throw MelderError
				frameIntoFrameCopy -> startFrame = fromFrame;   // TODO: remove
				for (integer iframe = fromFrame; iframe <= toFrame; iframe ++) {
					if (errorFlag)
						return;   // abort this thread
					Melder_thisThread_setCurrentElement (iframe);   // for progress bars in inner functions; otherwise not needed
					frameIntoFrameCopy -> currentFrame = iframe;
					frameIntoFrameCopy -> getInputFrame ();
					if (! frameIntoFrameCopy -> inputFrameToOutputFrame ())
						frameIntoFrameCopy -> framesErrorCount ++;   // TODO: remove
					frameIntoFrameCopy -> saveOutputFrame ();
					if (frameIntoFrameCopy -> updateStatus)   // TODO: remove
						frameIntoFrameCopy -> status -> frameIntoFrameInfo [frameIntoFrameCopy -> currentFrame] =
								frameIntoFrameCopy -> frameAnalysisInfo;   // TODO: remove
				}
				globalFrameErrorCount += frameIntoFrameCopy -> framesErrorCount;   // TODO: remove
			} catch (MelderError) {
				errorFlag = true;   // convert the MelderError to an error flag temporarily (after building up notification)
				return;   // abort the thread
			}
		};
		MelderThread_run (& errorFlag, my output -> nx, 40, analyseFrames);
	} catch (MelderError) {
		if (frameIntoFrame -> updateStatus)   // TODO: remove
			my status -> showStatus ();   // TODO: remove
		Melder_throw (me, U"The Sampled analysis could not be done.");
	}
	if (frameIntoFrame -> updateStatus)   // TODO: remove
		my status -> showStatus ();   // TODO: remove
	return globalFrameErrorCount;   // TODO: remove

	/*
		TODO: once globalFrameErrorCount has gone, we will be able to use the MelderThread macros
	*/
}

/*
	Performs timing of a number of scenarios for multi-threading.
	This timing is performed on the LPC analysis with the Burg algorithm on a sound file of a given duration
	and a sampling frequency of 11000 Hz.
	The workspace for the Burg algorithm needs more memory for its analyses than the other LPC algorithms (it needs
	n samples for the windowed sound frame and at least 2 vectors of length n for buffering).
	It varies the number of threads from 1 to the maximum number of concurrency available on the hardware.
	It varies, for each number of threads separately, the frame sizes (50, 100, 200, 400, 800, 1600, 3200)
	The data is represented in the info window as a space separated table with 4 columns:
	duration(s) nThread nFrames/thread toLPC(s)
	Saving this data, except for the last line, as a csv file and next reading this file as a Table,
	the best way to show the results would be
	Table > Scatter plot: "nFrames/thread", 0, 0, toLPC(s), 0, 0, nThread, 8, "yes"
*/
void SampledIntoSampled_timeMultiThreading (double soundDuration) {
	/*
		Save current multi-threading situation
	*/
	try {
		autoVEC framesPerThread { 1, 10, 20, 30, 40, 50, 70, 100, 200, 400, 800, 1600, 3200 };
		const integer maximumNumberOfThreads = MelderThread_getNumberOfProcessors ();
		autoSound me = Sound_createSimple (1_integer, soundDuration, 5500.0);
		for (integer i = 1; i <= my nx; i++) {
			const double time = my x1 + (i - 1) * my dx;
			my z [1] [i] = sin(2.0 * NUMpi * 377.0 * time) + NUMrandomGauss (0.0, 0.1);
		}
		bool useMultiThreading = true, extraAnalysisInfo = false;
		const int predictionOrder = 10;
		const double effectiveAnalysisWidth = 0.025, dt = 0.05, preEmphasisFrequency = 50;
		autoMelderProgress progress (U"Test multi-threading times...");
		Melder_clearInfo ();
		MelderInfo_writeLine (U"duration(s) nThread nFrames/thread toLPC(s)");
		integer numberOfThreads = maximumNumberOfThreads;
		double totalTime = 0.0;
		for (integer nThread = 1; nThread <= maximumNumberOfThreads; nThread ++) {
			const integer numberOfConcurrentThreadsToUse = nThread;
			for (integer index = 1; index <= framesPerThread.size; index ++) {
				const integer numberOfFramesPerThread = framesPerThread [index];
				MelderThread_debugMultithreading (useMultiThreading, nThread,
						numberOfFramesPerThread, extraAnalysisInfo);
				Melder_stopwatch ();
				autoLPC lpc = Sound_to_LPC_burg (me.get(), predictionOrder, effectiveAnalysisWidth, dt, preEmphasisFrequency);
				double t = Melder_stopwatch ();
				MelderInfo_writeLine (soundDuration, U" ", nThread, U" ", numberOfFramesPerThread, U" ", t);
				totalTime += t;
			}
			MelderInfo_drain ();
			try {
				Melder_progress (((double) nThread) / maximumNumberOfThreads, U"Number of threads: ", nThread);
			} catch (MelderError) {
				numberOfThreads = nThread;
				Melder_clearError ();
				break;
			}
		}
		MelderInfo_writeLine (U"Total time ", totalTime, U" seconds");
		MelderInfo_close ();
	} catch (MelderError) {
		Melder_throw (U"Could not perform timing.");
	}
}

/* End of file SampledIntoSampled.cpp */
