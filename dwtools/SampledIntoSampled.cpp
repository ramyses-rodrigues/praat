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

integer SampledIntoSampled_analyseThreaded (mutableSampledIntoSampled me)
{
	try {
		SampledFrameIntoSampledFrame frameIntoFrame = my frameIntoFrame.get();

		const integer numberOfFrames = my output -> nx;
		
		std::atomic<integer> globalFrameErrorCount (0);
		
		if (MelderThread_getUseMultithreading ()) {
			integer numberOfThreadsNeeded, numberOfFramesPerThread;
			MelderThread_getInfo (numberOfFrames, & numberOfThreadsNeeded, & numberOfFramesPerThread);

			/*
				We need to reserve all the working memory for each thread beforehand.
			*/
			const integer numberOfThreadsToUse = Melder_clippedRight (MelderThread_getMaximumNumberOfConcurrentThreads (),
				NUMrandom_maximumNumberOfParallelThreads);
			const integer numberOfThreads = std::min (numberOfThreadsToUse, numberOfThreadsNeeded);
			frameIntoFrame -> maximumNumberOfFrames = numberOfFramesPerThread;
			OrderedOf<structSampledFrameIntoSampledFrame> workThreads;
			for (integer ithread = 1; ithread <= numberOfThreads; ithread ++) {
				autoSampledFrameIntoSampledFrame frameIntoFrameCopy = Data_copy (frameIntoFrame);
				workThreads. addItem_move (frameIntoFrameCopy.move());
			}

			/*
				The following cannot be an `autovector`, because autovectors don't destroy their elements.
				So it has to be std::vector.
				Also, the default initialization of a std::thread my not be guaranteed to be all zeroes.
			*/
			std::vector<std::thread> threads (1+numberOfThreads);
			integer numberOfThreadsInRun;
			try {
				const integer numberOfThreadRuns = Melder_iroundUp ((double) numberOfThreadsNeeded / numberOfThreads);
				//TRACE
				trace (numberOfThreadRuns, U" ", numberOfFrames, U" ", numberOfThreadsNeeded, U" ", numberOfThreads);
				const integer numberOfFramesInRun = numberOfThreads * numberOfFramesPerThread;
				const integer remainingThreads = numberOfThreadsNeeded % numberOfThreads;
				const integer numberOfThreadsInLastRun = ( remainingThreads == 0 ? numberOfThreads : remainingThreads );
				for (integer irun = 1; irun <= numberOfThreadRuns; irun ++) {
					numberOfThreadsInRun = ( irun < numberOfThreadRuns ? numberOfThreads : numberOfThreadsInLastRun );
					const integer lastFrameInRun = ( irun < numberOfThreadRuns ? numberOfFramesInRun * irun : numberOfFrames);
					for (integer ithread = 1; ithread <= numberOfThreadsInRun; ithread ++) {
						SampledFrameIntoSampledFrame frameIntoFrameCopy = workThreads.at [ithread];
						const integer startFrame = numberOfFramesInRun * (irun - 1) + 1 + (ithread - 1) * numberOfFramesPerThread;
						const integer endFrame = ( ithread == numberOfThreadsInRun ? lastFrameInRun : startFrame + numberOfFramesPerThread - 1 );
						frameIntoFrameCopy -> startFrame = startFrame;
						frameIntoFrameCopy -> currentNumberOfFrames = endFrame - startFrame + 1;
						
						auto analyseFrames = [&globalFrameErrorCount] (int threadNumber, SampledFrameIntoSampledFrame fifthread, integer fromFrame, integer toFrame) {
							NUMrandom_setChannel (threadNumber);
							fifthread -> inputFramesToOutputFrames (fromFrame, toFrame);
							globalFrameErrorCount += fifthread -> framesErrorCount;
						};

						threads [ithread] = std::thread (analyseFrames, ithread - 1, frameIntoFrameCopy, startFrame, endFrame);
					}
					for (integer ithread = 1; ithread <= numberOfThreadsInRun; ithread ++) {
						threads [ithread]. join ();
						SampledFrameIntoSampledFrame frameIntoFrameCopy = workThreads.at [ithread];
					}
				}
			} catch (MelderError) {
				for (integer ithread = 1; ithread <= numberOfThreadsInRun; ithread ++)
					if (threads [ithread]. joinable ())
						threads [ithread]. join ();
				Melder_clearError ();
				throw;
			}
			my globalFrameErrorCount = globalFrameErrorCount;
		} else {
			frameIntoFrame -> inputFramesToOutputFrames (1, numberOfFrames); // no threading
			globalFrameErrorCount = frameIntoFrame -> framesErrorCount;
		}
		if (frameIntoFrame -> updateStatus)
			my status -> showStatus ();
		return globalFrameErrorCount;
	} catch (MelderError) {
		Melder_throw (me, U"The Sampled analysis could not be done.");
	}
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
						numberOfFramesPerThread, numberOfFramesPerThread, extraAnalysisInfo);
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
