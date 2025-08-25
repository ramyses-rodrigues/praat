#ifndef _SoundFrameIntoSampledFrame2_h_
#define _SoundFrameIntoSampledFrame2_h_
/* SoundFrameIntoSampledFrame2.h
 *
 * Copyright (C) 2024-2025 David Weenink
 *
 * This code is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

#include "melder.h"
#include "Sound_extensions.h"
#include "SampledFrameIntoSampledFrame2.h"

inline integer getSoundFrameSize2_uneven (double approximatePhysicalAnalysisWidth, double samplingPeriod) {
	const double halfFrameDuration = 0.5 * approximatePhysicalAnalysisWidth;
	const integer halfFrameSamples = Melder_ifloor (halfFrameDuration / samplingPeriod);
	return 2 * halfFrameSamples + 1;
}

inline integer getSoundFrameSize2 (double physicalAnalysisWidth, double samplingPeriod) {
	Melder_assert (physicalAnalysisWidth > 0.0);
	Melder_assert (samplingPeriod > 0.0);
	const double numberOfSamples_real = round (physicalAnalysisWidth / samplingPeriod);
	return (integer) numberOfSamples_real;
}

inline double getPhysicalAnalysisWidth2 (double effectiveAnalysisWidth, kSound_windowShape windowShape) {
	const double physicalAnalysisWidth = ( (windowShape == kSound_windowShape::RECTANGULAR ||
		windowShape == kSound_windowShape::TRIANGULAR || windowShape == kSound_windowShape::HAMMING ||
		windowShape == kSound_windowShape::HANNING) ? effectiveAnalysisWidth : 2.0 * effectiveAnalysisWidth);
	return physicalAnalysisWidth;
}

Thing_define (SoundFrameIntoSampledFrame2, SampledFrameIntoSampledFrame2) {
	
	constSound inputSound;
	double physicalAnalysisWidth; 			// depends on the effectiveAnalysiswidth and the window window shape
	integer soundFrameSize; 				// determined by the physicalAnalysisWidth and the samplingFrequency of the Sound
	autoSound frameAsSound;
	double soundFrameExtremum;				// the largest amplitude in the inputSound frame either positive or negative
	autoVEC windowFunction;					// the actual window used of size soundFrameSize
	VEC soundFrame;
	kSound_windowShape  windowShape;		// Type: Rectangular, triangular, hamming, etc..
	bool subtractFrameMean = true;				// if true, the frame mean will be subtracted before the windowing operation
	
	void initBasicSoundFrame (constSound input, mutableSampled output, double effectiveAnalysisWidth, kSound_windowShape windowShape) {
		SoundFrameIntoSampledFrame2_Parent :: initBasic (input, output);
		inputSound = input;
		physicalAnalysisWidth = getPhysicalAnalysisWidth2 (effectiveAnalysisWidth, windowShape);
		soundFrameSize = getSoundFrameSize2 (physicalAnalysisWidth, inputSound -> dx);
	}
	
	void getInputFrame (integer currentFrame) override {
		const double midTime = Sampled_indexToX (output, currentFrame);
		integer soundFrameBegin = Sampled_xToNearestIndex (inputSound, midTime - 0.5 * physicalAnalysisWidth); // approximation
		
		for (integer isample = 1; isample <= soundFrame.size; isample ++, soundFrameBegin ++) {
			soundFrame [isample] = ( soundFrameBegin > 0 && soundFrameBegin <= inputSound -> nx ? inputSound -> z [1] [soundFrameBegin] : 0.0 );
		}
		if (subtractFrameMean)
			centre_VEC_inout (soundFrame, nullptr);
		soundFrameExtremum = NUMextremum_u (soundFrame);
		soundFrame  *=  windowFunction.get();
	}
	
	void initHeap () override {
		//SoundFrameIntoSampledFrame2_Parent :: initHeap ();
		windowFunction = raw_VEC (soundFrameSize);
		windowShape_into_VEC (windowShape, windowFunction.get());
		frameAsSound = Sound_create (1_integer, 0.0, soundFrameSize * input -> dx, soundFrameSize, input -> dx, 0.5 * input -> dx); //
		soundFrame = frameAsSound -> z.row (1);
		Melder_assert (soundFrame.size == soundFrameSize);
	}
};

#endif /* _SoundFrameIntoSampledFrame2_h_ */
 
