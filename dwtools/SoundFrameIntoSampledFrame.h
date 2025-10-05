#ifndef _SoundFrameIntoSampledFrame_h_
#define _SoundFrameIntoSampledFrame_h_
/* SoundFrameIntoSampledFrame.h
 *
 * Copyright (C) 2024-2025 David Weenink
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

#include "melder.h"
#include "NUMFourier.h"
#include "Sound_extensions.h"
#include "SampledFrameIntoSampledFrame.h"
#include "Spectrum.h"

inline integer getSoundFrameSize2_odd (double approximatePhysicalAnalysisWidth, double samplingPeriod) {
	const double halfFrameDuration = 0.5 * approximatePhysicalAnalysisWidth;
	const integer halfFrameSamples = Melder_ifloor (halfFrameDuration / samplingPeriod);
	return 2 * halfFrameSamples + 1;
}

inline integer getSoundFrameSize2 (double physicalAnalysisWidth, double samplingPeriod) {
	Melder_assert (physicalAnalysisWidth > 0.0);
	Melder_assert (samplingPeriod > 0.0);
	return Melder_iround (physicalAnalysisWidth / samplingPeriod);
}

inline double getPhysicalAnalysisWidth2 (double effectiveAnalysisWidth, kSound_windowShape windowShape) {
	const double physicalAnalysisWidth = ( windowShape == kSound_windowShape::RECTANGULAR ||
		windowShape == kSound_windowShape::TRIANGULAR || windowShape == kSound_windowShape::HAMMING ||
		windowShape == kSound_windowShape::HANNING ? effectiveAnalysisWidth : 2.0 * effectiveAnalysisWidth )
	;
	return physicalAnalysisWidth;
}

Thing_define (SoundFrameIntoSampledFrame, SampledFrameIntoSampledFrame) {
	
	constSound inputSound;
	double physicalAnalysisWidth; 			// depends on the effectiveAnalysiswidth and the window window shape
	integer soundFrameSize; 				// determined by the physicalAnalysisWidth and the samplingFrequency of the Sound
	autoSound frameAsSound;
	double soundFrameExtremum;				// the largest amplitude in the inputSound frame either positive or negative
	autoVEC windowFunction;					// the actual window used of size soundFrameSize
	VEC soundFrame;
	kSound_windowShape windowShape;			// Type: Rectangular, triangular, hamming, etc..
	bool subtractFrameMean = true;			// if true, the frame mean will be subtracted before the windowing operation
	bool wantSpectrum = false;				// the spectrum of the frameAsSound;
	autoSpectrum spectrum;
	integer fftInterpolationFactor = 1;		// 0 = DFT, 1 = FFT, 2, 4, 8 FFT with extra zero's
	integer numberOfFourierSamples;
	autoVEC fourierSamples;					// size = numberOfFourierSamples
	autoNUMFourierTable fourierTable;		// of dimension numberOfFourierSamples;
	
	void initBasicSoundFrameIntoSampledFrame (constSound input, mutableSampled output, double effectiveAnalysisWidth, 
		kSound_windowShape windowShape)
	{
		SoundFrameIntoSampledFrame_Parent :: initBasic (input, output);
		inputSound = input;
		our windowShape = windowShape;
		physicalAnalysisWidth = getPhysicalAnalysisWidth2 (effectiveAnalysisWidth, windowShape);
	}
		
	void copyBasic (constSampledFrameIntoSampledFrame other2) override {
		constSoundFrameIntoSampledFrame other = reinterpret_cast<constSoundFrameIntoSampledFrame> (other2);
		SoundFrameIntoSampledFrame_Parent :: copyBasic (other);
		our inputSound = other -> inputSound;
		our physicalAnalysisWidth = other -> physicalAnalysisWidth;
		our windowShape = other -> windowShape;
		our subtractFrameMean = other -> subtractFrameMean;
		our wantSpectrum = other -> wantSpectrum;
		our fftInterpolationFactor = other -> fftInterpolationFactor;
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
		if (wantSpectrum)
			soundFrameIntoSpectrum ();
	}
	
	void initHeap () override {
		SoundFrameIntoSampledFrame_Parent :: initHeap ();
		soundFrameSize = getSoundFrameSize2 (physicalAnalysisWidth, inputSound -> dx);
		windowFunction = raw_VEC (soundFrameSize);   // TODO: move out of thread repetition
		windowShape_into_VEC (windowShape, windowFunction.get());   // TODO: move out of thread repetition
		frameAsSound = Sound_create (1_integer, 0.0, soundFrameSize * input -> dx, soundFrameSize, input -> dx, 0.5 * input -> dx); //
		soundFrame = frameAsSound -> z.row (1);
		Melder_assert (soundFrame.size == soundFrameSize);
		if (wantSpectrum) {
			numberOfFourierSamples = frameAsSound -> nx;
			if (fftInterpolationFactor > 0) {
				numberOfFourierSamples = Melder_iroundUpToPowerOfTwo (numberOfFourierSamples);
				for (integer imultiply = fftInterpolationFactor; imultiply > 1; imultiply --)
					numberOfFourierSamples *= 2;
			}
			fourierSamples = raw_VEC (numberOfFourierSamples);
			const integer numberOfFrequencies = numberOfFourierSamples / 2 + 1;
			fourierTable = NUMFourierTable_create (numberOfFourierSamples);
			spectrum = Spectrum_create (0.5 / frameAsSound -> dx, numberOfFrequencies);
			spectrum -> dx = 1.0 / (frameAsSound -> dx * numberOfFourierSamples);
		}
	}
	
	void soundFrameToForwardFourierTransform () {
		const integer numberOfChannels = frameAsSound -> ny;
		if (numberOfChannels == 1)
			fourierSamples.part (1, soundFrameSize)  <<=  frameAsSound -> z.row (1);
		else {
			/*
				Multiple channels: take the average.
				*/
			for (integer ichan = 1; ichan <= numberOfChannels; ichan ++)
				fourierSamples.part (1, soundFrameSize)  +=  frameAsSound -> z.row (ichan);
			fourierSamples.part (1, soundFrameSize)  *=  1.0 / numberOfChannels;
		}
		fourierSamples.part (soundFrameSize + 1, numberOfFourierSamples)  <<=  0.0;
		NUMfft_forward (fourierTable.get(), fourierSamples.get());
	}
	
	void soundFrameIntoSpectrum () {
		
		soundFrameToForwardFourierTransform ();
		
		const VEC re = spectrum -> z.row (1);
		const VEC im = spectrum -> z.row (2);
		const integer numberOfFrequencies = spectrum -> nx;
		const double scaling = output -> dx;
		re [1] = fourierSamples [1] * scaling;
		im [1] = 0.0;
		for (integer i = 2; i < numberOfFrequencies; i ++) {
			re [i] = fourierSamples [i + i - 2] * scaling;   // fourierSamples [2], [4], ...
			im [i] = fourierSamples [i + i - 1] * scaling;   // fourierSamples [3], [5], ...
		}
		if ((numberOfFourierSamples & 1) != 0) {
			if (numberOfFourierSamples > 1) {
				re [numberOfFrequencies] = fourierSamples [numberOfFourierSamples - 1] * scaling;
				im [numberOfFrequencies] = fourierSamples [numberOfFourierSamples] * scaling;
			}
		} else {
			re [numberOfFrequencies] = fourierSamples [numberOfFourierSamples] * scaling;
			im [numberOfFrequencies] = 0.0;
		}
	}
};

#endif /* _SoundFrameIntoSampledFrame_h_ */
 
