#ifndef _SoundFrames_h_
#define _SoundFrames_h_
/* SoundFrames.h
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

inline integer getSoundFrameSize_odd (double approximatePhysicalAnalysisWidth, double samplingPeriod) {
	const double halfFrameDuration = 0.5 * approximatePhysicalAnalysisWidth;
	const integer halfFrameSamples = Melder_ifloor (halfFrameDuration / samplingPeriod);
	return 2 * halfFrameSamples + 1;
}

inline integer getSoundFrameSize (double physicalAnalysisWidth, double samplingPeriod) {
	Melder_assert (physicalAnalysisWidth > 0.0);
	Melder_assert (samplingPeriod > 0.0);
	return Melder_iround (physicalAnalysisWidth / samplingPeriod);
}

inline double getPhysicalAnalysisWidth (double effectiveAnalysisWidth, kSound_windowShape windowShape) {
	const double physicalAnalysisWidth = ( windowShape == kSound_windowShape::RECTANGULAR ||
		windowShape == kSound_windowShape::TRIANGULAR || windowShape == kSound_windowShape::HAMMING ||
		windowShape == kSound_windowShape::HANNING ? effectiveAnalysisWidth : 2.0 * effectiveAnalysisWidth )
	;
	return physicalAnalysisWidth;
}

Thing_define (SoundFrames, Thing) {

	constSound inputSound;
	double t1, dt;
	integer numberOfFrames;
	double physicalAnalysisWidth;			// depends on the effectiveAnalysiswidth and the window window shape
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

	void init (constSound input, double effectiveAnalysisWidth, double timeStep, 
		kSound_windowShape windowShape, bool subtractFrameMean, bool wantSpectrum, 
		integer fftInterpolationFactor)
	{
		our physicalAnalysisWidth = getPhysicalAnalysisWidth (effectiveAnalysisWidth, windowShape);
		if (timeStep == 0.0) {
			// calculate output_dt
		}
		our dt = timeStep;
		Sampled_shortTermAnalysis (inputSound, physicalAnalysisWidth, dt, & numberOfFrames, & t1);
		initCommon (input, windowShape, subtractFrameMean, wantSpectrum, fftInterpolationFactor);
	}
	
	void initWithSampled (constSound input, Sampled output, double effectiveAnalysisWidth,
		kSound_windowShape windowShape, bool subtractFrameMean, bool wantSpectrum,
		integer fftInterpolationFactor)
	{
		Melder_require (input -> xmin == output -> xmin && input -> xmax == output -> xmax,
			U"The domains of Sound ", input, U" and Sampled ", output , U" should be equal.");
		our t1 = output -> x1;
		our numberOfFrames = output -> nx;
		our dt = output -> dx;
		// check
		const double physicalAnalysisWidth2 = getPhysicalAnalysisWidth (effectiveAnalysisWidth, windowShape);
		integer numberOfFrames_2;
		double t1_2;
		Sampled_shortTermAnalysis (inputSound, physicalAnalysisWidth, dt, & numberOfFrames_2, & t1_2);
		Melder_assert (numberOfFrames_2 == numberOfFrames);
		Melder_assert (t1_2 == t1);
		initCommon (input, windowShape, subtractFrameMean, wantSpectrum, fftInterpolationFactor);
	}
	
	void initCommon (constSound input, kSound_windowShape windowShape, bool subtractFrameMean, bool wantSpectrum, 
		integer fftInterpolationFactor)
	{
		inputSound = input;
		our windowShape = windowShape;
		our subtractFrameMean = subtractFrameMean;
		our wantSpectrum = wantSpectrum;

		soundFrameSize = getSoundFrameSize (physicalAnalysisWidth, inputSound -> dx);
		windowFunction = raw_VEC (soundFrameSize);   // TODO: move out of thread repetition
		windowShape_into_VEC (windowShape, windowFunction.get());
		frameAsSound = Sound_create (1_integer, 0.0, soundFrameSize * inputSound -> dx, soundFrameSize,
			inputSound -> dx, 0.5 * inputSound -> dx);
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
	
	void getFrame (integer iframe) {
		const double midTime = t1 + (iframe - 1) * dt;
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

	void soundFrameToForwardFourierTransform () {
		Melder_assert (wantSpectrum);
		const integer numberOfChannels = frameAsSound -> ny;

		fourierSamples.part (1, soundFrameSize)  <<=  frameAsSound -> z.row (1);
		if (numberOfChannels > 1) {
			/*
				Multiple channels: take the average.
			*/
			for (integer ichan = 2; ichan <= numberOfChannels; ichan ++)
				fourierSamples.part (1, soundFrameSize)  +=  frameAsSound -> z.row (ichan);
			fourierSamples.part (1, soundFrameSize)  *=  1.0 / numberOfChannels;
		}
		fourierSamples.part (soundFrameSize + 1, numberOfFourierSamples)  <<=  0.0;
		NUMfft_forward (fourierTable.get(), fourierSamples.get());
	}

	void soundFrameIntoSpectrum () {
		Melder_assert (wantSpectrum);
		soundFrameToForwardFourierTransform ();

		const VEC re = spectrum -> z.row (1);
		const VEC im = spectrum -> z.row (2);
		const integer numberOfFrequencies = spectrum -> nx;
		const double scaling = frameAsSound -> dx;
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

autoSoundFrames SoundFrame_create (constSound input, double effectiveAnalysisWidth,
	double timeStep, kSound_windowShape windowShape, bool subtractFrameMean, bool wantSpectrum, 
	integer fftInterpolationFactor)
{
	try {
		autoSoundFrames me = Thing_new (SoundFrames);
		my init (input, effectiveAnalysisWidth, timeStep, windowShape, subtractFrameMean,
			wantSpectrum, fftInterpolationFactor);
		return me;
	} catch (MelderError) {
		Melder_throw (U"SoundFrames could not be created.");
	}
}

#endif /* _SoundFrames_h_ */
