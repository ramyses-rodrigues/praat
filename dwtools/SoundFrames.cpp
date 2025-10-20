/* SoundFrames.cpp
 *
 * Copyright (C) 2025 David Weenink
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

#include "SoundFrames.h"

Thing_implement (SoundFrames, Thing, 0);

void structSoundFrames :: init (constSound input, double effectiveAnalysisWidth, double timeStep, 
	kSound_windowShape windowShape, bool subtractFrameMean, bool wantSpectrum, 
	integer fftInterpolationFactor)
{
	our inputSound = input;
	our physicalAnalysisWidth = getPhysicalAnalysisWidth (effectiveAnalysisWidth, windowShape);
	if (timeStep == 0.0) {
		// calculate output_dt
	}
	our dt = timeStep;
	Sampled_shortTermAnalysis (inputSound, physicalAnalysisWidth, dt, & our numberOfFrames, & our t1);
	initCommon (windowShape, subtractFrameMean, wantSpectrum, fftInterpolationFactor);
}
	
void structSoundFrames :: initWithSampled (constSound input, constSampled output, double effectiveAnalysisWidth,
	kSound_windowShape windowShape, bool subtractFrameMean, bool wantSpectrum,
	integer fftInterpolationFactor)
{
	Melder_require (input -> xmin == output -> xmin && input -> xmax == output -> xmax,
		U"The domains of Sound ", input, U" and Sampled ", output , U" should be equal.");
	our inputSound = input;
	our t1 = output -> x1;
	our numberOfFrames = output -> nx;
	our dt = output -> dx;
	our physicalAnalysisWidth = getPhysicalAnalysisWidth (effectiveAnalysisWidth, windowShape);	
	initCommon (windowShape, subtractFrameMean, wantSpectrum, fftInterpolationFactor);
}

void structSoundFrames :: initCommon (kSound_windowShape windowShape, 
	bool subtractFrameMean, bool wantSpectrum, integer fftInterpolationFactor)
{
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

VEC structSoundFrames :: getFrame (integer iframe) {
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
	return soundFrame;
}

void structSoundFrames :: soundFrameToForwardFourierTransform () {
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

void structSoundFrames :: soundFrameIntoSpectrum () {
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

autoSoundFrames SoundFrames_create (constSound input, constSampled output, double effectiveAnalysisWidth,
	kSound_windowShape windowShape, bool subtractFrameMean, bool wantSpectrum, integer fftInterpolationFactor)
{
	try {
		autoSoundFrames me = Thing_new (SoundFrames);
		my initWithSampled (input, output, effectiveAnalysisWidth, windowShape, subtractFrameMean,
			wantSpectrum, fftInterpolationFactor);
		return me;
	} catch (MelderError) {
		Melder_throw (U"SoundFrames (with Sampled) could not be created.");
	}
}


autoSoundFrames SoundFrames_create (constSound input, double effectiveAnalysisWidth,
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


/* End of file SoundFrames.cpp */
