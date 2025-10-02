/* Sound_to_Formant_mt.cpp
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this work. If not, see <http://www.gnu.org/licenses/>.
 */

#include <thread>
#include "Sound_to_Formant_mt.h"
#include "Sound_extensions.h"
#include "LPC_and_Formant.h"
#include "SampledIntoSampled.h"

Thing_implement (SoundFrameIntoFormantFrame, SampledFrameIntoSampledFrame, 0);

void structSoundFrameIntoFormantFrame :: initBasicSoundFrameIntoFormantFrame (constSound inputSound, mutableFormant outputFormant) {
	SoundFrameIntoFormantFrame_Parent :: initBasic (inputSound, outputFormant);
	
}

void structSoundFrameIntoFormantFrame :: initHeap () {
	SoundFrameIntoFormantFrame_Parent :: initHeap ();
	soundIntoLPC -> initHeap ();
	lpcIntoFormant -> initHeap ();
}

bool structSoundFrameIntoFormantFrame :: inputFrameIntoOutputFrame (integer iframe) {
	bool result = soundIntoLPC -> inputFrameIntoOutputFrame (iframe);
	if (result)
		result = lpcIntoFormant -> inputFrameIntoOutputFrame (iframe);
	return result;
}

autoSoundFrameIntoFormantFrame SoundFrameIntoFormantFrame_create (SoundFrameIntoLPCFrame soundIntoLPC, LPCFrameIntoFormantFrame lpcIntoFormant) {
	try {
		autoSoundFrameIntoFormantFrame me = Thing_new (SoundFrameIntoFormantFrame);
		my soundIntoLPC.adoptFromAmbiguousOwner (soundIntoLPC);
		my lpcIntoFormant.adoptFromAmbiguousOwner (lpcIntoFormant);
		return me;
	} catch (MelderError) {
		Melder_throw (U"Cannot create SoundFrameIntoFormantFrame.");
	}
}

/*
	Precondition:
		Sound already has the 'right' sampling frequency and has been pre-emphasized
*/


static void Sound_to_Formant_common (constSound inputSound, double& dt, double numberOfFormants_real, double maximumFrequency,
	double effectiveAnalysisWidth, double preEmphasisFrequency,	double safetyMargin, 
	autoFormant& outputFormant, autoLPC& outputLPC, autoSound& sound)
{
	try {
		if (dt<= 0.0)
			dt = effectiveAnalysisWidth / 4.0;
		autoSound resampled = Sound_resampleAndOrPreemphasize (inputSound, maximumFrequency, 50, preEmphasisFrequency);
		const integer numberOfPoles = 2.0 * numberOfFormants_real;
		const double physicalAnalysisWidth = getPhysicalAnalysisWidth2 (effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
		integer numberOfFrames;
		double t1;
		Sampled_shortTermAnalysis (sound.get(), physicalAnalysisWidth, dt, & numberOfFrames, & t1);
		const integer numberOfFormants = ( safetyMargin == 0.0 ? numberOfPoles : (numberOfPoles + 1) / 2 );
		autoFormant formant = Formant_create (sound -> xmin, sound -> xmax, numberOfFrames, dt, t1, numberOfFormants);
		for (integer iframe = 1; iframe <= numberOfFrames; iframe ++) {
			Formant_Frame formantFrame = &formant -> frames [iframe];
			Formant_Frame_init (formantFrame, numberOfFormants);
		}
		outputFormant = formant.move();
		autoLPC lpc = LPC_createCompletelyInitialized (sound -> xmin, sound -> xmax, formant -> nx, formant -> dx,
			formant -> x1, numberOfPoles, sound -> dx);
		outputLPC = lpc.move();
		sound = resampled.move();
	} catch (MelderError) {
		Melder_throw (U"Cannot create Formant or LPC or Sound.");
	}
	
}
static autoFormant createFormant_common (constSound me, double dt, integer numberOfPoles, double effectiveAnalysisWidth,
	double safetyMargin)
{
	integer numberOfFrames;
	double t1;
	const double physicalAnalysisWidth = getPhysicalAnalysisWidth2 (effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
	Sampled_shortTermAnalysis (me, physicalAnalysisWidth, dt, & numberOfFrames, & t1);
	const integer numberOfFormants = numberOfFormantsFromNumberOfCoefficients2 (numberOfPoles, safetyMargin);
	autoFormant formant = Formant_create (my xmin, my xmax, numberOfFrames, dt, t1, numberOfFormants);
	return formant;
}

autoFormant Sound_to_Formant_burg_mt (constSound me, double dt, double numberOfFormants, double maximumFrequency,
	double effectiveAnalysisWidth, double preEmphasisFrequency, double safetyMargin)
{
	try {
		autoSound sound;
		autoFormant formant;
		autoLPC lpc;
		Sound_to_Formant_common (me, dt, numberOfFormants, maximumFrequency, effectiveAnalysisWidth, preEmphasisFrequency,safetyMargin,
			formant, lpc, sound);
		autoSoundFrameIntoLPCFrameBurg soundIntoLPC = SoundFrameIntoLPCFrameBurg_create (me, lpc.get(), effectiveAnalysisWidth,
				kSound_windowShape ::GAUSSIAN_2);
		autoLPCFrameIntoFormantFrame lpcIntoFormant = LPCFrameIntoFormantFrame_create (lpc.get(), formant.get(), safetyMargin);
		autoSoundFrameIntoFormantFrame soundIntoFormant = SoundFrameIntoFormantFrame_create (soundIntoLPC.releaseToAmbiguousOwner(),
			lpcIntoFormant.releaseToAmbiguousOwner());
		SampledIntoSampled_mt (soundIntoFormant.get(), 40);
		return formant;
	} catch (MelderError) {
		Melder_throw (U"Could not create Formant (burg).");
	}
}

autoFormant Sound_and_LPC_to_Formant (constSound me, constLPC lpc, double effectiveAnalysisWidth, double preEmphasisFrequency, 
	double safetyMargin, double k_stdev, integer itermax, double tol, double location, bool wantlocation)
{
	try {
		const double maximumFrequency = 1.0 / lpc -> samplingPeriod;
		autoSound sound = Sound_resampleAndOrPreemphasize (me, maximumFrequency, 50, preEmphasisFrequency);
		autoFormant formant = createFormant_common (sound.get(), lpc -> dx, lpc -> maxnCoefficients, effectiveAnalysisWidth, safetyMargin);
		autoLPC outputLPC = Data_copy (lpc);
		autoLPCFrameAndSoundFrameIntoLPCFrameRobust lpcFrameAndSoundFrameIntoLPCFrame = Thing_new (LPCFrameAndSoundFrameIntoLPCFrameRobust);
		lpcFrameAndSoundFrameIntoLPCFrame -> initBasicLPCFrameAndSoundFrameIntoLPCFrameRobust (lpc, sound.get(),
			outputLPC.get(), effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2, k_stdev, itermax, tol, wantlocation);
		autoLPCFrameIntoFormantFrame lpcIntoFormant = LPCFrameIntoFormantFrame_create (outputLPC.get(), formant.get(), safetyMargin);
		autoSoundFrameIntoFormantFrame soundIntoFormant = SoundFrameIntoFormantFrame_create (lpcFrameAndSoundFrameIntoLPCFrame.releaseToAmbiguousOwner(), lpcIntoFormant.releaseToAmbiguousOwner());
		SampledIntoSampled_mt (soundIntoFormant.get(), 40);
		return formant;
	} catch (MelderError) {
		Melder_throw (U"Could not create Formant from Sound and LPC.");
	}
}

autoFormant Sound_to_Formant_robust_mt (constSound inputSound, double dt, double numberOfFormants, double maximumFrequency,
	double effectiveAnalysisWidth, double preEmphasisFrequency, double safetyMargin, double k_stdev, integer itermax, double tol,
	double location, bool wantlocation)
{
	try {
		autoSound sound;
		autoFormant outputFormant;
		autoLPC outputLPC;		
		Sound_to_Formant_common (inputSound, dt, numberOfFormants, maximumFrequency, effectiveAnalysisWidth, preEmphasisFrequency,
			safetyMargin, outputFormant, outputLPC, sound);
		const kSound_windowShape windowShape = kSound_windowShape::GAUSSIAN_2;
		
		Sound_into_LPC_robust (sound.get(), outputLPC.get(), effectiveAnalysisWidth, k_stdev, itermax, tol, wantlocation);
		LPC_into_Formant (outputLPC.get(), outputFormant.get(), safetyMargin);
		return outputFormant;
	} catch (MelderError) {
		Melder_throw (inputSound, U": no robust Formant created.");
	}
}

autoFormant Sound_to_Formant_robust (Sound me, double dt_in, double numberOfFormants, double maximumFrequency,
	double effectiveAnalysisWidth, double preEmphasisFrequency, double safetyMargin,
	double numberOfStandardDeviations, integer maximumNumberOfIterations, double tolerance,
	double location, bool wantlocation)
{
	const double dt = dt_in > 0.0 ? dt_in : effectiveAnalysisWidth / 4.0;
	const double nyquist = 0.5 / my dx;
	const integer predictionOrder = Melder_ifloor (2 * numberOfFormants);
	try {
		autoSound sound;
		if (maximumFrequency <= 0.0 || fabs (maximumFrequency / nyquist - 1.0) < 1.0e-12)
			sound = Data_copy (me);   // will be modified
		else
			sound = Sound_resample (me, maximumFrequency * 2.0, 50);

		autoLPC lpc = Sound_to_LPC_auto (sound.get(), predictionOrder, effectiveAnalysisWidth, dt, preEmphasisFrequency);
		autoLPC lpcRobust = LPC_and_Sound_to_LPC_robust (lpc.get(), sound.get(), effectiveAnalysisWidth, preEmphasisFrequency,
			numberOfStandardDeviations, maximumNumberOfIterations, tolerance, wantlocation);
		autoFormant thee = LPC_to_Formant (lpcRobust.get(), safetyMargin);
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": no robust Formant created.");
	}
}


/*
void Sound_into_Formant_robust_mt (constSound me, mutableFormant thee, double effectiveAnalysisWidth, integer numberOfPoles, double safetyMargin,
	double k_stdev, integer itermax, double tol, double location, bool wantlocation)
{
	autoSoundFrameIntoFormantFrameRobust ws = SoundFrameIntoFormantFrameRobust_create (me, thee,
		effectiveAnalysisWidth, kSound_windowShape :: GAUSSIAN_2, k_stdev, itermax, tol, location, wantlocation, numberOfPoles, safetyMargin);
	
	
	SampledToSampled_analyseThreaded (ws.get());
}

autoFormant Sound_to_Formant_robust_mt (constSound me, double dt_in, double numberOfFormants, double maximumFrequency,
	double effectiveAnalysisWidth, double preEmphasisFrequency, double safetyMargin, double k_stdev, integer itermax, double tol,
	double location, bool wantlocation)
{
	const double dt = dt_in > 0.0 ? dt_in : effectiveAnalysisWidth / 4.0;
	try {
		autoSound sound = Sound_resampleAndOrPreemphasize (me, maximumFrequency, 50, preEmphasisFrequency);
		integer numberOfFrames;
		double t1;
		const double physicalAnalysisWidth = getPhysicalAnalysisWidth (effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
		Sampled_shortTermAnalysis (sound.get(), physicalAnalysisWidth, dt, & numberOfFrames, & t1);
		const integer numberOfPoles = numberOfPolesFromNumberOfFormants (numberOfFormants);
		const integer numberOfFormants = numberOfFormantsFromNumberOfCoefficients (numberOfPoles, safetyMargin);
		autoFormant formant = Formant_create (my xmin, my xmax, numberOfFrames, dt, t1, numberOfFormants);
		Sound_into_Formant_robust_mt (sound.get(), formant.get(), effectiveAnalysisWidth, numberOfPoles, safetyMargin, k_stdev,
			itermax, tol, location, wantlocation);
		return formant;
	} catch (MelderError) {
		Melder_throw (me, U": no robust Formant created.");
	}
}
*/
/* End of file Sound_to_Formant_mt.cpp */
