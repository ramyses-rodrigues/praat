/* Sound_and_LPC2.cpp
 *
 * Copyright (C) 1994-2025 David Weenink
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

/*
 djmw 20020625 GPL header
*/

#include "SampledIntoSampled2.h"
#include "Sound_and_LPC2.h"
#include "Sound_extensions.h"
#include "Spectrum.h"
#include "NUM2.h"

Thing_implement (SoundFrameIntoLPCFrame2, SoundFrameIntoSampledFrame2, 0);

void structSoundFrameIntoLPCFrame2 :: initBasicSoundFrameAndLPC (constSound input, mutableLPC outLPC, 
	double effectiveAnalysisWidth, kSound_windowShape windowShape)
{
	initBasicSoundFrame (input, outLPC, effectiveAnalysisWidth, windowShape);
	outputLPC = outLPC;
	order = outLPC -> maxnCoefficients;
	orderp1 = order + 1;
}

void structSoundFrameIntoLPCFrame2 :: initHeap () {
	SoundFrameIntoLPCFrame2_Parent :: initHeap ();
	a = zero_VEC (orderp1);
}

/*********************** Autocorrelation method *************************************************************/

Thing_implement (SoundFrameIntoLPCFrame2Auto, SoundFrameIntoLPCFrame2, 0);

void structSoundFrameIntoLPCFrame2Auto :: initHeap () {
	SoundFrameIntoLPCFrame2Auto_Parent :: initHeap ();
	r = raw_VEC (orderp1);
	rc = raw_VEC (orderp1);
}
		
bool structSoundFrameIntoLPCFrame2Auto :: inputFrameIntoOutputFrame (integer currentFrame) {
	LPC_Frame lpcFrame = & outputLPC -> d_frames [currentFrame];
	Melder_assert (lpcFrame -> nCoefficients > 0);
	integer frameAnalysisInfo = 0;

	VEC x = soundFrame;
	
	/*
		Compute the autocorrelations
	*/
	lpcFrame -> a.get()  <<=  0.0;
	lpcFrame -> gain = 0.0;
	for (integer i = 1; i <= orderp1; i ++)
		r [i] = NUMinner (x.part (1, x.size - i + 1), x.part (i, x.size));
	if (r [1] == 0.0) {
		/*
			The sound frame contains only zero's
		*/
		lpcFrame -> nCoefficients = 0;
		lpcFrame -> a.resize (lpcFrame -> nCoefficients); // maintain invariant
		frameAnalysisInfo = 1;
		return false;
	}
	a [1] = 1.0;
	a [2] = rc [1] = - r [2] / r [1];
	double gain = r [1] + r [2] * rc [1];
	lpcFrame -> gain = gain;
	integer i = 1;
	for (i = 2; i <= order; i ++) {
		longdouble s = 0.0;
		for (integer j = 1; j <= i; j ++)
			s += r [i - j + 2] * a [j];
		rc [i] = - s / gain;
		for (integer j = 2; j <= i / 2 + 1; j ++) {
			const double at = a [j] + rc [i] * a [i - j + 2];
			a [i - j + 2] += rc [i] * a [j];
			a [j] = at;
		}
		a [i + 1] = rc [i];
		gain += rc [i] * s;
		if (gain <= 0.0) {
			frameAnalysisInfo = 2;
			break;
		}
		lpcFrame -> gain = gain;
	}
	-- i;
	lpcFrame -> a.part (1, i)  <<=  a.part (2, i + 1);
	lpcFrame -> a.resize (i);
	lpcFrame -> nCoefficients = lpcFrame -> a.size; // maintain invariant
	return true;
}

autoSoundFrameIntoLPCFrame2Auto SoundFrameIntoLPCFrame2Auto_create (constSound input, mutableLPC outLPC, double effectiveAnalysisWidth, kSound_windowShape windowShape) {
	try {
		autoSoundFrameIntoLPCFrame2Auto me = Thing_new (SoundFrameIntoLPCFrame2Auto);
		my initBasicSoundFrameAndLPC (input, outLPC, effectiveAnalysisWidth, windowShape);
		return me;
	} catch (MelderError) {
		Melder_throw (U"Cannot create SoundFrameIntoLPCFrame2Auto");
	}
}

/*********************** Covariance method *************************************************************/

Thing_implement (SoundFrameIntoLPCFrame2Covar, SoundFrameIntoLPCFrame2, 0);

void structSoundFrameIntoLPCFrame2Covar :: initBasicSoundFrameAndLPC (constSound input, mutableLPC outLPC, double effectiveAnalysisWidth, kSound_windowShape windowShape) 
{
	SoundFrameIntoLPCFrame2Covar_Parent :: initBasicSoundFrameAndLPC (input, outLPC, effectiveAnalysisWidth, windowShape);
	order2 = order * (order + 1);
}

void structSoundFrameIntoLPCFrame2Covar :: initHeap () {
	SoundFrameIntoLPCFrame2_Parent :: initHeap ();
	b = raw_VEC (order2);
	grc = raw_VEC (order);
	beta = raw_VEC (order);
	cc = raw_VEC (orderp1);
}

bool structSoundFrameIntoLPCFrame2Covar :: inputFrameIntoOutputFrame (integer currentFrame) {
	LPC_Frame lpcFrame = & outputLPC -> d_frames [currentFrame];
	const integer n = soundFrameSize, m = order;
	integer frameAnalysisInfo = 0;
	
	if (lpcFrame -> nCoefficients == 0) {
		frameAnalysisInfo = 6;
		return false;
	}		
	constVEC x = soundFrame;
	/*
		Compute the covariances
	*/
	constVECVU xi = x.part (m + 1, n), xim1 = x.part (m, n - 1);
	double gain = NUMinner (xi, xi);
	cc [1] = NUMinner (xi, xim1);
	cc [2] = NUMinner (xim1, xim1);

	if (gain == 0.0) {
		frameAnalysisInfo = 1;
		lpcFrame -> nCoefficients = 0;
		lpcFrame -> gain = gain;
		lpcFrame -> a.resize (lpcFrame -> nCoefficients); //maintain invariant
		return false;
	}

	b [1] = 1.0;
	for (integer i = 2; i <= b.size; i ++)
		b [i] = 0.0;
	beta [1] = cc [2];
	a [1] = 1.0;
	a [2] = grc [1] = -cc [1] / cc [2];
	lpcFrame -> gain = gain += grc [1] * cc [1];
	integer i = 2;
	for (i = 2; i <= m; i ++) { // 130
		for (integer j = 1; j <= i; j ++)
			cc [i - j + 2] = cc [i - j + 1] + x [m - i + 1] * x [m - i + j] - x [n - i + 1] * x [n - i + j];
		
		// cc[1]=0.0; for (integer j = m + 1; j <= n; j ++) cc [1] += x [j - i] * x [j];
		cc [1] = NUMinner (x.part (m + 1 - i, n - i), x.part (m + 1, n)); //30
			
		b [i * (i + 1) / 2] = 1.0;
		for (integer j = 1; j <= i - 1; j ++) { // 70
			if (beta [j] < 0.0) {
				frameAnalysisInfo = 2;
				goto end;
			} else if (beta [j] == 0.0)
				continue;

			double gam = 0.0;
			for (integer k = 1; k <= j; k ++)
				gam += cc [k + 1] * b [j * (j - 1) / 2 + k]; // 50
			gam /= beta [j];
			for (integer k = 1; k <= j; k ++)
				b [i * (i - 1) / 2 + k] -= gam * b [j * (j - 1) / 2 + k]; // 60
		}

		beta [i] = 0.0;
		for (integer j = 1; j <= i; j ++)
			beta [i] += cc [j + 1] * b [i * (i - 1) / 2 + j]; // 80
		if (beta [i] <= 0.0) {
			frameAnalysisInfo = 3;
			break;
		}
		double s = 0.0;
		for (integer j = 1; j <= i; j ++)
			s += cc [j] * a [j]; // 100
		grc [i] = -s / beta [i];

		for (integer j = 2; j <= i; j ++)
			a [j] += grc [i] * b [i * (i - 1) / 2 + j - 1]; // 110
		a [i + 1] = grc [i];
		s = grc [i] * grc [i] * beta [i];
		gain -= s;
		if (gain <= 0.0) {
			frameAnalysisInfo = 4;
			break;
		}
		lpcFrame -> gain = gain;
	}
end:
	const integer numberOfCoefficients = i - 1;
	lpcFrame -> a.resize (numberOfCoefficients);
	lpcFrame -> a.part (1, numberOfCoefficients)  <<=  a.part (2, i);
	lpcFrame -> nCoefficients = numberOfCoefficients; // maintain invariant
	return true;
}

/*********************** Burg method *************************************************************/

Thing_implement (SoundFrameIntoLPCFrame2Burg, SoundFrameIntoLPCFrame2, 0);

void structSoundFrameIntoLPCFrame2Burg :: initHeap () {
	SoundFrameIntoLPCFrame2Burg_Parent :: initHeap ();
	b1 = raw_VEC (soundFrameSize);
	b2 = raw_VEC (soundFrameSize);
	aa = raw_VEC (order);
}

double structSoundFrameIntoLPCFrame2Burg :: burg (VEC const& a, constVEC const& x, integer& frameAnalysisInfo) {
	const integer n = x.size, m = a.size;
	a   <<=  0.0; // always safe
	if (n <= 2) {
		a [1] = -1.0;
		return ( n == 2 ? 0.5 * (x [1] * x [1] + x [2] * x [2]) : x [1] * x [1] );
	}

	// (3)

	double p = NUMinner (x, x);

	if (p == 0.0) {
		frameAnalysisInfo = 1;
		return 0.0;
	}
	// (9)

	b1 [1] = x [1];
	for (integer j = 2; j <= n - 1; j ++)
		b1 [j] = b2 [j - 1] = x [j];
	b2 [n - 1] = x [n];

	longdouble xms = p / n;
	for (integer i = 1; i <= m; i ++) {
		// (7)

		/*
			longdouble num = 0.0, denum = 0.0;
			for (integer j = 1; j <= n - i; j ++) {
				num += b1 [j] * b2 [j];
				denum += b1 [j] * b1 [j] + b2 [j] * b2 [j];
			}
		*/
		VEC b1part = b1.part (1, n - i), b2part = b2.part (1, n - i);
		const double num = NUMinner (b1part, b2part);
		const double denum = NUMinner (b1part, b1part) + NUMinner (b2part, b2part);
		
		if (denum <= 0.0) {
			frameAnalysisInfo = 1;
			return 0.0;	// warning ill-conditioned
		}
		a [i] = 2.0 * num / denum;

		// (10)

		xms *= 1.0 - a [i] * a [i];

		// (5)

		for (integer j = 1; j <= i - 1; j ++)
			a [j] = aa [j] - a [i] * aa [i - j];

		if (i < m) {

			// (8) Watch out: i -> i+1

			for (integer j = 1; j <= i; j ++)
				aa [j] = a [j];
			for (integer j = 1; j <= n - i - 1; j ++) {
				b1 [j] -= aa [i] * b2 [j];
				b2 [j] = b2 [j + 1] - aa [i] * b1 [j + 1];
			}
		}
	}
	return double (xms);
}

bool structSoundFrameIntoLPCFrame2Burg :: inputFrameIntoOutputFrame (integer currentFrame) {
	LPC_Frame lpcFrame = & outputLPC -> d_frames[currentFrame];
	integer frameAnalysisInfo = 0;
	lpcFrame -> gain = burg (lpcFrame -> a.get(), soundFrame, frameAnalysisInfo);
	if (lpcFrame -> gain <= 0.0) {
		lpcFrame -> nCoefficients = 0;
		lpcFrame -> a.resize (lpcFrame -> nCoefficients); // maintain invariant
		return false;
	} else {
		lpcFrame -> gain *= soundFrame.size;
		for (integer i = 1; i <= lpcFrame -> nCoefficients; i ++)
			lpcFrame -> a [i] = - lpcFrame -> a [i];
		return true;
	}
}

void checkLPCAnalysisParameters_e (double sound_dx, integer sound_nx, double physicalAnalysisWidth, integer predictionOrder) {
	volatile const double physicalDuration = sound_dx * sound_nx;
	Melder_require (physicalAnalysisWidth <= physicalDuration,
		U"Your sound is shorter than two window lengths. "
		"Either your sound is too short or your window is too long.");
	// we round the minimum duration to be able to use asserterror in testing scripts.
	conststring32 minimumDurationRounded = Melder_fixed (predictionOrder * sound_dx , 5);
	const integer approximateNumberOfSamplesPerWindow = Melder_iroundDown (physicalAnalysisWidth / sound_dx);
	Melder_require (approximateNumberOfSamplesPerWindow > predictionOrder,
		U"Analysis window duration too short. For a prediction order of ", predictionOrder,
		U", the analysis window duration should be greater than ", minimumDurationRounded,
		U" s. Please increase the analysis window duration or lower the prediction order.");
}

static void Sound_and_LPC_require_equalDomainsAndSamplingPeriods (constSound me, constLPC thee) {
	Melder_require (my xmin == thy xmin && thy xmax == my xmax,
			U"The domains of the Sound and the LPC should be equal.");
	Melder_require (my dx == thy samplingPeriod,
			U"The sampling periods of the Sound and the LPC should be equal.");
}

static autoLPC LPC_createFullFromAnalysisSpecifications (constSound me, int predictionOrder, double physicalAnalysisWidth, double dt) {
	try {
		checkLPCAnalysisParameters_e (my dx, my nx, physicalAnalysisWidth, predictionOrder);		
		integer numberOfFrames;
		double t1;
		Sampled_shortTermAnalysis (me, physicalAnalysisWidth, dt, & numberOfFrames, & t1);
		autoLPC thee = LPC_create (my xmin, my xmax, numberOfFrames, dt, t1, predictionOrder, my dx);
		for (integer iframe = 1; iframe <= numberOfFrames; iframe ++)
			LPC_Frame_init (& thy d_frames [iframe], thy maxnCoefficients);
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": LPC not created from specification.");
	}
}

void Sound_into_LPC_auto2 (constSound me, mutableLPC thee, double effectiveAnalysisWidth) {
	Sound_and_LPC_require_equalDomainsAndSamplingPeriods (me, thee);
	autoSoundFrameIntoLPCFrame2Auto frameIntoFrame = SoundFrameIntoLPCFrame2Auto_create (me, thee, effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
	SampledIntoSampled_mt (frameIntoFrame.get(), 40);
}

autoLPC Sound_to_LPC_auto2 (constSound me, int predictionOrder, double effectiveAnalysisWidth, double dt, double preEmphasisFrequency) {
	try {
		const double physicalAnalysisWidth = getPhysicalAnalysisWidth2 (effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
		checkLPCAnalysisParameters_e (my dx, my nx, physicalAnalysisWidth, predictionOrder);
		autoSound emphasized = Sound_resampleAndOrPreemphasize (me, 0.0, 0, preEmphasisFrequency);
		autoLPC thee = LPC_createFullFromAnalysisSpecifications (emphasized.get(), predictionOrder, physicalAnalysisWidth, dt);
		Sound_into_LPC_auto2 (emphasized.get(), thee.get(), effectiveAnalysisWidth);
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": no LPC (auto) created.");
	}
}

#if 0
void Sound_into_LPC_covar (constSound me, mutableLPC thee, double effectiveAnalysisWidth) {
	Sound_and_LPC2_require_equalDomainsAndSamplingPeriods (me, thee);
	autoSoundFrameIntoLPCFrameCovar ws = SoundFrameIntoLPCFrameCovar_create (
		me, thee, effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
	autoSoundIntoLPCStatus status = SoundIntoLPCStatus_create (thy nx);
	autoSampledIntoSampled sis = SampledIntoSampled_create (me, thee, ws.move(), status.move());
	SampledIntoSampled_analyseThreaded (sis.get());
}

autoLPC Sound_to_LPC_covar (constSound me, int predictionOrder, double effectiveAnalysisWidth, double dt, double preEmphasisFrequency) {
	try {
		const double physicalAnalysisWidth = getPhysicalAnalysisWidth (effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
		checkLPCAnalysisParameters_e (my dx, my nx, physicalAnalysisWidth, predictionOrder);
		autoSound emphasized = Sound_resampleAndOrPreemphasize (me, 0.0, 0, preEmphasisFrequency);
		autoLPC thee = LPC_createEmptyFromAnalysisSpecifications (emphasized.get(), predictionOrder, physicalAnalysisWidth, dt);
		Sound_into_LPC_covar (emphasized.get(), thee.get(), effectiveAnalysisWidth);
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": no LPC (covar) created.");
	}
}

void Sound_into_LPC_burg (constSound me, mutableLPC thee, double effectiveAnalysisWidth) {
	Sound_and_LPC2_require_equalDomainsAndSamplingPeriods (me, thee);
	autoSoundFrameIntoLPCFrameBurg ws = SoundFrameIntoLPCFrameBurg_create (
		me, thee, effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
	autoSoundIntoLPCStatus status = SoundIntoLPCStatus_create (thy nx);
	autoSampledIntoSampled sis = SampledIntoSampled_create (me, thee, ws.move(), status.move());
	SampledIntoSampled_analyseThreaded (sis.get());
}

autoLPC Sound_to_LPC_burg (constSound me, int predictionOrder, double effectiveAnalysisWidth, double dt, double preEmphasisFrequency) {
	try {
		const double physicalAnalysisWidth = getPhysicalAnalysisWidth (effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
		checkLPCAnalysisParameters_e (my dx, my nx, physicalAnalysisWidth, predictionOrder);
		autoSound emphasized = Sound_resampleAndOrPreemphasize (me, 0.0, 0, preEmphasisFrequency);
		autoLPC thee = LPC_createEmptyFromAnalysisSpecifications (emphasized.get(), predictionOrder, physicalAnalysisWidth, dt);
		Sound_into_LPC_burg (emphasized.get(), thee.get(), effectiveAnalysisWidth);
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": no LPC (burg) created.");
	}
}

void Sound_into_LPC_marple (constSound me, mutableLPC thee, double effectiveAnalysisWidth, double tol1, double tol2) {
	Sound_and_LPC2_require_equalDomainsAndSamplingPeriods (me, thee);
	autoSoundFrameIntoLPCFrameMarple ws = SoundFrameIntoLPCFrameMarple_create (
		me, thee, effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2, tol1, tol2);
	autoSoundIntoLPCStatus status = SoundIntoLPCStatus_create (thy nx);
	autoSampledIntoSampled sis = SampledIntoSampled_create (me, thee, ws.move(), status.move());
	SampledIntoSampled_analyseThreaded (sis.get());
}

autoLPC Sound_to_LPC_marple (constSound me, int predictionOrder, double effectiveAnalysisWidth, double dt, double preEmphasisFrequency,
	double tol1, double tol2)
{
	try {
		const double physicalAnalysisWidth = getPhysicalAnalysisWidth (effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
		checkLPCAnalysisParameters_e (my dx, my nx, physicalAnalysisWidth, predictionOrder);
		autoSound emphasized = Sound_resampleAndOrPreemphasize (me, 0.0, 0, preEmphasisFrequency);
		autoLPC thee = LPC_createEmptyFromAnalysisSpecifications (emphasized.get(), predictionOrder, physicalAnalysisWidth, dt);
		Sound_into_LPC_marple (emphasized.get(), thee.get(), effectiveAnalysisWidth, tol1, tol2);
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": no LPC (marple) created.");
	}
}


void LPC_and_Sound_into_LPC_robust (constLPC inputLPC, constSound sound, mutableLPC outputlpc, double effectiveAnalysisWidth, double k_stdev,
	integer itermax, double tol, bool wantlocation)
{
	try {
		Sound_and_LPC2_require_equalDomainsAndSamplingPeriods (sound, inputLPC);	
		const double physicalAnalysisWidth = getPhysicalAnalysisWidth (effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
		double location = 0.0;
		checkLPCAnalysisParameters_e (sound -> dx, sound -> nx, physicalAnalysisWidth, outputlpc -> maxnCoefficients);
		autoLPCAndSoundFramesIntoLPCFrameRobust ws = LPCAndSoundFramesIntoLPCFrameRobust_create (inputLPC, sound, outputlpc,
			effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2, k_stdev, itermax, tol, location, wantlocation);
		autoSoundIntoLPCStatus status = SoundIntoLPCStatus_create (outputlpc -> nx);
		autoSampledIntoSampled sis = SampledIntoSampled_create (sound, outputlpc, ws.move(), status.move());
		SampledIntoSampled_analyseThreaded (sis.get());
	} catch (MelderError) {
		Melder_throw (sound, U": no LPC (robust) calculated.");
	}	
}

autoLPC LPC_and_Sound_to_LPC_robust (constLPC thee, constSound me, double effectiveAnalysisWidth, 
	double preEmphasisFrequency, double k_stdev, integer itermax, double tol, bool wantlocation)
{
	try {
		Sound_and_LPC2_require_equalDomainsAndSamplingPeriods (me, thee);
		autoSound emphasized = Sound_resampleAndOrPreemphasize (me, 0.0, 0, preEmphasisFrequency);
		autoLPC result = LPC_create (thy xmin, thy xmax, thy nx, thy dx, thy x1, thy maxnCoefficients, thy samplingPeriod);
		LPC_and_Sound_into_LPC_robust (thee, emphasized.get(), result.get(), effectiveAnalysisWidth, k_stdev, itermax, tol, wantlocation);
		return result;
	} catch (MelderError) {
		Melder_throw (me, U": no robust LPC created.");
	}	
}

void Sound_into_LPCrobust_common (constSound me, mutableLPC outputlpc, autoSoundFrameIntoLPCFrame soundIntoLPCany, double effectiveAnalysisWidth, 
	double k_stdev, integer itermax, double tol, bool wantlocation)
{
	Sound_and_LPC2_require_equalDomainsAndSamplingPeriods (me, outputlpc);
	autoLPC inputLPC = Data_copy (outputlpc);
	autoLPCAndSoundFramesIntoLPCFrameRobust lpcAndSoundIntoLPC = LPCAndSoundFramesIntoLPCFrameRobust_create (inputLPC.get(), me, outputlpc, effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2, k_stdev, itermax, tol, 0.0, wantlocation);
	autoSoundFrameIntoLPCFrameRobust soundIntoLPCrobust = SoundFrameIntoLPCFrameRobust_create (soundIntoLPCany.move(),
		lpcAndSoundIntoLPC.move());
	autoSoundIntoLPCStatus status = SoundIntoLPCStatus_create (outputlpc -> nx); // TODO adapt
	autoSampledIntoSampled sis = SampledIntoSampled_create (me, outputlpc, soundIntoLPCrobust.move(), status.move());
	SampledIntoSampled_analyseThreaded (sis.get());
}

autoLPC Sound_to_LPC_robust (constSound me, int predictionOrder, double effectiveAnalysisWidth, double dt, double preEmphasisFrequency,
	double k_stdev, integer itermax, double tol, bool wantlocation)
{
	try {
		const double physicalAnalysisWidth = getPhysicalAnalysisWidth (effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
		autoSound emphasized = Sound_resampleAndOrPreemphasize (me, 0.0, 0, preEmphasisFrequency);
		autoLPC thee = LPC_createEmptyFromAnalysisSpecifications (emphasized.get(), predictionOrder, physicalAnalysisWidth, dt);
		autoLPC output = Data_copy (thee.get());
		Sound_into_LPC_auto (emphasized.get(), thee.get(), effectiveAnalysisWidth);
		LPC_and_Sound_into_LPC_robust (thee.get(), emphasized.get(), output.get(), effectiveAnalysisWidth, k_stdev,
			itermax,  tol,  wantlocation);
		return output;
	} catch (MelderError) {
		Melder_throw (me, U": no LPC (robust) created.");
	}
}

#endif

/*********************** (inverse) filtering ******************************/

static void LPC_Frame_Sound_filter (constLPC_Frame me, mutableSound thee, integer channel) {
	const VEC y = thy z.row (channel);
	for (integer i = 1; i <= thy nx; i ++) {
		const integer m = ( i > my nCoefficients ? my nCoefficients : i - 1 );   // ppgb: what is m?
		for (integer j = 1; j <= m; j ++)
			y [i] -= my a [j] * y [i - j];
	}
}

autoSound LPC_Sound_filterInverse (constLPC me, constSound thee) {
	try {
		Melder_require (my samplingPeriod == thy dx,
			U"The sampling frequencies should be equal.");
		Melder_require (my xmin == thy xmin && thy xmax == my xmax,
			U"The domains of LPC and Sound should be equal.");
		
		autoSound him = Data_copy (thee);
		VEC source = his z.row (1);
		VEC sound = thy z.row (1);
		for (integer isamp = 1; isamp <= his nx; isamp ++) {
			const double sampleTime = Sampled_indexToX (him.get(), isamp);
			const integer frameNumber = Sampled_xToNearestIndex (me, sampleTime);
			if (frameNumber < 1 || frameNumber > my nx) {
				source [isamp] = 0.0;
				continue;
			}
			const LPC_Frame frame = & my d_frames [frameNumber];
			const integer maximumFilterDepth = frame -> nCoefficients;
			const integer maximumSoundDepth = isamp - 1;
			const integer usableDepth = std::min (maximumFilterDepth, maximumSoundDepth);
			for (integer icoef = 1; icoef <= usableDepth; icoef ++)
				source [isamp] += frame -> a [icoef] * sound [isamp - icoef];
		}
		return him;
	} catch (MelderError) {
		Melder_throw (thee, U": not inverse filtered.");
	}
}

/*
	Gain used as a constant amplitude multiplier within a frame of duration my dx.
	future alternative: convolve gain with a smoother.
*/
autoSound LPC_Sound_filter (constLPC me, constSound thee, bool useGain) {
	try {
		const double xmin = std::max (my xmin, thy xmin);
		const double xmax = std::min (my xmax, thy xmax);
		Melder_require (xmin < xmax,
			U"Domains of Sound [", thy xmin, U",", thy xmax, U"] and LPC [",
			my xmin, U",", my xmax, U"] should overlap."
		);
		/*
			Resample the sound if the sampling frequencies do not match.
		*/
		autoSound source;
		if (my samplingPeriod != thy dx) {
			source = Sound_resample (thee, 1.0 / my samplingPeriod, 50);
			thee = source.get();   // reference copy; remove at end
		}

		autoSound him = Data_copy (thee);

		const integer ifirst = std::max (1_integer, Sampled_xToHighIndex (thee, xmin));
		const integer ilast = std::min (Sampled_xToLowIndex (thee, xmax), thy nx);
		for (integer isamp = ifirst; isamp <= ilast; isamp ++) {
			const double sampleTime = Sampled_indexToX (him.get(), isamp);
			const integer frameNumber = Sampled_xToNearestIndex (me, sampleTime);
			if (frameNumber < 1 || frameNumber > my nx) {
				his z [1] [isamp] = 0.0;
				continue;
			}
			const LPC_Frame frame = & my d_frames [frameNumber];
			const integer maximumFilterDepth = frame -> nCoefficients;
			const integer maximumSourceDepth = isamp - 1;
			const integer usableDepth = std::min (maximumFilterDepth, maximumSourceDepth);
			for (integer icoef = 1; icoef <= usableDepth; icoef ++)
				his z [1] [isamp] -= frame -> a [icoef] * his z [1] [isamp - icoef];
		}
		/*
			Make samples before first frame and after last frame zero.
		*/
		for (integer isamp = 1; isamp < ifirst; isamp ++)
			his z [1] [isamp] = 0.0;
		for (integer isamp = ilast + 1; isamp <= his nx; isamp ++)
			his z [1] [isamp] = 0.0;

		if (useGain) {
			for (integer isamp = ifirst; isamp <= ilast; isamp ++) {
				const double sampleTime = Sampled_indexToX (him.get(), isamp);
				const double realFrameNumber = Sampled_xToIndex (me, sampleTime);
				const integer leftFrameNumber = Melder_ifloor (realFrameNumber);
				const integer rightFrameNumber = leftFrameNumber + 1;
				const double phase = realFrameNumber - leftFrameNumber;
				if (rightFrameNumber < 1 || leftFrameNumber > my nx)
					his z [1] [isamp] = 0.0;
				else if (rightFrameNumber == 1)
					his z [1] [isamp] *= sqrt (my d_frames [1]. gain) * phase;
				else if (leftFrameNumber == my nx)
					his z [1] [isamp] *= sqrt (my d_frames [my nx]. gain) * (1.0 - phase);
				else 
					his z [1] [isamp] *=
							phase * sqrt (my d_frames [rightFrameNumber]. gain) +
							(1.0 - phase) * sqrt (my d_frames [leftFrameNumber]. gain);
			}
		}
		return him;
	} catch (MelderError) {
		Melder_throw (thee, U": not filtered.");
	}
}

void LPC_Sound_filterWithFilterAtTime_inplace (constLPC me, mutableSound thee, integer channel, double time) {
	integer frameIndex = Sampled_xToNearestIndex (me, time);
	Melder_clip (1_integer, & frameIndex, my nx);   // constant extrapolation
	if (channel > thy ny)
		channel = 1;
	Melder_require (frameIndex > 0 && frameIndex <= my nx,
		U"Frame should be in the range [1, ", my nx, U"].");

	if (channel > 0)
		LPC_Frame_Sound_filter (& my d_frames [frameIndex], thee, channel);
	else
		for (integer ichan = 1; ichan <= thy ny; ichan ++)
			LPC_Frame_Sound_filter (& my d_frames [frameIndex], thee, ichan);
}

autoSound LPC_Sound_filterWithFilterAtTime (constLPC me, constSound thee, integer channel, double time) {
	try {
		autoSound him = Data_copy (thee);
		LPC_Sound_filterWithFilterAtTime_inplace (me, him.get(), channel, time);
		return him;
	} catch (MelderError) {
		Melder_throw (thee, U": not filtered.");
	}
}

void LPC_Sound_filterInverseWithFilterAtTime_inplace (constLPC me, mutableSound thee, integer channel, double time) {
	try {
		integer frameIndex = Sampled_xToNearestIndex (me, time);
		Melder_clip (1_integer, & frameIndex, my nx);   // constant extrapolation
		if (channel > thy ny)
			channel = 1;
		LPC_Frame lpc = & my d_frames [frameIndex];
		autoVEC workspace = raw_VEC (lpc -> nCoefficients);
		if (channel > 0)
			VECfilterInverse_inplace (thy z.row (channel), lpc -> a.get(), workspace.get());
		else
			for (integer ichan = 1; ichan <= thy ny; ichan ++)
				VECfilterInverse_inplace (thy z.row (ichan), lpc -> a.get(), workspace.get());
	} catch (MelderError) {
		Melder_throw (thee, U": not inverse filtered.");
	}
}

autoSound LPC_Sound_filterInverseWithFilterAtTime (constLPC me, constSound thee, integer channel, double time) {
	try {
		autoSound him = Data_copy (thee);
		LPC_Sound_filterInverseWithFilterAtTime_inplace (me, him.get(), channel, time);
		return him;
	} catch (MelderError) {
		Melder_throw (thee, U": not inverse filtered.");
	}
}

/* End of file Sound_and_LPC2.cpp */
