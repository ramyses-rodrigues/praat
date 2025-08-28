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

static void Sound_to_LPC_common_e (constSound inputSound, int predictionOrder, double effectiveAnalysisWidth, double dt,
	double preEmphasisFrequency, autoSound& outputSound, autoLPC& outputLPC)
{
	try {
		const double physicalAnalysisWidth = getPhysicalAnalysisWidth2 (effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
		checkLPCAnalysisParameters_e (inputSound -> dx, inputSound -> nx, physicalAnalysisWidth, predictionOrder);
		integer numberOfFrames;
		double t1;
		autoSound sound = Sound_resampleAndOrPreemphasize (inputSound, 0.0, 0, preEmphasisFrequency);
		outputSound = sound.move();
		Sampled_shortTermAnalysis (outputSound.get(), physicalAnalysisWidth, dt, & numberOfFrames, & t1);
		autoLPC lpc = LPC_create (outputSound -> xmin, outputSound -> xmax, numberOfFrames, dt, t1, predictionOrder, outputSound -> dx);
		for (integer iframe = 1; iframe <= numberOfFrames; iframe ++) {
			LPC_Frame lpcFrame = & lpc -> d_frames [iframe];
			LPC_Frame_init (lpcFrame, lpc -> maxnCoefficients);
		}
		outputLPC = lpc.move();
	} catch (MelderError) {
		Melder_throw (inputSound, U": Sound and/or LPC not created from specification.");
	}
}

static autoLPC LPC_createFullFromAnalysisSpecifications (constSound me, int predictionOrder, double physicalAnalysisWidth, double dt) {
	try {
		checkLPCAnalysisParameters_e (my dx, my nx, physicalAnalysisWidth, predictionOrder);		
		integer numberOfFrames;
		double t1;
		Sampled_shortTermAnalysis (me, physicalAnalysisWidth, dt, & numberOfFrames, & t1);
		autoLPC thee = LPC_create (my xmin, my xmax, numberOfFrames, dt, t1, predictionOrder, my dx);
		for (integer iframe = 1; iframe <= numberOfFrames; iframe ++) {
			LPC_Frame lpcFrame = & thy d_frames [iframe];
			LPC_Frame_init (lpcFrame, thy maxnCoefficients);
		}
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": LPC not created from specification.");
	}
}

Thing_implement (SoundFrameIntoLPCFrame2, SoundFrameIntoSampledFrame2, 0);

void structSoundFrameIntoLPCFrame2 :: initBasicSoundFrameAndLPC (constSound inputSound, mutableLPC outputLPC, 
	double effectiveAnalysisWidth, kSound_windowShape windowShape)
{
	initBasicSoundFrame (inputSound, outputLPC, effectiveAnalysisWidth, windowShape);
	our outputLPC = outputLPC;
	order = outputLPC -> maxnCoefficients;
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
	frameAnalysisInfo = 0;

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

autoSoundFrameIntoLPCFrame2Auto SoundFrameIntoLPCFrame2Auto_create (constSound inputSound, mutableLPC outputLPC, double effectiveAnalysisWidth, kSound_windowShape windowShape) {
	try {
		autoSoundFrameIntoLPCFrame2Auto me = Thing_new (SoundFrameIntoLPCFrame2Auto);
		my initBasicSoundFrameAndLPC (inputSound, outputLPC, effectiveAnalysisWidth, windowShape);
		return me;
	} catch (MelderError) {
		Melder_throw (U"Cannot create SoundFrameIntoLPCFrame2Auto");
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

/*********************** Covariance method *************************************************************/

Thing_implement (SoundFrameIntoLPCFrame2Covar, SoundFrameIntoLPCFrame2, 0);

void structSoundFrameIntoLPCFrame2Covar :: initBasicSoundFrameAndLPC (constSound inputSound, mutableLPC outputLPC, double effectiveAnalysisWidth, kSound_windowShape windowShape) 
{
	SoundFrameIntoLPCFrame2Covar_Parent :: initBasicSoundFrameAndLPC (inputSound, outputLPC, effectiveAnalysisWidth, windowShape);
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
	frameAnalysisInfo = 0;
	
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

autoSoundFrameIntoLPCFrame2Covar SoundFrameIntoLPCFrame2Covar_create (constSound input, mutableLPC outputLPC, double effectiveAnalysisWidth, kSound_windowShape windowShape) {
	try {
		autoSoundFrameIntoLPCFrame2Covar me = Thing_new (SoundFrameIntoLPCFrame2Covar);
		my initBasicSoundFrameAndLPC (input, outputLPC, effectiveAnalysisWidth, windowShape);
		return me;
	} catch (MelderError) {
		Melder_throw (U"Cannot create SoundFrameIntoLPCFrame2Covar");
	}
}

void Sound_into_LPC_covar2 (constSound me, mutableLPC thee, double effectiveAnalysisWidth) {
	Sound_and_LPC_require_equalDomainsAndSamplingPeriods (me, thee);
	autoSoundFrameIntoLPCFrame2Covar frameIntoFrame = SoundFrameIntoLPCFrame2Covar_create (me, thee, effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
	SampledIntoSampled_mt (frameIntoFrame.get(), 40);
}

autoLPC Sound_to_LPC_covar2 (constSound me, int predictionOrder, double effectiveAnalysisWidth, double dt, double preEmphasisFrequency) {
	try {
		const double physicalAnalysisWidth = getPhysicalAnalysisWidth2 (effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
		checkLPCAnalysisParameters_e (my dx, my nx, physicalAnalysisWidth, predictionOrder);
		autoSound emphasized = Sound_resampleAndOrPreemphasize (me, 0.0, 0, preEmphasisFrequency);
		autoLPC thee = LPC_createFullFromAnalysisSpecifications (emphasized.get(), predictionOrder, physicalAnalysisWidth, dt);
		Sound_into_LPC_covar2 (emphasized.get(), thee.get(), effectiveAnalysisWidth);
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": no LPC (covar) created.");
	}
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

autoSoundFrameIntoLPCFrame2Burg SoundFrameIntoLPCFrame2Burg_create (constSound inputSound, mutableLPC outputLPC, double effectiveAnalysisWidth, kSound_windowShape windowShape) {
	try {
		autoSoundFrameIntoLPCFrame2Burg me = Thing_new (SoundFrameIntoLPCFrame2Burg);
		my initBasicSoundFrameAndLPC (inputSound, outputLPC, effectiveAnalysisWidth, windowShape);
		return me;
	} catch (MelderError) {
		Melder_throw (U"Cannot create SoundFrameIntoLPCFrame2Burg");
	}
}

void Sound_into_LPC_burg2 (constSound me, mutableLPC thee, double effectiveAnalysisWidth) {
	Sound_and_LPC_require_equalDomainsAndSamplingPeriods (me, thee);
	autoSoundFrameIntoLPCFrame2Burg frameIntoFrame = SoundFrameIntoLPCFrame2Burg_create (me, thee, effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
	SampledIntoSampled_mt (frameIntoFrame.get(), 40);
}

autoLPC Sound_to_LPC_burg2 (constSound me, int predictionOrder, double effectiveAnalysisWidth, double dt, double preEmphasisFrequency) {
	try {
		const double physicalAnalysisWidth = getPhysicalAnalysisWidth2 (effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
		checkLPCAnalysisParameters_e (my dx, my nx, physicalAnalysisWidth, predictionOrder);
		autoSound emphasized = Sound_resampleAndOrPreemphasize (me, 0.0, 0, preEmphasisFrequency);
		autoLPC thee = LPC_createFullFromAnalysisSpecifications (emphasized.get(), predictionOrder, physicalAnalysisWidth, dt);
		Sound_into_LPC_burg2 (emphasized.get(), thee.get(), effectiveAnalysisWidth);
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": no LPC (burg) created.");
	}
}


/*********************** Marple method *************************************************************/

Thing_implement (SoundFrameIntoLPCFrame2Marple, SoundFrameIntoLPCFrame2, 0);

void structSoundFrameIntoLPCFrame2Marple :: initHeap () {
	SoundFrameIntoLPCFrame2Marple_Parent :: initHeap ();
	c = raw_VEC (orderp1);
	d = raw_VEC (orderp1);
	r = raw_VEC (orderp1);
}

bool structSoundFrameIntoLPCFrame2Marple :: inputFrameIntoOutputFrame (integer iframe) {
	const integer mmax = order, n = soundFrame.size;
	LPC_Frame lpcFrame = & outputLPC -> d_frames [iframe];
	VEC x = soundFrame;
	
	frameAnalysisInfo = 0;
	VEC c = a.get(); // yes 'a'
	VEC a = lpcFrame -> a.get();

	double gain = 0.0, e0 = 2.0 * NUMsum2 (x);
	integer m = 1;
	if (e0 == 0.0) {
		lpcFrame -> nCoefficients = 0;
		lpcFrame -> a.resize (lpcFrame -> nCoefficients); // maintain invariant
		lpcFrame -> gain = gain;
		frameAnalysisInfo = 1;
		return false;
	}
	double q1 = 1.0 / e0;
	double q2 = q1 * x [1], q = q1 * x [1] * x [1], w = q1 * x [n] * x [n];
	double v = q, u = w;
	double den = 1.0 - q - w;
	double q4 = 1.0 / den, q5 = 1.0 - q, q6 = 1.0 - w;
	double h = q2 * x [n], s = h;
	gain = e0 * den;
	q1 = 1.0 / gain;
	c [1] = q1 * x [1];
	d [1] = q1 * x [n];
	double s1 = 0.0;
	for (integer k = 1; k <= n - 1; k ++)
		s1 += x [k + 1] * x [k];
	r [1] = 2.0 * s1;
	a [1] = - q1 * r [1];
	gain *= (1.0 - a [1] * a [1]);
	while (m < mmax) {
		const double eOld = gain;
		double f = x [m + 1], b = x [n - m]; // n-1 ->n-m
		for (integer k = 1; k <= m; k ++) {
			// n-1 -> n-m
			f += x [m + 1 - k] * a [k];
			b += x [n - m + k] * a [k];
		}
		q1 = 1.0 / gain;
		q2 = q1 * f;
		const double q3 = q1 * b;
		for (integer k = m; k >= 1; k--) {
			c [k + 1] = c [k] + q2 * a [k];
			d [k + 1] = d [k] * q3 * a [k];
		}
		c [1] = q2;
		d [1] = q3;
		const double q7 = s * s;
		double y1 = f * f;
		const double y2 = v * v;
		const double y3 = b * b;
		const double y4 = u * u;
		double y5 = 2.0 * h * s;
		q += y1 * q1 + q4 * (y2 * q6 + q7 * q5 + v * y5);
		w += y3 * q1 + q4 * (y4 * q5 + q7 * q6 + u * y5);
		h = s = u = v = 0.0;
		for (integer k = 0; k <= m; k ++) {
			h += x [n - m + k] * c [k + 1];
			s += x [n - k] * c [k + 1];
			u += x [n - k] * d [k + 1];
			v += x [k + 1] * c [k + 1];
		}
		q5 = 1.0 - q;
		q6 = 1.0 - w;
		den = q5 * q6 - h * h;
		if (den <= 0.0) {
			frameAnalysisInfo = 2;
			goto end; // 2: ill-conditioning
		}
		q4 = 1.0 / den;
		q1 *= q4;
		const double alf = 1.0 / (1.0 + q1 * (y1 * q6 + y3 * q5 + 2.0 * h * f * b));
		gain *= alf;
		y5 = h * s;
		double c1 = q4 * (f * q6 + b * h);
		double c2 = q4 * (b * q5 + h * f);
		const double c3 = q4 * (v * q6 + y5);
		const double c4 = q4 * (s * q5 + v * h);
		const double c5 = q4 * (s * q6 + h * u);
		const double c6 = q4 * (u * q5 + y5);
		for (integer k = 1; k <= m; k ++)
			a [k] = alf * (a [k] + c1 * c [k + 1] + c2 * d [k + 1]);
		for (integer k = 1; k <= m / 2 + 1; k ++) {
			s1 = c [k];
			const double s2 = d [k], s3 = c [m + 2 - k], s4 = d [m + 2 - k];

			c [k] += c3 * s3 + c4 * s4;
			d [k] += c5 * s3 + c6 * s4;
			if (m + 2 - k == k)
				continue;
			c [m + 2 - k] += c3 * s1 + c4 * s2;
			d [m + 2 - k] += c5 * s1 + c6 * s2;
		}
		m ++;
		c1 = x [n + 1 - m];
		c2 = x [m];
		double delta = 0.0;
		for (integer k = m - 1; k >= 1; k--) {
			r [k + 1] = r [k] - x [n + 1 - k] * c1 - x [k] * c2;
			delta += r [k + 1] * a [k];
		}
		s1 = 0.0;
		for (integer k = 1; k <= n - m; k ++)
			s1 += x [k + m] * x [k];
		r [1] = 2.0 * s1;
		delta += r [1];
		q2 = - delta / gain;
		a [m] = q2;
		for (integer k = 1; k <= m / 2; k ++) {
			s1 = a [k];
			a [k] += q2 * a [m - k];
			if (k == m - k)
				continue;
			a [m - k] += q2 * s1;
		}
		y1 = q2 * q2;
		gain *= 1.0 - y1;
		if (y1 >= 1.0) {
			frameAnalysisInfo = 3;
			break; // |a [m]| > 1
		}
		if (gain < e0 * tol1) {
			frameAnalysisInfo = 4;
			break;
		}
		if (eOld - gain < eOld * tol2) {
			frameAnalysisInfo = 5;
			break;
		}
	}
end:
	lpcFrame -> gain = gain * 0.5;   // because e0 is twice the energy
	lpcFrame -> a.resize (m);
	lpcFrame -> nCoefficients = m;   // maintain invariant
	return frameAnalysisInfo == 0 || frameAnalysisInfo == 4 || frameAnalysisInfo == 5;
}


autoSoundFrameIntoLPCFrame2Marple SoundFrameIntoLPCFrame2Marple_create (constSound inputSound, mutableLPC outputLPC, double effectiveAnalysisWidth, kSound_windowShape windowShape, double tol1, double tol2) {
	try {
		autoSoundFrameIntoLPCFrame2Marple me = Thing_new (SoundFrameIntoLPCFrame2Marple);
		my initBasicSoundFrameAndLPC (inputSound, outputLPC, effectiveAnalysisWidth, windowShape);
		my tol1 = tol1;
		my tol2 = tol2;
		return me;
	} catch (MelderError) {
		Melder_throw (U"Cannot create SoundFrameIntoLPCFrame2Burg");
	}
}

void Sound_into_LPC_marple2 (constSound me, mutableLPC thee, double effectiveAnalysisWidth, double tol1, double tol2) {
	Sound_and_LPC_require_equalDomainsAndSamplingPeriods (me, thee);
	autoSoundFrameIntoLPCFrame2Marple frameIntoFrame = SoundFrameIntoLPCFrame2Marple_create (
		me, thee, effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2, tol1, tol2);
	SampledIntoSampled_mt (frameIntoFrame.get(), 40);
}

autoLPC Sound_to_LPC_marple2 (constSound me, int predictionOrder, double effectiveAnalysisWidth, double dt, 
	double preEmphasisFrequency, double tol1, double tol2)
{
	try {
		const double physicalAnalysisWidth = getPhysicalAnalysisWidth2 (effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
		checkLPCAnalysisParameters_e (my dx, my nx, physicalAnalysisWidth, predictionOrder);
		autoSound emphasized = Sound_resampleAndOrPreemphasize (me, 0.0, 0, preEmphasisFrequency);
		autoLPC thee = LPC_createFullFromAnalysisSpecifications (emphasized.get(), predictionOrder, physicalAnalysisWidth, dt);
		Sound_into_LPC_marple2 (emphasized.get(), thee.get(), effectiveAnalysisWidth, tol1, tol2);
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": no LPC (marple) created.");
	}
}

/*********************** Robust method (LPC & Sound) *************************************************************/

Thing_implement (LPCAndSoundFramesIntoLPCFrameRobust2, SoundFrameIntoLPCFrame2, 0);


void structLPCAndSoundFramesIntoLPCFrameRobust2 :: initHeap () {
	LPCAndSoundFramesIntoLPCFrameRobust2_Parent :: initHeap ();
	error = raw_VEC (soundFrameSize);
	sampleWeights = raw_VEC (soundFrameSize);
	coefficients = raw_VEC (order);
	covariancesw = raw_VEC (order);
	covarmatrixw = raw_MAT (order, order);
	svd = SVD_create (order, order);
	SVD_setTolerance (svd.get(), tolSVD);
	computedSVDworksize = SVD_getWorkspaceSize (svd.get());
	svdwork1 = raw_VEC (computedSVDworksize);
	svdwork2 = raw_VEC (order);
	filterMemory = raw_VEC (order);
	huberwork = raw_VEC (soundFrameSize);
}

void structLPCAndSoundFramesIntoLPCFrameRobust2 :: resize () {
	if (currentPredictionOrder == svd -> numberOfColumns)
		return;
	Melder_assert (currentPredictionOrder <= order);
	coefficients.resize (currentPredictionOrder);
	covariancesw.resize (currentPredictionOrder);
	covarmatrixw.resize (currentPredictionOrder, currentPredictionOrder);
	SVD_resizeWithinOldBounds (svd.get(), order, order, currentPredictionOrder, currentPredictionOrder);
}

void structLPCAndSoundFramesIntoLPCFrameRobust2 :: setSampleWeights () {
	const double kstdev = k_stdev * scale;
	for (integer isamp = 1 ; isamp <= error.size; isamp ++) {
		const double absDiff = fabs (error [isamp] - location);
		sampleWeights [isamp] = ( absDiff <= kstdev ? 1.0 : kstdev / absDiff );
	}
}

void structLPCAndSoundFramesIntoLPCFrameRobust2 :: setCovariances () {
	MATVU covar = MATVU (covarmatrixw.part (1, currentPredictionOrder, 1, currentPredictionOrder));
	for (integer i = 1; i <= currentPredictionOrder; i ++) {
		for (integer j = i; j <= currentPredictionOrder; j ++) {
			longdouble cv1 = 0.0;
			/*
				The following inner loop will need the most CPU time of all the robust calculations
				
					for (integer k = my currentPredictionOrder + 1; k <= s.size; k ++)
						cv1 += s [k - j] * s [k - i] *  my sampleWeights [k];
				
				The following code speeds it up from 23.24% to 18.32% of the total CPU time used
				(sound with 44100 Hz sampling frequency, 0.025 s window length)
			*/
			const double *skmj = & soundFrame [currentPredictionOrder - j];
			const double *skmi = & soundFrame [currentPredictionOrder - i];
			const double *sw = & sampleWeights [currentPredictionOrder];
			for (integer k = 1; k <= soundFrame.size - currentPredictionOrder; k ++)
				cv1 += *++skmj * *++skmi * *++sw;
			covar [i] [j] = covar [j] [i] = (double) cv1;
		}
		longdouble cv2 = 0.0;
		for (integer k = currentPredictionOrder + 1; k <= soundFrame.size; k ++)
			cv2 += soundFrame [k - i] * soundFrame [k] *  sampleWeights [k];
		covariancesw [i] = - cv2;
	}
}

void structLPCAndSoundFramesIntoLPCFrameRobust2 :: solvelpc () {
	svd -> u.all()  <<=  covarmatrixw.all();
	svdwork2.resize (currentPredictionOrder);
	SVD_compute (svd.get(), svdwork1.get());
	SVD_solve_preallocated (svd.get(), covariancesw.get(), coefficients.get(), svdwork2.get());
	coefficients.resize (currentPredictionOrder); // maintain invariant
}

bool structLPCAndSoundFramesIntoLPCFrameRobust2 :: inputFrameIntoOutputFrame (integer currentFrame) {
	LPC_Frame inputLPCFrame = & inputLPC -> d_frames [currentFrame];
	LPC_Frame outputLPCFrame = & outputLPC -> d_frames [currentFrame];
	currentPredictionOrder = inputLPCFrame -> nCoefficients;
	for (integer i = 1; i <= inputLPCFrame -> nCoefficients; i ++)
		outputLPCFrame -> a[i] = inputLPCFrame -> a [i];
	if (currentPredictionOrder == 0) // is empty frame ?
		return true;
	outputLPCFrame -> gain = inputLPCFrame -> gain;
	
	VEC inout_a = outputLPCFrame -> a.part (1, currentPredictionOrder);
	iter = 0;
	scale = 1e308;
	bool farFromScale = true;
	resize ();
	filterMemory.resize (currentPredictionOrder);
	frameAnalysisInfo = 0;
	do {
		const double previousScale = scale;
		error.all()  <<=  soundFrame;
		VECfilterInverse_inplace (error.get(), inout_a, filterMemory.get());
		NUMstatistics_huber (error.get(), & location, wantlocation, & scale, wantscale, k_stdev, tol1, huber_iterations, huberwork.get());
		setSampleWeights ();

		setCovariances ();
		/*
			Solve C a = [-] c
		*/
		try {
			solvelpc ();
			inout_a  <<=  coefficients.all();
			farFromScale = ( fabs (scale - previousScale) > std::max (tol1 * fabs (scale), NUMeps) );
		} catch (MelderError) {
			Melder_clearError (); // No change could be made
			frameAnalysisInfo = 2; // solvelpc in error
			inputLPCFrame -> copy (outputLPCFrame);
			return false;
		}
	} while (++ iter < itermax && farFromScale);
	frameAnalysisInfo = 3; // maximum number of iterations
	return true;
}

autoLPCAndSoundFramesIntoLPCFrameRobust2 LPCAndSoundFramesIntoLPCFrameRobust2_create (constLPC inputLPC,
	constSound inputSound, mutableLPC outputLPC, double effectiveAnalysisWidth, kSound_windowShape windowShape, double k_stdev, integer itermax, double tol1, bool wantlocation)
{
	try {
		autoLPCAndSoundFramesIntoLPCFrameRobust2 me = Thing_new (LPCAndSoundFramesIntoLPCFrameRobust2);
		my initBasicSoundFrameAndLPC (inputSound, outputLPC, effectiveAnalysisWidth, windowShape);
		my inputLPC = inputLPC;
		my k_stdev = k_stdev;
		my itermax = itermax;
		my tol1 = tol1;
		my wantlocation = wantlocation;
		return me;
	} catch (MelderError) {
		Melder_throw (U"Cannot create SoundFrameIntoLPCFrame2Robust");
	}
}

void LPC_and_Sound_into_LPC_robust2 (constLPC inputLPC, constSound inputSound, mutableLPC outputLPC, double effectiveAnalysisWidth,
	double preEmphasisFrequency, double k_stdev, integer itermax, double tol, bool wantlocation)
{
	Sound_and_LPC_require_equalDomainsAndSamplingPeriods (inputSound, inputLPC);	
	const double physicalAnalysisWidth = getPhysicalAnalysisWidth2 (effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2);
	checkLPCAnalysisParameters_e (inputSound -> dx, inputSound -> nx, physicalAnalysisWidth, outputLPC -> maxnCoefficients);
	autoLPCAndSoundFramesIntoLPCFrameRobust2 frameIntoFrame = LPCAndSoundFramesIntoLPCFrameRobust2_create (inputLPC, 
		inputSound, outputLPC, effectiveAnalysisWidth, kSound_windowShape::GAUSSIAN_2, k_stdev, itermax, tol, wantlocation);
	SampledIntoSampled_mt (frameIntoFrame.get(), 40);
}

autoLPC LPC_and_Sound_to_LPC_robust2 (LPC inputLPC, constSound inputSound, double effectiveAnalysisWidth, 
	double preEmphasisFrequency, double k_stdev, integer itermax, double tol, bool wantlocation)
{
	try {
		Sound_and_LPC_require_equalDomainsAndSamplingPeriods (inputSound, inputLPC);
		autoSound emphasizedSound = Sound_resampleAndOrPreemphasize (inputSound, 0.0, 0, preEmphasisFrequency);
		autoLPC outputLPC = Data_copy (inputLPC);
		LPC_and_Sound_into_LPC_robust2 (inputLPC, emphasizedSound.get(), outputLPC.get(), effectiveAnalysisWidth,
			preEmphasisFrequency, k_stdev, itermax, tol, wantlocation);
		return outputLPC;
	} catch (MelderError) {
		Melder_throw (U"");
	}
}

/*********************** Robust method (Sound) *************************************************************/

Thing_implement (SoundFrameIntoLPCFrameRobust2, SoundFrameIntoLPCFrame2Auto, 0);

void structSoundFrameIntoLPCFrameRobust2 :: initHeap () {
	SoundFrameIntoLPCFrameRobust2_Parent :: initHeap ();
	lpcAndSoundFramesIntoLPCframe -> initHeap ();
}

bool structSoundFrameIntoLPCFrameRobust2 :: inputFrameIntoOutputFrame (integer iframe) {
	bool firstPart = SoundFrameIntoLPCFrameRobust2_Parent :: inputFrameIntoOutputFrame (iframe);
	bool secondPart = true;
	if (firstPart) {
		lpcAndSoundFramesIntoLPCframe -> soundFrame = our soundFrame;
		secondPart = lpcAndSoundFramesIntoLPCframe -> inputFrameIntoOutputFrame (iframe);
	}
	return firstPart && secondPart;
}

autoSoundFrameIntoLPCFrameRobust2 SoundFrameIntoLPCFrameRobust2_create (constSound inputSound, mutableLPC outputLPC, mutableLPC inputLPC,
	double effectiveAnalysisWidth, kSound_windowShape windowShape, double k_stdev, integer itermax, double tol, bool wantlocation)
{
	try {
		autoSoundFrameIntoLPCFrameRobust2 me = Thing_new (SoundFrameIntoLPCFrameRobust2);
		// TODO other initBasic ???
		my initBasicSoundFrameAndLPC (inputSound, outputLPC, effectiveAnalysisWidth, windowShape);
		my lpcAndSoundFramesIntoLPCframe = LPCAndSoundFramesIntoLPCFrameRobust2_create (inputLPC, inputSound, outputLPC,
			effectiveAnalysisWidth, windowShape, k_stdev, itermax, tol, wantlocation);
		my lpcAndSoundFramesIntoLPCframe -> initBasicSoundFrameAndLPC (inputSound, inputLPC, effectiveAnalysisWidth, windowShape);
		return me;
	} catch (MelderError) {
		Melder_throw (U"Cannot create SoundFrameIntoLPCFrame2Robust");
	}
}

void Sound_into_LPC_robust2 (constSound inputSound, mutableLPC outputLPC, double effectiveAnalysisWidth, kSound_windowShape windowShape, double k_stdev, integer itermax, double tol, bool wantlocation)
{
	autoLPC inputLPC = Data_copy (outputLPC);
	autoSoundFrameIntoLPCFrameRobust2 frameIntoFrame = SoundFrameIntoLPCFrameRobust2_create (inputSound, outputLPC, inputLPC.get(),
		effectiveAnalysisWidth, windowShape, k_stdev, itermax, tol, wantlocation);
	SampledIntoSampled_mt (frameIntoFrame.get(), 40);
}

autoLPC Sound_to_LPC_robust2 (constSound me, int predictionOrder, double effectiveAnalysisWidth, double dt, 
	double preEmphasisFrequency, kSound_windowShape windowShape, double k_stdev, integer itermax, double tol, bool wantlocation)
{
	try {
		autoSound emphasizedSound;
		autoLPC outputLPC;
		Sound_to_LPC_common_e (me, predictionOrder, effectiveAnalysisWidth, dt, preEmphasisFrequency, emphasizedSound, outputLPC);
		Sound_into_LPC_robust2 (emphasizedSound.get(), outputLPC.get(), effectiveAnalysisWidth, windowShape, k_stdev, itermax, tol, wantlocation);
		return outputLPC;
	} catch (MelderError) {
		Melder_throw (me, U": no LPC (robust) created.");
	}
}

#if 0



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
