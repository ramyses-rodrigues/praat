#ifndef _NUMrandom_h_
#define _NUMrandom_h_
/* NUMrandom.h
 *
 * Copyright (C) 1992-2018,2020,2025 Paul Boersma
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

/********** Random numbers **********/

void NUMrandom_initializeSafelyAndUnpredictably ();
void NUMrandom_initializeWithSeedUnsafelyButPredictably (uint64 seed);

constexpr integer NUMrandom_maximumNumberOfParallelThreads = 40;
// The number of random channels will be three more than this:
//    thread 0: the thread in which Praat commands and Praat scripts run
//    threads 1 through 40: extra threads for parallellization in analyses (the master thread stays 0)
//    thread 41: user-interface randomization (not really a thread)
//    thread 42: spare "thread"
constexpr integer NUMrandom_numberOfChannels = NUMrandom_maximumNumberOfParallelThreads + 3;

void NUMrandom_setChannel (integer channelNumber);
// To be set by each parallel thread to a value between 0 and NUMrandom_numberOfChannels - 1.
integer NUMrandom_getChannel ();
// Rarely useful, but who knows. Only to be called after NUMrandom_setChannel() has been called.

double NUMrandomFraction ();

double NUMrandomUniform (double lowest, double highest);

integer NUMrandomInteger (integer lowest, integer highest);

bool NUMrandomBernoulli (double probability);
double NUMrandomBernoulli_real (double probability);

double NUMrandomGauss (double mean, double standardDeviation);

double NUMrandomPoisson (double mean);

uint32 NUMhashString (conststring32 string);

/* End of file NUMrandom.h */
#endif
