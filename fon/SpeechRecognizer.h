#ifndef _SpeechRecognizer_h_
#define _SpeechRecognizer_h_
/* SpeechRecognizer.h
 *
 * Copyright (C) 2025,2026 Anastasia Shchupak
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

#include "Sound.h"
#include <set>

struct structSpeechRecognizer;
struct whisper_context;
struct whisper_vad_context;
struct whisper_vad_segments;
struct diarize_context;

struct autoWhisperContext {
	whisper_context *ptr;

	autoWhisperContext (whisper_context *p = nullptr) : ptr(p) {}
	~autoWhisperContext ();

	autoWhisperContext (const autoWhisperContext&) = delete;
	autoWhisperContext& operator= (const autoWhisperContext&) = delete;

	autoWhisperContext (autoWhisperContext&& other) noexcept : ptr (other.ptr) {
		other.ptr = nullptr;
	}
	autoWhisperContext& operator= (autoWhisperContext&& other) noexcept;

	[[nodiscard]]
	whisper_context *get () const { return ptr; }
};

struct autoWhisperVadContext {
	whisper_vad_context *ptr;

	autoWhisperVadContext (whisper_vad_context *p = nullptr) : ptr (p) {}
	~autoWhisperVadContext ();

	autoWhisperVadContext (const autoWhisperVadContext&) = delete;
	autoWhisperVadContext& operator= (const autoWhisperVadContext&) = delete;

	autoWhisperVadContext (autoWhisperVadContext&& other) noexcept : ptr (other.ptr) {
		other.ptr = nullptr;
	}
	autoWhisperVadContext& operator= (autoWhisperVadContext&& other) noexcept;

	[[nodiscard]]
	whisper_vad_context *get () const { return ptr; }
};

struct autoWhisperVadSegments {
	whisper_vad_segments *ptr;

	autoWhisperVadSegments (whisper_vad_segments *p = nullptr) : ptr (p) {}
	~autoWhisperVadSegments ();

	autoWhisperVadSegments (const autoWhisperVadSegments&) = delete;
	autoWhisperVadSegments& operator= (const autoWhisperVadSegments&) = delete;

	autoWhisperVadSegments (autoWhisperVadSegments && other) noexcept : ptr (other.ptr) {
		other.ptr = nullptr;
	}
	autoWhisperVadSegments& operator= (autoWhisperVadSegments && other) noexcept;

	[[nodiscard]]
	whisper_vad_segments *get () const { return ptr; }
};

struct autoDiarizeContext {
	diarize_context *ptr;

	autoDiarizeContext (diarize_context *p = nullptr) : ptr(p) {}
	~autoDiarizeContext ();

	autoDiarizeContext (const autoDiarizeContext&) = delete;
	autoDiarizeContext& operator= (const autoDiarizeContext&) = delete;

	autoDiarizeContext (autoDiarizeContext&& other) noexcept : ptr(other.ptr) {
		other.ptr = nullptr;
	}
	autoDiarizeContext& operator= (autoDiarizeContext&& other) noexcept;

	[[nodiscard]]
	diarize_context *get () const { return ptr; }
};

/*
	Default transcription parameters for the UI.
*/
constexpr conststring32 theSpeechRecognizerDefaultModelName = U"ggml-base.bin";
constexpr conststring32 theSpeechRecognizerDefaultLanguageName = U"Autodetect language";
constexpr bool theSpeechRecognizerDefaultIncludeWords = true;
constexpr bool theSpeechRecognizerDefaultIncludeDiarization = false;
constexpr bool theSpeechRecognizerDefaultUseVad = true;

struct SileroVadParams {
	double speechProbabilityThreshold = 0.5;   // probability threshold to decide that sound is speech
	double minNonSpeechDuration = 0.1;   // min duration of a non-speech segment
	double minSpeechDuration = 0.25;   // min duration of a speech segment
	double speechPad = 0.03;   // padding added before and after each speech segment
};
/*
	Default Silero-VAD parameters for the UI.
*/
constexpr conststring32 theVadDefaultThreshold = U"0.5";
constexpr conststring32 theVadDefaultMinNonSpeechDuration = U"0.1";
constexpr conststring32 theVadDefaultMinSpeechDuration = U"0.25";
constexpr conststring32 theVadDefaultSpeechPad = U"0.03";
constexpr conststring32 theVadDefaultNonSpeechLabel = U"non-speech";
constexpr conststring32 theVadDefaultSpeechLabel = U"speech";

struct DiarizationParams {
	integer numSpeakers = 0;
	integer minSpeakers = 0;
	integer maxSpeakers = 0;
	bool allowSpeakersOverlap = true;
	double clusterThreshold = 0.7045654963945799;
	double segmentationStep = 0.1;
};
/*
	Default diarization parameters for the UI.
*/
constexpr conststring32 theDiarizationNumSpeakers = U"0 (= auto)";
constexpr conststring32 theDiarizationMinSpeakers = U"0";
constexpr conststring32 theDiarizationMaxSpeakers = U"0 (= unlimited)";
constexpr bool theDiarizationAllowSpeakersOverlap = true;
constexpr conststring32 theDiarizationDefaultNonSpeechLabel = U"";
constexpr conststring32 theDiarizationDefaultSpeechLabel = U"speech";
constexpr conststring32 theDiarizationClusterThreshold = U"0.7045654963945799";
constexpr conststring32 theDiarizationSegmentationStep = U"0.1";


struct SpeechSegment {
	autostring32 text;
	double tmin;
	double tmax;
};

struct WhisperTranscription {
	SpeechSegment fullTranscription;
	autovector <SpeechSegment> words;
	autovector <SpeechSegment> sentences;
};

inline std::set <structSpeechRecognizer *> theLivingSpeechRecognizers;
#include "SpeechRecognizer_def.h"

/*
	Functions to access lists of models and languages from everywhere.
*/
constSTRVEC theCurrentSpeechRecognizerModelNames ();
constSTRVEC theSpeechRecognizerLanguageNames ();

/*
	Class SpeechRecognizer functions.
*/
autoSpeechRecognizer SpeechRecognizer_create (conststring32 modelName, conststring32 languageName);
WhisperTranscription SpeechRecognizer_recognize (constSpeechRecognizer me, constSound sound,
		bool useVad, SileroVadParams const& sileroVadParams);

/*
	Silero-VAD functions.
*/
autovector <SpeechSegment> doSileroVad (constSound sound, SileroVadParams const& sileroVadParams,
		conststring32 nonSpeechLabel, conststring32 speechLabel);

/*
	Diarization functions.
*/
/*
	Run speaker diarization on Sound. Return a vector of timelines, one timeline per speaker.
	Each timeline covers [sound->xmin, sound->xmax] as alternating silent/active intervals,
	labeled `nonSpeechLabel` and `speechLabel` respectively.
*/
autovector <autovector <SpeechSegment>> doDiarization (constSound sound, DiarizationParams const& diarizationParams,
		conststring32 nonSpeechLabel, conststring32 speechLabel);


/* End of file SpeechRecognizer.h */
#endif
