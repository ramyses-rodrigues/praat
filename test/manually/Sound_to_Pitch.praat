# test/manually/Sound_to_Pitch.praat
# Paul Boersma, 1 August 2025

sound1000 = Create Sound from formula: "sineWithNoise", 1, 0, 1000, 44100, "1/2 * sin(2*pi*377*x) + randomGauss(0,1)"
stopwatch
pitch1000 = noprogress To Pitch (raw autocorrelation): 0.001, 75, 600, 15, "no", 0.03, 0.45, 0.01, 0.35, 0.14
writeInfoLine: stopwatch
sound100 = Create Sound from formula: "sineWithNoise", 1, 0, 100, 44100, "1/2 * sin(2*pi*377*x) + randomGauss(0,1)"
stopwatch
pitch100 = noprogress To Pitch (raw autocorrelation): 0.001, 75, 600, 15, "no", 0.03, 0.45, 0.01, 0.35, 0.14
appendInfoLine: stopwatch
sound10 = Create Sound from formula: "sineWithNoise", 1, 0, 10, 44100, "1/2 * sin(2*pi*377*x) + randomGauss(0,1)"
stopwatch
pitch10 = noprogress To Pitch (raw autocorrelation): 0.001, 75, 600, 15, "no", 0.03, 0.45, 0.01, 0.35, 0.14
appendInfoLine: stopwatch

removeObject: sound10, sound100, sound1000, pitch10, pitch100, pitch1000

timeStep = 0.01 ; seconds
margin = 0.035 ; seconds
numberOfReplications = 100
maximumNumberOfFrames = 50
time# = zero# (maximumNumberOfFrames)
for replication to numberOfReplications
	for desiredNumberOfFrames from 1 to maximumNumberOfFrames
		durationOfFrames = desiredNumberOfFrames * timeStep
		durationOfSound = durationOfFrames + margin
		sound = Create Sound from formula: "sineWithNoise", 1, 0, durationOfSound, 44100, "1/2 * sin(2*pi*377*x) + randomGauss(0,1)"
		stopwatch
		pitch = noprogress To Pitch (raw autocorrelation): timeStep, 75, 600, 15, "no", 0.03, 0.45, 0.01, 0.35, 0.14
		time# [desiredNumberOfFrames] += stopwatch
		resultingNumberOfFrames = Get number of frames
		assert resultingNumberOfFrames = desiredNumberOfFrames   ; 'desiredNumberOfFrames' 'resultingNumberOfFrames'
		removeObject: sound, pitch
	endfor
endfor

for numberOfFrames from 1 to maximumNumberOfFrames
	appendInfoLine: numberOfFrames, " ", round (time# [numberOfFrames] / numberOfReplications * 1e6), "μ"
endfor

appendInfoLine: "Mean: ", round (mean (time#) / numberOfReplications * 1e6), "μ"