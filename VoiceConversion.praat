#! praat
# 

form Select audio
	sentence Recorded_audio ookhetweer.wav
	real Pitch 70 (= Hz, mean)
	real Pitch_SD 15 (= % of mean)
	real Duration 1.3 (= mult. factor)
	real HNR 5 (= dB SNR)
	real Bubbles 5 (= rate)
	real Bubbles_SNR 15 (= dB SNR)
	real Jitter 5 (= %)
	real Shimmer 10 (= %)
	positive Voicing_floor_(dB) 15 (= below maximum)
	boolean Help 0
endform

########################################################################
# 
# VoiceConversion.praat
#
# Change the input speech to resemble Tracheoesophageal speech.
# Changes the Pitch (F0) and pitch movements, duration. Filtered noise
# is added as well as filtered "bubble" sounds.
# Increase the Jitter and Shimmer of a speech recording to the
# number given. Cannot reduce Jitter or Shimmer.
# Note that Jitter and Shimmer are ill-defined in anything but
# sustained vowels.
# 
# Uses the To PointProcess (periodic, cc) to calculate the jitter
# and To PointProcess (periodic, peaks): 60, 300, "yes", "yes"
# to change the timing of the periods.
# 
# Periods are moved with Overlap-and-Add
#
# Shimmer is adapted using additive noise over an intensity tier and
# adapting each period individually. Periods are determined with the 
# To PointProcess (periodic, peaks) pulses.
#
########################################################################
#
# Copyright (C) 2016-2017 NKI-AVL, R. J. J. H. van Son
# R.v.Son@nki.nl
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Full license text is available at:
# http://www.gnu.org/licenses/gpl-3.0.html
#
########################################################################
# 
# Input parameters (<=0 means "do not change"):
#
# Input file    A file name (with full path). If a Sound object is selected, that will be used instead
# Pitch         Average pitch of the new speech in Hz [F'(t) = Fnew/Fold * F(t)]
# Pitch SD      Standard deviation of the Pitch of the new speech in Hz (compresses pitch movements)
#               [SD'(t) = SDnew/SDold * (F(t) - Faverage) + Faverage]
# Duration      Factor with which to multiply the duration
# HNR           Signal to Noise ratio of new noise added to obtain the HNR given
# Bubbles       Rate of bubble sounds added (per second). Select random bubbles from bubbles.wav&bubbles.TextGrid
# Bubbles SNR   Signal to Noise ratio of bubble sounds added (use bubbles.wav)
# Jitter        New jitter in %
# Shimmer       New Shimmer in %
# Voicing floor Lowest level of sound still considered voiced, in dB below the maximum
# 
# Help          Print this text and exit
# 
# Output:
# The input sound converted according to the specifications
#
# Print debugging information
debug = 1

#
# Output:
# A Praat Sound object with the transformed speech
#
# Example:
# praat VoiceConversion.praat Speech/Example1.wav 70 15 1.3 5 5 15 5 10 15 no
#
# The Help text
#
if help
	clearinfo
	printline Help text
	printline
	printline Input parameters (<=0 means "do not change"):
	printline Input file'tab$''tab$'A file name (with full path). If a Sound object is selected, that will be used instead
	printline Pitch'tab$''tab$''tab$'Average pitch of the new speech in Hz [F'(t) = Fnew/Fold * F(t)]
	printline Pitch SD'tab$''tab$'Standard deviation of the Pitch of the new speech in Hz (compresses pitch movements)
	printline 'tab$''tab$''tab$''tab$'[SD'(t) = SDnew/SDold * (F(t) - Faverage) + Faverage]
	printline Duration'tab$''tab$'Factor with which to multiply the duration
	printline HNR'tab$''tab$''tab$''tab$'Signal to Noise ratio of new noise added to obtain the HNR given
	printline Bubbles'tab$''tab$''tab$'Rate of bubble sounds added (per second). Select random bubbles from bubbles.wav&bubbles.TextGrid
	printline Bubbles SNR'tab$''tab$'Signal to Noise ratio of bubble sounds added (use bubbles.wav)
	printline Jitter'tab$''tab$''tab$'New jitter in %
	printline Shimmer'tab$''tab$''tab$'New Shimmer in %
	printline Voicing floor'tab$'Lowest level of sound still considered voiced, in dB below the maximum
	printline
	printline Help'tab$''tab$''tab$'Print this text and exit
	printline
	printline Output:
	printline The input sound converted according to the specifications
	exit
endif
#
#
#

if numberOfSelected("Sound") > 0
	.recordedSound = selected("Sound")
elsif recorded_audio$ <> "" and fileReadable(recorded_audio$) 
	.recordedSound = Read from file: recorded_audio$
	Rename: "RecordedSpeech"
endif

bubblesAudioName$ = "bubbles.wav"

.thresshold = -voicing_floor

# Scale intensity:
select .recordedSound
global.setIntensity = Get intensity (dB)

call convert_speechOverlapAndAdd  .recordedSound .thresshold jitter shimmer pitch pitch_SD duration hNR bubbles bubbles_SNR

# Definitions of functions

# Main functions
procedure convert_speechOverlapAndAdd .recordedSound .thresshold .jitter .shimmer .pitch .pitch_SD .durationFactor .newHNR .bubble_rate .bubble_snr
	call change_pitch_duration .recordedSound .pitch .pitch_SD .durationFactor
	.newPitchSound = selected("Sound")
	
	call extractVoicingParameters .newPitchSound .thresshold
	.recordedTextGrid = selected("TextGrid")
	.recordedPulses = selected("PointProcess")
	.recordedInt = selected("Intensity")
	
	select .newPitchSound
	.recordedPulsesPeaks = To PointProcess (periodic, peaks): 60, 300, "yes", "yes"

	call create_additive_noise .newPitchSound .newHNR .recordedTextGrid
	.additiveNoise = selected("Sound")
	
	call add_bubbles .newPitchSound .bubble_rate .bubble_snr .recordedTextGrid 'bubblesAudioName$'
	.additiveBubbles = selected("Sound")
 	
	# Change Jitter, use CC to determine jitter and Peaks to change the periods
	selectObject: .recordedPulsesPeaks
	.newPointProcess = Copy: "New_Pulses"
	call set_jitter .jitter .newPointProcess .recordedPulses
	call test_overlap_add .newPitchSound .recordedPulsesPeaks .recordedTextGrid .newPointProcess .shimmer
	.newSound = selected("Sound")
	
	# Debug tests
	if debug
		# Old numbers
		selectObject: .recordedPulses
		.old_jitter = Get jitter (local): 0, 0, 0.0001, 0.03, 2
		selectObject: .newPointProcess
		.newPPjitter = Get jitter (local): 0, 0, 0.0001, 0.03, 2
		selectObject: .recordedSound
		plus .recordedPulses
		.old_amplitude = To AmplitudeTier (period): 0, 0, 0.0001, 0.03, 2
		.old_shimmer = Get shimmer (local): 0.0001, 0.03, 2
		
		
		selectObject: .newSound
		.pointP = To PointProcess (periodic, cc): 60, 300
		.new_jitter = Get jitter (local): 0, 0, 0.0001, 0.03, 2

		selectObject: .newSound
		plus .pointP
		.new_amplitude = To AmplitudeTier (period): 0, 0, 0.0001, 0.03, 2
		.new_shimmer = Get shimmer (local): 0.0001, 0.03, 2

		appendInfoLine:  "New Jitter: '.new_jitter:1%' ('.old_jitter:1%' ~> '.newPPjitter:1%')"
		appendInfoLine:  "New Shimmer: '.new_shimmer:1%' ('.old_shimmer:1%')"
		
		selectObject: .old_amplitude, .pointP, .new_amplitude
		Remove
	endif
	
	# Add noise to result
	call add_sounds .newSound .additiveNoise
	.resultNoise = selected("Sound")
	Rename: "NewSpeech"
	
	# Add bubbles to result
	call add_sounds .resultNoise .additiveBubbles
	.result = selected("Sound")
	Rename: "NewSpeech"
	
	# Clean up
	selectObject: .newPitchSound, .recordedTextGrid, .recordedPulses, .recordedInt, .newPointProcess, .recordedPulsesPeaks, .newSound, .additiveNoise, .additiveBubbles, .resultNoise
	Remove
	
	selectObject: .result
endproc

procedure change_pitch_duration .sound .pitch .pitchFraction .durationFactor
	select .sound
	.duration = Get total duration
	
	.manipulation = To Manipulation: 0.01, 70, 300
	.pitchTier = Extract pitch tier
	.currentPitch = Get mean (points): 0, 0
	.pitch_SD = .pitchFraction / 100 * .pitch
	
	select .manipulation
	.durationTier = Extract duration tier
	
	# Change duration
	if .durationFactor > 0
		.numPoints = Get number of points
		if .numPoints <= 0
			Add point: 0, 1
		endif
		Formula: "self*'.durationFactor'"
		select .manipulation
		plus .durationTier
		Replace duration tier
	endif
	
	if .pitch > 0
		select .pitchTier
		.factor = (.pitch / .currentPitch)
		Multiply frequencies: 0, .duration, .factor
		.currentSD = Get standard deviation (points): 0, 0
		
		if .currentSD > 0
			.factor = .pitch_SD / .currentSD
			Formula: "'.pitch' + (self - '.pitch') * '.factor'"
		endif

		select .manipulation
		plus .pitchTier
		Replace pitch tier
	endif
	
	.newSound = -1
	if .currentPitch > 0 or .durationFactor > 0
		select .manipulation
		.newSound = Get resynthesis (overlap-add)
	else
		select .sound 
		.newSound = Copy: "New Sound"
	endif
	select .manipulation
	plus .pitchTier
	plus .durationTier
	Remove
	
	select .newSound
endproc
	
procedure extractVoicingParameters .recordedSound .thresshold
	# The lowest level of voiced sounds
	select .recordedSound
	.pointPcc = To PointProcess (periodic, cc): 60, 300
	Rename: "RecordedPulses"
	.textGrid = To TextGrid (vuv): 0.02, 0.01
	Rename: "RecordedVoicing"
	.numIntervals = Get number of intervals: 1

	# Correct voicing boundaries
	select .recordedSound
	.intensity = To Intensity: 100, 0, "yes"
	Rename: "RecordedIntensity"
	.silences = To TextGrid (silences): .thresshold, 0.1, 0.05, "silent", "sounding"

	# Start boundaries
	for .i to .numIntervals
		select .textGrid
		.label$ = Get label of interval: 1, .i
		if .label$ = "V"
			.start = Get starting point: 1, .i
			.end = Get end point: 1, .i
			
			# Starting point of voiced interval
			select .silences
			.s = Get interval at time: 1, .start
			.sLabel$ = Get label of interval: 1, .s
			if .sLabel$ = "silent"
				.sStart = Get starting point: 1, .s
				.sEnd = Get end point: 1, .s
				select .textGrid
				if .sEnd < .end
					Set interval text: 1, .i, "U"
					# Shift boundaries: Insert&Remove
					Insert boundary: 1, .sEnd
					Set interval text: 1, .i+1, "V"
					if .i > 1
						Set interval text: 1, .i, ""
						Remove left boundary: 1, .i
					endif
				else
					# Low intensity, unvoiced
					Set interval text: 1, .i, "U"
				endif
			endif
		endif
	endfor
	
	# End boundaries
	for .i to .numIntervals
		select .textGrid
		.label$ = Get label of interval: 1, .i
		if .label$ = "V"
			.start = Get starting point: 1, .i
			.end = Get end point: 1, .i
			
			# Starting point of voiced interval
			select .silences
			.s = Get interval at time: 1, .end
			.sLabel$ = Get label of interval: 1, .s
			if .sLabel$ = "silent"
				.sStart = Get starting point: 1, .s
				.sEnd = Get end point: 1, .s
				select .textGrid
				if .sStart > .start
					Set interval text: 1, .i, "U"
					# Shift boundaries: Insert&Remove
					Insert boundary: 1, .sStart
					Set interval text: 1, .i, "V"
					if .i > 1
						Set interval text: 1, .i+1, ""
						Remove right boundary: 1, .i+1
					endif
				else
					# Low intensity, unvoiced
					Set interval text: 1, .i, "U"
				endif
			endif
		endif
	endfor
	
	select .silences
	Remove
	
	select .textGrid
	plus .pointPcc
	plus .intensity
endproc

procedure create_additive_noise .sound .newHNR .vuvTextGrid
	select .sound
	.duration = Get total duration
	.sampleFreq = Get sampling frequency

	.additiveNoiseSound = -1
	if .newHNR > 0
		# Determine noise level
		select .sound
		.originalIntensity = Get intensity (dB)
		.additiveNoise = .originalIntensity - .newHNR
		
		# Get filter
		select .sound
		.downsampled = Resample: 10000, 50
		.lpc = To LPC (autocorrelation): 10, 0.025, 0.005, 50
		plus .downsampled
		.source = Filter (inverse)
		.sourceInt = To Intensity: 70, 0, "yes"
		.sourceIntTier = To IntensityTier (peaks)

		# Create additive noise
		if .additiveNoise > 0
			.noise = Create Sound from formula: "WhiteNoise", 1, 0, .duration, .sampleFreq, "randomGauss(0,0.1)"
			plus .lpc
			.filteredNoise = Filter: "no"
			plus .sourceIntTier
			.additiveNoiseSoundTMP = Multiply: "yes"
			call set_VUV_to_zero .additiveNoiseSoundTMP .vuvTextGrid U
			Scale intensity: .additiveNoise
			.additiveNoiseSound = Resample: .sampleFreq, 50
			
			selectObject: .noise, .filteredNoise, .additiveNoiseSoundTMP
			Remove
		endif
		
		selectObject: .downsampled, .lpc, .source, .sourceInt, .sourceIntTier
		Remove
		
	endif
	
	if .additiveNoiseSound <= 0
		.additiveNoiseSound = Create Sound from formula: "AdditiveNoise", 1, 0, .duration, .sampleFreq, "0"
	endif
	
	select .additiveNoiseSound
endproc

procedure add_sounds .sound1 .sound2
	selectObject: .sound1, .sound2
	.stereo = Combine to stereo
	.addedSound = Convert to mono
	
	select .stereo
	Remove
	
	selectObject: .addedSound
endproc

# 
# Set Jitter to a specified number
# 
# Ti = ti - ti-1 (interval i)
# Jitter (absolute)  is Sum[ abs(Ti - Ti-1) ] / N-1
# Jitter = Jitter (absolute) / mean(Ti)
# 
# For a Normal distribution
# E(|X|) = sqrt(2/pi) * stdev(X)
# 
# E(Ti - Ti-1) = 0
# E(Ti^2) = var(Ti) + E(Ti)^2
# E(Ti*Ti-1) = cor(Ti, Ti-1) + E(Ti)^2
# var(Ti - Ti-1) = E(Ti^2 - 2*Ti*Ti-1 + Ti-1^2)
#                = 2*E(Ti^2) - 2*E(Ti*Ti-1)
#                = 2*[ var(Ti) * (1 - cor(Ti, Ti-1)) ]
# 
# Combine, assuming a Normal distribution:
# Jitter = E(|Ti - Ti-1|) / E(Ti)
#        = sqrt(2/pi) * stdev(Ti - Ti-1) / mean(Ti)
#        = sqrt(2/pi * var(Ti - Ti-1)) / mean(Ti)
#        = sqrt[ 4/pi * ( var(Ti) * (1 - cor(Ti, Ti-1)) ) ] / mean(Ti)
# 
# Change Ti -> T'i; Jitter -> a*Jitter while keeping mean(Ti) = mean(T'i) constant
# ei = (ti + ti-2)/2 - ti-1
# Jitter' = a * Jitter
#         = a * sqrt[ 4/pi * var(Ti - Ti-1) ] / mean(Ti)
# 
# => var(T'i - T'i-1) = a^2 * var(Ti - Ti-1) 
#                     = a^2 * E[ (Ti - Ti-1)^2 ]
#                     = a^2 * E[ (ti - ti-1 - ti-1 + ti-2)^2 ]
#                     = a^2 * 2 * E[ ((ti + ti-2)/2 - ti-1)^2 ]
#                     = a^2 * 2 * E[ ei^2 ]
#                     = 2 * E[ (a*ei)^2 ]
#                     = 2 * var(ei')
# 
# Generalizing, var(T'i - T'i-1) = 2*(var(ti-1) + var(ti) + var(ti+1))
# To increase Jitter -> Jitter'
# 1) Determine var(Ti - Ti-1) = (Jitter * mean(Ti))^2 * pi / 2
# 2) Calculate var(T'i - T'i-1) = (Jitter' * mean(T'i))^2 * pi / 2
# 3) Determine var to add: 
#    add_var(Ti - Ti-1) = var(T'i - T'i-1) - var(Ti - Ti-1)
# 4) Var of Noise to add per ti: add_var(ti) = add_var(Ti - Ti-1)/(2*3)
# 5) Sd of Noise to add per ti: add_sd(ti) = sqrt(add_var(ti))
# 
# .newJitter is in %
# Converts .pulses into pulses with new Jitter
procedure set_jitter .newJitter .pulses .pulsesCC
	
	if .pulses > 0 and .newJitter > 0
		.newJitter /= 100
		select .pulses
		# Use CC to determine real jitter
		if .pulsesCC > 0
			select .pulsesCC
		endif
		.current_jitter = Get jitter (local): 0, 0, 0.0001, 0.03, 2
		.current_abs_jitter = Get jitter (local, absolute): 0, 0, 0.0001, 0.03, 2
		.current_mean_period = Get mean period: 0, 0, 0.0001, 0.03, 2
		.current_stdev_period = Get stdev period: 0, 0, 0.0001, 0.03, 2

		if .newJitter > .current_jitter
			.current_var = .current_abs_jitter**2 * pi/2
			.end_var = (.newJitter * .current_mean_period)**2 * pi/2
			# The variance to add per boundary (total / (2*3))
			.add_var = (.end_var - .current_var) / 6
			.stdev_e = sqrt(.add_var)
			
			# Keep the original pulses just is case the order of the pulses might change
			select .pulses
			.origPulses = Copy: "Original_Pulses"
			.numPoints = Get number of points
			
			# New jitter
			# Change jitter by moving the ti according to 
			# t'i = ti - randomGauss(0, stdev(e'))
			for .p from 1 to .numPoints
				select .origPulses
				.t = Get time from index: .p
				.new_t = .t - randomGauss(0, .stdev_e)
		
				# Remove current point 
				select .pulses
				.r = Get nearest index: .t
				Remove point: .r
				Add point: .new_t
			endfor
			
			select .origPulses
			Remove
		else
			pause New jitter: '.newJitter' must be larger than current jitter '.current_jitter:4'
		endif
		
		# Calculate new jitter
		select .pulses
		.jitter_new = Get jitter (local): 0, 0, 0.0001, 0.03, 2
		.jitter_new *= 100
		.current_jitter *= 100
	endif
	
	select .pulses
endproc


# We cannot use the shimmer of a sentence, so we can only "add" shimmer
#
# .new_shimmer is in %
# .sound: Source Sound
# .pulses: PointProcess
# .voicing: VUV TextGrid
# .new_shimmer: New shimmer in %
procedure increase_shimmer .sound .pulses .voicing .newShimmer
	if .newShimmer > 0
		.newShimmer /= 100
		.shimmer_new = 0
		
		# Create Amplitude Tier and get current shimmer
		select .sound
		.duration = Get total duration
		plus .pulses
		.current_amplitude = To AmplitudeTier (period): 0, 0, 0.0001, 0.03, 2
		.current_shimmer = Get shimmer (local): 0.0001, 0.03, 2
		select .current_amplitude
		.numPoints = Get number of points
		.ampreal = Down to TableOfReal
		.sumamp = 0
		.n = 0
		for .p from 1 to .numPoints
			select .ampreal
			.tmp = Get value: .p, 2
			if .tmp > 0
				.sumamp += .tmp
				.n += 1
			endif
		endfor
		.meanAmp = .sumamp / .n
		
		# Sd must be multiplied with the amplitude
		if .newShimmer > .current_shimmer
			.new_var = (.newShimmer**2 - .current_shimmer**2) * .meanAmp**2 * pi / 2
		else
			.new_var = .newShimmer**2 * .meanAmp**2 * pi / 2
		endif
		if .new_var > 0	
			.new_sd = sqrt(.new_var / 2)
		else
			.new_sd = 0
		endif
		
		.new_amplitude = Create AmplitudeTier: "New_Amplitude", 0, .duration
		for .p from 1 to .numPoints
			select .ampreal
			.t = Get value: .p, 1
			.a = Get value: .p, 2
			if .a = undefined
				.a = 0
			endif
			if .a > 0
				.new_a = .a - randomGauss(0, .new_sd)
				if .new_a < 0 
					.new_a = 0
				endif
				
				# Add new value
				select .new_amplitude
				Add point: .t, .new_a / .a
			else
				Add point: .t, .a
			endif
		endfor
		
		# Set unvoiced parts to 1
		select .new_amplitude
		Add point: 0, 1
	
		select .voicing
		.numIntervals = Get number of intervals: 1
		for .i from 1 to .numIntervals
			select .voicing
			.t = Get end point: 1, .i
			select .new_amplitude
			Add point: .t, 1
		endfor
		
		select .ampreal
		plus .current_amplitude
		Remove
		
		# Overlay shimmer over sound
		select .sound
		plus .new_amplitude
		.new_sound = Multiply
		Rename: "NewSound_Shimmer"
		
		select .new_sound
		plus .pulses
		.shimmer_new = Get shimmer (local): 0, 0, 0.0001, 0.02, 1.3, 1.6
		.shimmer_new *= 100
		.current_shimmer *= 100
		.newShimmer *= 100
	
		select .new_amplitude
		Remove
		
	else
		select .sound
		.new_sound = Copy: "NewSound_Shimmer"
	endif
	
	select .new_sound
endproc

# Make a copy of the source to the target matching the pulses in source and target
# Copies fragments around pulses in sourcePulses under the direction of the 
# corresponding pulses in targetPulses using the Overlap&Add method (Gaussian window)
#
# Ignores voiceless parts, ie, intervals between pulses > .maxInt
# For voices, .maxInt should be ~0.02 (F0 > 50Hz). For other sounds, e.g., bubbles, this
# should be increased to fit the whole sound between pulses.
# 
# Midpoint between the pulses, periods add up to a factor of ~1.04. 
# At the pulses themselves, it adds up to ~1.12 (summed left and right)
# 
procedure overlap_add .sourceSound .sourcePulses .targetSound .targetPulses .maxInt
	# Maximum interval between pulses (maximum pitch period)
	if .maxInt <= 0
		.maxInt = 0.02
	endif
	.margin = 8*.maxInt
	select .sourceSound
	.sourceName$ = replace_regex$(selected$(), " ", "_", 0)
	select .targetSound
	.targetName$ = replace_regex$(selected$(), " ", "_", 0)

	# Iterate over target pulses
	select .targetPulses
	.numPulses = Get number of points
	for .p to .numPulses
		# Target
		select .targetPulses
		.tTarget = Get time from index: .p
		.pLeft = Get interval: .tTarget - 0.001
		.pRight = Get interval: .tTarget + 0.001
		# Source
		select .sourcePulses
		.q = Get nearest index: .tTarget
		.tSource = Get time from index: .q
		.qLeft = Get interval: .tSource - 0.001
		.qRight = Get interval: .tSource + 0.001
		# Gaussian window parameters (FWHM Left and Right)
		# FWHM = 2*sqrt(2*ln(2)) * c
		.cL = min(.pLeft,.qLeft)/(2*sqrt(2*ln(2)))
		.cR = min(.pRight,.qRight)/(2*sqrt(2*ln(2)))
		if not( .cL = undefined or .cL > .maxInt/(2*sqrt(2*ln(2))) or .cR = undefined or .cR > .maxInt/(2*sqrt(2*ln(2))) )
			# Copy one window
			select .targetSound
			Formula (part): .tTarget-.margin, .tTarget+.margin, 1, 1, "if x<.tTarget then self + '.sourceName$'((x - .tTarget) + .tSource)*exp(-1*(((x - .tTarget)/.cL)^2)/2) else self + '.sourceName$'((x - .tTarget) + .tSource)*exp(-1*(((x - .tTarget)/.cR)^2)/2) endif"
		endif
	endfor

endproc

# Test overlap_add
procedure test_overlap_add .sourceAudio .sourcePulses .vuvTextGrid .targetPulses .newShimmer
	# Use overlap-add to add new source intervals
	# Copy only voiced pulses
	call set_VUV_to_zero .targetPulses .vuvTextGrid U
	# Create a copy of the old source with voiced parts zeroed
	select .sourceAudio
	.testSource = Copy: "OaAsound"
	call set_VUV_to_zero .testSource .vuvTextGrid V

	# Copy the voiced parts of the new source to the zeroed voiced parts of the old source
	call overlap_add .sourceAudio .sourcePulses .testSource .targetPulses 0.02
	call increase_shimmer .testSource .targetPulses .vuvTextGrid .newShimmer
	.newSound = selected("Sound")
	Scale intensity: global.setIntensity

	select .testSource
	Remove
	
	select .newSound
	
endproc

# Set intervals matching a label text to Zero or remove the pulses
# Works on Sound and Pulses
procedure set_VUV_to_zero .sound .vuvTextGrid .zeroLabel$
	select .sound
	.objectType$ = selected$()
	.objectType$ = extractWord$ (.objectType$, "")
	select .vuvTextGrid
	.numIntervals = Get number of intervals: 1
	# Zero out VU intervals
	for .i to .numIntervals
		select .vuvTextGrid
		.vuvLabel$ = Get label of interval: 1, .i
		.start = Get starting point: 1, .i
		.end = Get end point: 1, .i
		if .vuvLabel$ = .zeroLabel$
			select .sound
			if .objectType$ = "Sound"
				Set part to zero: .start, .end, "at nearest zero crossing"
			elsif .objectType$ = "PointProcess"
				Remove points between: .start, .end
			else
				printline Unsupported object type for set_VUV_to_zero
			endif
		endif
	endfor
	select .sound
endproc


# 
# Add bubbles
#
#
# Add bubbles
# Select a random puls in the bubbles and add it to a random puls in the target
# 
# Creates a sound with only the bubbles
# 
procedure add_bubbles .sound .rate .snr .vuvTextGrid .bubblesAudioName$
	# Get filter
	select .sound
	.targetIntensity = Get intensity (dB)
	.targetDuration = Get total duration
	.tagetSamplingFrequency = Get sampling frequency
	.targetNumBubbles = .rate * .targetDuration
	.downsampled = Resample: 10000, 50
	.lpc = To LPC (autocorrelation): 10, 0.025, 0.005, 50
	plus .downsampled
	.source = Filter (inverse)
	.sourceInt = To Intensity: 70, 0, "yes"
	.sourceIntTier = To IntensityTier (peaks)
	select .sourceInt
	plus .downsampled
	plus .source
	Remove
	
	# Create an empty sound to receive the bubbles
	.bubblesAudio = Read from file: .bubblesAudioName$
	.bubblesTextGridName$ = replace_regex$(.bubblesAudioName$, "\.[a-z0-9]{2,}$", ".TextGrid", 0)
	.bubblesTextGrid = Read from file: .bubblesTextGridName$
	select .bubblesAudio
	.sourceName$ = replace_regex$(selected$(), " ", "_", 0)
	.bubblesSamplingFrequency = Get sampling frequency
	.bubblesIntensity = Get intensity (dB)
	.bubbleSound = Create Sound: "Bubbles", 0, .targetDuration, .bubblesSamplingFrequency, "0"
	
	# Fill the new Bubbles
	select .bubblesTextGrid
	.numIntervals = Get number of intervals: 1
	.bubblesFound = 0
	while .bubblesFound < .targetNumBubbles
		.i = randomInteger(1, .numIntervals)
		select .bubblesTextGrid
		.label$ = Get label of interval: 1, .i
		if .label$ = "sounding"
			.bubblesFound += 1
			.startPoint = Get starting point: 1, .i
			.endPoint = Get end point: 1, .i
			.midPoint = (.startPoint + .endPoint)/2
			.bubbleDuration = .endPoint - .startPoint
			
			# Get random insertion point
			.t = randomUniform (0.001, .targetDuration-0.001)
			.targetStart = .t - .bubbleDuration/2
			.targetEnd = .t + .bubbleDuration/2
			select .bubbleSound
			Formula (part): .targetStart, .targetEnd, 1, 1, "self + '.sourceName$'((x - .t) + .midPoint)"
		endif
	endwhile

	# Convert selected bubbles to scaled source
	select .bubbleSound
	.resampledBubbleSound = Resample: .tagetSamplingFrequency, 50
	plus .sourceIntTier
	.scaledBubbleSource = Multiply: "yes"
	call set_VUV_to_zero .scaledBubbleSource .vuvTextGrid U
	
	# The measured Intensity of the few selected bubbles can be too low. Correct for scaling
	select .scaledBubbleSource
	.bubbleSoundIntensity = Get intensity (dB)
	.attenuation = .bubblesIntensity - .bubbleSoundIntensity
	
	# Scale bubble sounds
	select .scaledBubbleSource
	Scale intensity: .targetIntensity - .snr - .attenuation
	
	select .scaledBubbleSource
	plus .lpc
	.filteredBubbles = Filter: "no"
	Rename: "FilteredBubbleNoise"
	.additiveBubblesSound = Resample: .tagetSamplingFrequency, 50

	# Clean up
	select .resampledBubbleSound
	plus .scaledBubbleSource
	plus .filteredBubbles
	plus .lpc
	plus .sourceIntTier
	plus .bubblesAudio
	plus .bubblesTextGrid
	plus .bubbleSound
	Remove

	select .additiveBubblesSound
endproc

procedure add_single_bubble .sourceAudio .sourcePulses .sourceI .targetAudio .targetPulses .targetI
	.margin = 1
	select .sourceAudio
	.sourceName$ = replace_regex$(selected$(), " ", "_", 0)
	select .targetAudio
	.targetName$ = replace_regex$(selected$(), " ", "_", 0)

	# Target
	select .targetPulses
	.tTarget = Get time from index: .targetI
	.pLeft = Get interval: .tTarget - 0.001
	.pRight = Get interval: .tTarget + 0.001
	
	# Source
	select .sourcePulses
	.tSource = Get time from index: .sourceI
	.qLeft = Get interval: .tSource - 0.001
	.qRight = Get interval: .tSource + 0.001
	
	# Gaussian window parameters (FWHM Left and Right)
	# FWHM = 2*sqrt(2*ln(2)) * c
	.c = (.qLeft+.qRight)/(2*sqrt(2*ln(2)))
	if not( .cL = undefined or .cR = undefined)
		# Copy one window
		select .targetAudio
		Formula (part): .tTarget-.margin, .tTarget+.margin, 1, 1, "self + '.sourceName$'((x - .tTarget) + .tSource)*exp(-1*(((x - .tTarget)/.c)^2)/2)"
	endif
endproc

