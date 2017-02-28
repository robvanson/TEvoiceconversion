########################################################################

VoiceConversion.praat

Change the input speech to resemble Tracheoesophageal speech.
Changes the Pitch (F0) and pitch movements, duration. Filtered noise
is added as well as filtered "bubble" sounds.
Increase the Jitter and Shimmer of a speech recording to the
number given. Cannot reduce Jitter or Shimmer.
Note that Jitter and Shimmer are ill-defined in anything but
sustained vowels.

Uses the To PointProcess (periodic, cc) to calculate the jitter
and To PointProcess (periodic, peaks): 60, 300, "yes", "yes"
to change the timing of the periods.

Periods are moved with Overlap-and-Add

Shimmer is adapted using additive noise over an intensity tier and
adapting each period individually. Periods are determined with the 
To PointProcess (periodic, peaks) pulses.

########################################################################

Copyright (C) 2016-2017 NKI-AVL, R. J. J. H. van Son
R.v.Son@nki.nl

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Full license text is available at:
http://www.gnu.org/licenses/gpl-3.0.html

########################################################################

Input parameters (<=0 means "do not change"):

Input file    A file name (with full path). If a Sound object is 
              selected, that will be used instead
Pitch         Average pitch of the new speech in Hz 
              [F'(t) = Fnew/Fold * F(t)]
Pitch SD      Standard deviation of the Pitch of the new speech in Hz 
              (compresses pitch movements)
              [SD'(t) = SDnew/SDold * (F(t) - Faverage) + Faverage]
Duration      Factor with which to multiply the duration
HNR           Signal to Noise ratio of new noise added to obtain the 
              HNR given
Bubbles       Rate of bubble sounds added (per second). Select random 
              bubbles from bubbles.wav&bubbles.TextGrid
Bubbles SNR   Signal to Noise ratio of bubble sounds added 
              (use bubbles.wav)
Jitter        New jitter in %
Shimmer       New Shimmer in %
Voicing floor Lowest level of sound still considered voiced, in dB 
              below the maximum

Help          Print this text and exit

Output:
The input sound converted according to the specifications

