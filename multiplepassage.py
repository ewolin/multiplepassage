#!/usr/bin/env python
# After a large (M7ish) quake, 
# obtain waveform data from a broadband station
# and plot a spectrogram (OK, actually a continuous wavelet transfomr)
# and optionally overlay theoretical Rayleigh wave arrival times
# to demonstrate multiple passages of surface waves
# (and make a pretty plot)
# written by Emily Wolin, 2021
#
# requires ObsPy:
# https://docs.obspy.org/
# https://github.com/obspy/obspy/wiki#installation

import obspy
from obspy.clients.fdsn import Client
from obspy.signal.tf_misfit import plot_tfr  
from obspy.geodetics import gps2dist_azimuth

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

#####################
# earthquake, station, and plotting parameters
# edit to suit your quake/station of interest
# easy math: 3600 s = 1 hr, 86400 s = 1 day

# Earthquake origin parameters
# example is for an M8.1 in Kermadec Islands
elat = 29.735 # earthquake latitude
elon = -177.282 # earthquake longitude
t0 = obspy.UTCDateTime('2021-03-04T19:28:32') # quake origin time
t_start = t0 - 0.5*86400 
t_end = t0 + 2*86400

# Station parameters
slat = 34.94591
slon = -106.4572
net = "IU"
sta = "ANMO"
loc = "00"
cha = "VHZ"

# Plot parameters
fmin = 1e-3 # minimum frequency for spectrum + spectrogram
fmax = 1e-2 # max frequency for spectrum + spectrogram
figtitle = "M8.1 Kermadec Earthquake, recorded on GSN IU.ANMO.00"
clim=1e3 # max value for color bar
ylim=(-1*clim,clim) # y limits for time series plot
plot_rayleigh = False # plot theoretical Rayleigh wave arrival times
#####################

# initialize client for requesting waveforms
client = Client("IRIS")

print('requesting data')
st = client.get_waveforms("IU", "ANMO", "00", "VHZ", t_start, t_end)
st.write('waveforms.mseed')
st = obspy.read('waveforms.mseed')


# Optional: decimate if a very-low-sample-rate stream is not available
#print('decimating')
#st.decimate(factor=5)
#st.decimate(factor=5)
#st.decimate(factor=2)
#st.decimate(factor=2)

# Optional: trim to start/end points
# useful if a web service does not return exactly the start/end requested
#st.select(channel='LHZ', location='00')
#st.trim(starttime=t1)
#st.trim(endtime=st[0].stats.starttime + 86400)

print(st)
print(st[0].stats.delta)

# Remove mean and apply a bandpass filter 
print('filtering')
st.detrend('demean')
st.filter('bandpass', freqmin=fmin, freqmax=fmax)

tr = st[0]

###############################
# Plot figure!

dt = tr.stats.starttime - t0 # time between start of window and quake origin

fig = plot_tfr(tr.data, dt=tr.stats.delta, fmin=fmin, fmax=fmax, w0=10,
         clim=1e3, cmap='magma', show=False)

# plot dispersion curves on top
T, u = np.loadtxt('t-vs-u.csv',  unpack=True, delimiter=',')
f = 1./T

# Optional: Plot predicted Rayleigh wave arrival times 
# in t-vs-u.csv, digitized from Oliver (1962) Fig 2.
# https://pubs.geoscienceworld.org/ssa/bssa/article/52/1/81/101314

x_m, meh, meh = gps2dist_azimuth(elat, elon, slat, slon)
x = x_m/1000.
earth_circ = 4e4 # circumference of Earth in km, roughly 

ax_tfr = fig.axes[1]

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

font_props = {'ha':'center',
              'va':'bottom',
              'fontsize':'large'}

if plot_rayleigh:
    for i in range(4):
        n_odd = 2*i+1
        c_odd = cycle[n_odd]
        t_arr = (x+i*earth_circ)/u - dt
        ax_tfr.plot(t_arr,f, label='R{}'.format(2*i+1), color=c_odd)
        t = ax_tfr.text(t_arr[0], 1.0e-2, 'R{}'.format(n_odd), color=c_odd, 
                        **font_props)
        print(2*i+1, x+i*earth_circ)
    
        n_even = 2*(i+1)
        c_even = cycle[n_even]
        t_arr_min = (earth_circ - x + i*earth_circ)/u - dt
        ax_tfr.plot(t_arr_min,f, ls=':', label='R{}'.format(n_even), 
                    color=c_even)
        t = ax_tfr.text(t_arr_min[0], 1e-2, 'R{}'.format(n_even), 
                        color=c_even, **font_props)
        print(2*(i+1), earth_circ - x + i*earth_circ, '--')#


# Add labels etc to axes
ax_seis = fig.axes[0]
ax_seis.set_ylim(-1.0e3,1.0e3)
ax_seis.set_xlabel('time (s) since {} UTC'.format(tr.stats.starttime.strftime('%Y-%m-%dT%H:%M')))


ax_freq = fig.axes[2]
ax_freq.set_ylabel('frequency (Hz)')
fig.set_size_inches(10,5)
fig.suptitle(figtitle)

plt.savefig('multiplepassage.png', dpi=300)
print('plot saved to multiplepassage.png')
