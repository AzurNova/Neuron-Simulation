#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

__author__ = 'Lucas'

from math import cos, sin, pi
from scipy.stats import cauchy
import matplotlib.pyplot as plt
from random import random
import os


def integrator(t, theta, h, tend):
    while t < tend:
        if tend - t < h:
            h = tend - t
        t, theta, order = rk4(t, theta, h)      # calculates rk4 with step size of dt
    return t, theta, order


def rk4(t, theta, h):
    theta_temp = [None]*n
    k1, order = derivs(t, theta)
    for i in range(0, n):
        theta_temp[i] = theta[i] + k1[i]*h/2
    k2, order = derivs(t+h/2, theta_temp)
    for i in range(0, n):
        theta_temp[i] = theta[i] + k2[i]*h/2
    k3, order = derivs(t+h/2, theta_temp)
    for i in range(0, n):
        theta_temp[i] = theta[i] + k3[i]*h
    k4, order = derivs(t+h, theta_temp)
    for i in range(0, n):
        theta[i] += ((k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6)*h
    return t + h, theta, order


# Definition of the Ermentrout-Kopell canonical model with
# additional synaptic current:
# θ' = (1-cos(θ))+(1+cos(θ))[η+Isyn]
# Isyn = k/N * ΣPn(θ)
# Pn(θ) = an(1-cos(θ))^n
def derivs(t, theta):
    thetadot = [0]*n
    order = [0, 0]
    couple_temp = [0]*n
    sum = 0
    for i in range(0, n):
        couple_temp[i] = an*pow((1-cos(theta[i])), sharp)
        sum += couple_temp[i]
    for i in range(0, n):
        thetadot[i] = (1-cos(theta[i]))+(1+cos(theta[i]))*(eta[i]+(k[i]/n)*(sum-couple_temp[i]))
        order[0] += cos(theta[i])
        order[1] += sin(theta[i])
    order[0] /= n
    order[1] /= n
    return thetadot, order


##########################
# Assign Constant Values #
##########################

n = 10000                          # number of equations/neurons
thetai = []                  # initial values of n dependent variables in an array
for i in range(n):
    thetai.append(random()*2*pi)
ti = 0                          # inital value of independent variable
tf = 30                        # final value of independent variable
dt = 0.001                        # calculation step size
tout = 0.001                      # output interval

#### Model Variables ####

eta0 = -0.5                   # center of Lorentzian distribution for excitability
k0 = 0.5                      # center of Lorentzian distribution for synaptic strength
deltaeta = 0.5                  # degree of heterogeneity for excitability
deltak = 0.1                      # degree of heterogeneity for synaptic strength
sharp = 2                       # sharpness of the synaptic function
an_values = {2: 2/3,            # normalization constants, pre-calculated
             3: 2/5,
             5: 8/63,
             9: 128/12155,
             15: 2048/9694845}
an = an_values[sharp]           # normalization constant for specific sharpness
dist1 = cauchy(k0, deltak)        # Lorentzian distribution
k = dist1.rvs(size=n).tolist()   # initializing the k-value for all the neurons
dist2 = cauchy(eta0, deltaeta)
eta = dist2.rvs(size=n).tolist()
print('Initial Values Loaded')

foldername = raw_input('Enter Folder Name: ')
os.makedirs(foldername)


def main():
    t = ti                                          # t is the time variable
    theta = [None]*n                                # initializing a list with n entries for each neuron
    t_output = [t]                                  # t_output collects data to be plotted
    load = 0
    loadstep = tout/(tf-ti)
    print("Calculating...%3.2f%%" % load, end='\r')
    for i in range(0, n):
        theta[i] = thetai[i]
    theta_output = [theta[:]]                       # theta_output collects data to be plotted
    order_output = [[sum([cos(x) for x in theta_output[0]])/n, sum([sin(y) for y in theta_output[0]])/n]]
    while t < tf:
        tend = t + tout
        if tend > tf:
            tend = tf
        h = dt
        t, theta, order = integrator(t, theta, h, tend)    # returns values at every tout interval
        t_output.append(t)
        theta_output.append(theta[:])
        order_output.append(order[:])
        load += loadstep
        print("Calculating...%3.2f%%" % (load*100), end='\r')
    print("Calculating...Done!")
    write_files(t_output, theta_output, order_output)
    display_results(t_output, theta_output, order_output)


def write_files(t, theta, order):
    output = open(os.path.join(foldername, 'values.txt'), 'w')
    output.write('n = ' + str(n) + '\n' +
                 'ti = ' + str(ti) + '\n' +
                 'tf = ' + str(tf) + '\n' +
                 'dt = ' + str(dt) + '\n' +
                 'tout = ' + str(tout) + '\n\n' +
                 'Model Variables\n' +
                 'eta0 = ' + str(eta0) + '\n' +
                 'k0 = ' + str(k0) + '\n' +
                 'deltaeta = ' + str(deltaeta) + '\n' +
                 'deltak = ' + str(deltak) + '\n' +
                 'sharp = ' + str(sharp) + '\n')
    output = open(os.path.join(foldername, 'time.txt'), 'w')
    output.write(str(t).strip('[]'))
#    output = open(os.path.join(foldername, 'theta.txt'), 'w')
#    output.write(str(theta).strip('[]'))
    output = open(os.path.join(foldername, 'order.txt'), 'w')
    output.write(str(order).strip('[]'))


def display_results(t, theta, order):
    #print(t)
    #print(theta)
    animatepoints(t, theta, order)


def animatepoints(t, theta, order):
    fig, (ax, ax2) = plt.subplots(1, 2, subplot_kw=dict(polar=True))
    ax2 = plt.subplot(1, 2, 2, polar=False)
    ax.set_yticklabels([])
    ax.set_title('Individual Neuron Simulation')
    ax2.set_title('Order Parameter Trajectory')
    r = [0.98]*len(theta[0])
    pausetime = (t[1]-t[0])/1000
    for i in range(0, len(t)):
        if i == 0:
            points, = ax.plot(theta[i], r, color='r', marker='.', linestyle='None')
            ax.set_rmax(1.0)
            ax.grid = True
            unpackorder = [[order[0][0]], [order[0][1]]]
            orderpoints, = ax2.plot(unpackorder[0], unpackorder[1], color='b')
            ax2.set_ylim([-1, 1])
            ax2.set_xlim([-1, 1])
        else:
            points.set_data(theta[i], r)
            unpackorder[0].append(order[i][0])
            unpackorder[1].append(order[i][1])
            orderpoints.set_data(unpackorder[0], unpackorder[1])
#        print(unpackorder)
        plt.pause(pausetime)
    plt.show()
    print('Plotting Done.')


if __name__ == '__main__': main()