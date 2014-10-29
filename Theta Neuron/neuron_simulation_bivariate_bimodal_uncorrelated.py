#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2014 Lucas Lin
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
    theta_temp = [None]*(2*halfn)
    k1, order = derivs(t, theta)
    for i in range(0, 2*halfn):
        theta_temp[i] = theta[i] + k1[i]*h/2
    k2, order = derivs(t+h/2, theta_temp)
    for i in range(0, 2*halfn):
        theta_temp[i] = theta[i] + k2[i]*h/2
    k3, order = derivs(t+h/2, theta_temp)
    for i in range(0, 2*halfn):
        theta_temp[i] = theta[i] + k3[i]*h
    k4, order = derivs(t+h, theta_temp)
    for i in range(0, 2*halfn):
        theta[i] += ((k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6)*h
    return t + h, theta, order


# Definition of the Ermentrout-Kopell canonical model with
# additional synaptic current:
# θ' = (1-cos(θ))+(1+cos(θ))[η+Isyn]
# Isyn = k/N * ΣPn(θ)
# Pn(θ) = an(1-cos(θ))^n
def derivs(t, theta):
    thetadot = [0]*(2*halfn)
    order = [0, 0]
    couple_temp = [0]*(2*halfn)
    sum = 0
    for i in range(0, 2*halfn):
        couple_temp[i] = an*pow((1-cos(theta[i])), sharp)
        sum += couple_temp[i]
    for i in range(0, 2*halfn):
        if i < halfn:
            thetadot[i] = (1-cos(theta[i]))+(1+cos(theta[i]))*(eta1[i]+(k1[i]/(2*halfn))*(sum-couple_temp[i]))
        else:
            thetadot[i] = (1-cos(theta[i]))+(1+cos(theta[i]))*(eta2[halfn-i]+(k2[halfn-i]/(2*halfn))*(sum-couple_temp[i]))
        order[0] += cos(theta[i])
        order[1] += sin(theta[i])
    order[0] /= (2*halfn)
    order[1] /= (2*halfn)
    return thetadot, order


##########################
# Assign Constant Values #
##########################

halfn = 5000                          # number of equations/neurons
thetai = []                  # initial values of n dependent variables in an array
for i in range(2*halfn):
    thetai.append(random()*2*pi)
ti = 0                          # inital value of independent variable
tf = 30                        # final value of independent variable
dt = 0.01                        # calculation step size
tout = 0.01                      # output interval

#### Model Variables ####

eta01 = 3                   # center of Lorentzian distribution for excitability
eta02 = -5
k01 = -3                      # center of Lorentzian distribution for synaptic strength
k02 = 5
deltaeta1 = 1                  # degree of heterogeneity for excitability
deltak1 = 1                      # degree of heterogeneity for synaptic strength
deltaeta2 = 0.5
deltak2 = 0.5
sharp = 2                       # sharpness of the synaptic function
an_values = {2: 2/3,            # normalization constants, pre-calculated
             3: 2/5,
             5: 8/63,
             9: 128/12155,
             15: 2048/9694845}
an = an_values[sharp]            # normalization constant for specific sharpness
dist1 = cauchy(k01, deltak1)       # Lorentzian distribution
dist2 = cauchy(k02, deltak2)
k1 = dist1.rvs(size=halfn).tolist()   # initializing the k-value for all the neurons
k2 = dist2.rvs(size=halfn).tolist()
dist3 = cauchy(eta01, deltaeta1)
dist4 = cauchy(eta02, deltaeta2)
eta1 = dist3.rvs(size=halfn).tolist()
eta2 = dist4.rvs(size=halfn).tolist()
print('Initial Values Loaded')

foldername = raw_input('Enter Folder Name: ')
os.makedirs(foldername)


def main():
    t = ti                                          # t is the time variable
    theta = [None]*2*halfn                                # initializing a list with n entries for each neuron
    t_output = [t]                                  # t_output collects data to be plotted
    load = 0
    loadstep = tout/(tf-ti)
    print("Calculating...%3.2f%%" % load, end='\r')
    for i in range(0, 2*halfn):
        theta[i] = thetai[i]
#    theta_output = [theta[:]]                       # theta_output collects data to be plotted
#    order_output = [[sum([cos(x) for x in theta_output[0]])/n, sum([sin(y) for y in theta_output[0]])/n]]
    order_output = [[sum([cos(x) for x in theta])/(2*halfn), sum([sin(y) for y in theta])/(2*halfn)]]
    while t < tf:
        tend = t + tout
        if tend > tf:
            tend = tf
        h = dt
        t, theta, order = integrator(t, theta, h, tend)    # returns values at every tout interval
        t_output.append(t)
#        theta_output.append(theta[:])
        order_output.append(order[:])
        load += loadstep
        print("Calculating...%3.2f%%" % (load*100), end='\r')
    print("Calculating...Done!")
    write_files(t_output, order_output)
    display_results(t_output, order_output)
#    write_files(t_output, theta_output, order_output)
#    display_results(t_output, theta_output, order_output)


def write_files(t, order):
    output = open(os.path.join(foldername, 'values.txt'), 'w')
    output.write('halfn = ' + str(halfn) + '\n' +
                 'ti = ' + str(ti) + '\n' +
                 'tf = ' + str(tf) + '\n' +
                 'dt = ' + str(dt) + '\n' +
                 'tout = ' + str(tout) + '\n\n' +
                 'Model Variables\n' +
                 'eta1 = ' + str(eta01) + '\n' +
                 'k1 = ' + str(k01) + '\n' +
                 'eta2 = ' + str(eta02) + '\n' +
                 'k2 = ' + str(k02) + '\n' +
                 'deltaeta1 = ' + str(deltaeta1) + '\n' +
                 'deltak1 = ' + str(deltak1) + '\n' +
                 'deltaeta2 = ' + str(deltaeta2) + '\n' +
                 'deltak2 = ' + str(deltak2) + '\n' +
                 'sharp = ' + str(sharp) + '\n')
    output = open(os.path.join(foldername, 'time.txt'), 'w')
    output.write(str(t).strip('[]'))
#    output = open(os.path.join(foldername, 'theta.txt'), 'w')
#    output.write(str(theta).strip('[]'))
    output = open(os.path.join(foldername, 'order.txt'), 'w')
    output.write(str(order).strip('[]'))


def display_results(t, order, **kwargs):
    #print(t)
    #print(theta)
    theta = None
    for name, value in kwargs.items():
        if name == "theta":
            theta = value
    if theta is None:
        animateorder(t, order)
    else:
        animatepoints(t, order, theta=theta)


def animateorder(t, order):
    fig, ax2 = plt.subplots()
    ax2.set_title('Order Parameter Trajectory')
    pausetime = (t[1]-t[0])/1000
    for i in range(0, len(t)):
        if i == 0:
            unpackorder = [[order[0][0]], [order[0][1]]]
            orderpoints, = ax2.plot(unpackorder[0], unpackorder[1], color='b')
            ax2.set_ylim([-1, 1])
            ax2.set_xlim([-1, 1])
        else:
            unpackorder[0].append(order[i][0])
            unpackorder[1].append(order[i][1])
            orderpoints.set_data(unpackorder[0], unpackorder[1])
        plt.pause(pausetime)
    plt.show()


def animatepoints(t, order, theta):
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