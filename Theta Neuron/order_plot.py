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

__author__ = 'Lucas'

import matplotlib.pyplot as plt
import os
import sys


def main():
    foldername = raw_input('Enter Folder Name: ')
    instant = raw_input('Instant? (y/n): ')
    if instant != "y" and instant != "n":
        print('Error: Not a valid answer!')
        sys.exit()
    input = open(os.path.join(foldername, 'time.txt'), 'r')
    t = [float(x) for x in input.read().split(', ')]
    input = open(os.path.join(foldername, 'order.txt'), 'r')
    order = [[float(x) for x in pair.split(', ')] for pair in input.read().split('], [')]
    animatepoints(t, order, instant)


def animatepoints(t, order, instant):
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
        if instant == 'y':
            continue
        plt.pause(pausetime)
    plt.show()

if __name__ == '__main__': main()
