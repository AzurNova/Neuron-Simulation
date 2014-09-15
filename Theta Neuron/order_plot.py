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