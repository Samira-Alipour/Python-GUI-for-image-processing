import numpy as np
import os, sys
import argparse
import logging
import time
import shutil
import pylab as pl

import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from matplotlib.widgets import Button

'''Draw polygon regions of interest (ROIs) in matplotlib images,
similar to Matlab's roipoly function.
See the file example.py for an application.
Created by Joerg Doepfert 2014 based on code posted by Daniel
Kornhauser.
'''




class roipoly:
    def __init__(self, fig=[], ax=[], roicolor='b'):
        if fig == []:
            fig = plt.gcf()

        if ax == []:
            ax = plt.gca()

        self.previous_point = []
        self.allxpoints = []
        self.allypoints = []
        self.start_point = []
        self.end_point = []
        self.line = None
        self.roicolor = roicolor
        self.fig = fig
        self.ax = ax
        # self.fig.canvas.draw()


        self.__ID1 = self.fig.canvas.mpl_connect(
            'motion_notify_event', self.__motion_notify_callback)
        self.__ID2 = self.fig.canvas.mpl_connect(
            'button_press_event', self.__button_press_callback)


        if sys.flags.interactive:
            plt.show(block=False)
        else:
            plt.show()



    def getMask(self, currentImage):
        ny, nx = np.shape(currentImage)
        poly_verts = [(self.allxpoints[0], self.allypoints[0])]
        for i in range(len(self.allxpoints) - 1, -1, -1):
            poly_verts.append((self.allxpoints[i], self.allypoints[i]))

        # Create vertex coordinates for each grid cell...
        # (<0,0> is at the top left of the grid in this system)
        x, y = np.meshgrid(np.arange(nx), np.arange(ny))
        x, y = x.flatten(), y.flatten()
        points = np.vstack((x, y)).T

        ROIpath = mplPath.Path(poly_verts)
        grid = ROIpath.contains_points(points).reshape((ny, nx))
        return grid

    def displayROI(self, **linekwargs):
        l = plt.Line2D(self.allxpoints +
                       [self.allxpoints[0]],
                       self.allypoints +
                       [self.allypoints[0]],
                       color=self.roicolor, **linekwargs)
        ax = plt.gca()
        ax.add_line(l)
        plt.draw()

    def displayMean(self, currentImage, **textkwargs):
        mask = self.getMask(currentImage)
        meanval = np.mean(np.extract(mask, currentImage))
        stdval = np.std(np.extract(mask, currentImage))
        string = "%.3f +- %.3f" % (meanval, stdval)
        plt.text(self.allxpoints[0], self.allypoints[0],
                 string, color=self.roicolor,
                 bbox=dict(facecolor='w', alpha=0.6), **textkwargs)

    def __motion_notify_callback(self, event):
        # print(event)

        if event.inaxes:
            ax = event.inaxes
            x, y = event.xdata, event.ydata
            if (event.button == None or event.button == 1) and self.line != None:  # Move line around
                self.line.set_data([self.previous_point[0], x],
                                   [self.previous_point[1], y])
                self.fig.canvas.draw()

    def __button_press_callback(self, event):
        print(event)

        if event.inaxes:
            x, y = event.xdata, event.ydata
            ax = event.inaxes
            if event.button == 1 and event.dblclick == False:  # If you press the left button, single click
                if self.line == None:  # if there is no line, create a line
                    self.line = plt.Line2D([x, x],
                                           [y, y],
                                           marker='o',
                                           color=self.roicolor)
                    self.start_point = [x, y]
                    self.previous_point = self.start_point
                    self.allxpoints = [x]
                    self.allypoints = [y]

                    ax.add_line(self.line)
                    self.fig.canvas.draw()
                    # add a segment
                else:  # if there is a line, create a segment
                    self.line = plt.Line2D([self.previous_point[0], x],
                                           [self.previous_point[1], y],
                                           marker='o', color=self.roicolor)
                    self.previous_point = [x, y]
                    self.allxpoints.append(x)
                    self.allypoints.append(y)

                    event.inaxes.add_line(self.line)
                    self.fig.canvas.draw()
            elif ((event.button == 1 and event.dblclick == True) or
                      (
                              event.button == 3 and event.dblclick == False)) and self.line != None:  # close the loop and disconnect
                self.fig.canvas.mpl_disconnect(self.__ID1)  # joerg
                self.fig.canvas.mpl_disconnect(self.__ID2)  # joerg

                self.line.set_data([self.previous_point[0],
                                    self.start_point[0]],
                                   [self.previous_point[1],
                                    self.start_point[1]])
                ax.add_line(self.line)
                self.fig.canvas.draw()
                self.line = None

                if sys.flags.interactive:
                    pass
                else:
                    # figure has to be closed so that code can continue
                    plt.close(self.fig)







class Mask:
    def __init__(self):
        self.im = []
        self.width = 0
        self.height = 0

    def argp(self):
        parser = argparse.ArgumentParser(description='Mask single files')
        parser.add_argument('-i', dest='Filename', help="Input filename with *.ras extension")
        parser.add_argument('-r', dest='Resolution', default='hr', help="Processing resolution; hr, lr, ...")
        parser.add_argument('-e', dest='Extension', default='.HC.ras', help="Extension of raster file")


        args = parser.parse_args()
        return(args)

    def getWidth(self, res):
        width, height = vx.getStats(res, filename='masterpar.xml')
        return width, height


    def SelectRoi(self, event):
        image = self.im
        MyROI = roipoly(fig=[], ax=self.ax, roicolor='r')
        # MyROI.displayROI()
        mask = MyROI.getMask(image)
        # print(np.shape(mask))

        # mask = np.invert(mask)
        # image = mask*image
        # image[image == 0.0] = np.nan
        # self.ax.imshow(image)
        # pl.show(block=True)



    def SaveRoi(self, event):
        print("dfd")

    def mask(self):
        args = self.argp()
        srcename = args.Filename
        ifgres = args.Resolution
        extension = args.Extension
        filename = srcename + extension
        image = vi.readRas(filename)[0]
        image = image.astype(float)
        image[image == 0.0] = np.nan
        self.im = image.astype('float')
        dstdir = srcename + ".org"
        shutil.copy(srcename, dstdir)

        self.width, self.height = self.getWidth(ifgres)

        self.fig, self.ax = plt.subplots()

        self.add_ax = self.fig.add_axes([0.85, 0.72, 0.1, 0.04])
        self.save_ax = self.fig.add_axes([0.85, 0.66, 0.1, 0.04])
        pl.ion()
        B1 = Button(self.add_ax, 'ROI', color='white', hovercolor='green')
        B2 = Button(self.save_ax, 'Save',color='white', hovercolor='green')

        # print(self.ax)
        self.ax.imshow(image)
        B1.on_clicked(self.SelectRoi)
        B2.on_clicked(self.SaveRoi)
        MyROI = roipoly(fig=self.fig, ax=self.ax, roicolor='r')

        plt.show(block=True)







        mask = MyROI.getMask(image)
        mask = np.invert(mask)
        image = mask*image
        image[image == 0.0] = np.nan
        pl.imshow(image)
        pl.show(block=True)

        data = vi.readBin(dstdir, self.width, 'float')
        data = data*mask
        vi.writeBin(srcename, data)

        cmd = 'raster', srcename , '-r' + ifgres
        vl.execute2(cmd)



        # cmd = 'raster', srcename, '-r', ifgres, "-c", str(0), str(18.84), "-e", ".HC.ras"
        # vl.execute2(cmd)
        #
        cmd = 'raster ' + srcename + ' -r ' + ifgres + " -c " + str(0) + " " + str(18.84) + " -e .HC.ras"
        vl.execute2(cmd, shell=True)










if __name__ == '__main__':
    try:
        Mask().mask()
    except Exception:
        logging.exception('')
        sys.exit(1)
