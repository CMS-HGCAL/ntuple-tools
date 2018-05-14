# tool to determine if a point (plane, x, y) is within the sensor pad according to the given detector geometry
# to be used for masking of the rechits according to different geometry scenarios
# based on work of O.Arzi (CERN summer student in 2017)
from __future__ import division
import matplotlib.path as mplPath
import re
import itertools
from pandas import read_csv
import numpy as np
from array import array

# default geometry file name
default_geometry_file_path = 'fullgeometry_.txt'

# plane/layer with sensor pads/cells and max/min radius (where area in between is fully covered with sensor pads)
class Plane:
    def __init__(self, name, max_rad, min_rad, cells, edges):
        self.name = name
        self.min_rad = min_rad
        self.max_rad = max_rad
        self.cells = cells
        self.edges = edges
    # method to check if a (x,y) point is within the cells/pads of the plane
    def contains(self, p_x, p_y):
        rad = (p_x ** 2 + p_y ** 2) ** 0.5
        if self.min_rad < rad < self.max_rad:
            return True
        else:
            return any([e.contains(p_x, p_y) for e in self.edges])
    # method to check if a rechit with (x,y) coordiantes is within the cells/pads of the plane
    def is_contained(self, rechit):
        p_x, p_y = 10 * rechit.x(), 10 * rechit.y()
        return self.contains(p_x, p_y)

# sensor pad/cell with polygon path (defined by coordinates) and tag (full/half/edge pad)
class Cell:
    def __init__(self, coordinates, full, large, edge):
        self.coordinates = coordinates
        self.full = full
        self.large = large
        self.edge = edge
    # method to check if a (x,y) point is within the cell/pad
    def contains(self, p_x, p_y):
        return self.coordinates.contains_point((p_x, p_y))

# utility class to handle the detector geometry description
class GeoUtil:
    def __init__(self, geometry_file_path = default_geometry_file_path):
        self.geometry = geometry_file_path
        self.planes = self.read_planes()
    # get the layer indecies with basic information
    def _get_plane_indices(self, geometry_file_path):
        indices = {line_no: line for line_no, line in enumerate(geometry_file_path) if re.search('[a-zA-Z]+', line.decode()) is None}
        return indices
    # pairs of iterables
    def _pairwise(self, iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = itertools.tee(iterable)
        next(b, None)
        return list(zip(a, b))
    # extract info about one plan from the geometry description, and return the plane (with necessary information)
    def _read_plane(self, plane_details, start, end):
        # read in the plane id, max and min radia
        plane_details = [float(i) for i in plane_details.split()]
        plane_id, max_rad, min_rad = int(plane_details[0]), plane_details[2], plane_details[3]
        # read in the sensor pad vertecies of the current plane
        cells, edges = [], []
        df = read_csv(self.geometry, sep=' ', skiprows=lambda x: x not in range(start + 1, end), nrows=end - start, names=['flag', 'num', 'x1', 'y1', 'x2', 'y2', 'x3', 'y3', 'x4', 'y4', 'x5', 'y5', 'x6', 'y6', 'x7', 'y7', ], header=None, index_col=False)
        # loop over the sensor pads
        for index, row in df.iterrows():
            # sensor pad flag
            full = True if 'F' in row.flag else False
            large = True if 'O' in row.flag else False
            edge = True if 'e' in row.flag else False
            # sensor pad coordinates
            coordinates = np.empty((0, 2))
            for i in range(1, row.num + 1):
                xi, yi = row[2 * i], row[2 * i + 1]
                coordinates = np.vstack((coordinates, (xi, yi)))
            # get polygon path defined by the sensor pad coordinates
            polyPath = mplPath.Path(coordinates)
            # save sensor pad/cell, with its properties
            curr_cel = Cell(polyPath, full, large, edge)
            cells.append(curr_cel)
            if edge:
                edges.append(curr_cel)
        # return plane, with its properties
        plane = Plane(plane_id, max_rad, min_rad, cells, edges)
        return plane
    # extracts the geometry information for all layers and returns the list of planes (with necessary information)
    def read_planes(self):
        # get the list of all layers with the coresponding verticies that define sensor hexagons
        with open(self.geometry, 'rb') as geometry_file:
            planes_indices = self._get_plane_indices(geometry_file) # dictionary of indecies with layer descriptions (indecies as keys)
            keys = sorted(planes_indices.keys())
            pairs = self._pairwise(keys)
            # read in individual layers from the geometry description
            planes = [self._read_plane(planes_indices[start], start, end) for start, end in pairs]
        # return the list
        print("Number of planes loaded: ",len(planes))
        return planes
    # check if a point (x, y) is contained in the polygon defined with hexagon vertecies for the chosen plane
    def point_inside_plane(self, plane_num, p_x, p_y):
        return self.planes[plane_num].contains(p_x, p_y)
    # check if a point (x, y) is contained in the given geometry
    def point_inside_geometry(self, p_x, p_y):
        return any([plane.contains(p_x, p_y) for plane in self.planes])
    # helper function to prepare the geometry file (strip away the spaces)
    def remove_spaces(self, geometry_file_path):
        readfile = open(geometry_file_path)
        writefile = open(default_geometry_file_path, "w+")
        r = lambda line: re.sub(' +', ' ', line)
        for line in readfile:
            writefile.write(r(line))
        readfile.close()
        writefile.close()

# helper method to print/plot cells (to be deprecated)
import ROOT
def print_to_pic(cells, outDir = "./", tag = "test_cells"):
    # supress info messages
    ROOT.gErrorIgnoreLevel = ROOT.kInfo+1
    # set default style values
    ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.05)
    # set image extensions
    imgTypes = ["pdf","png"]
    # create canvas
    canvas = ROOT.TCanvas(tag, tag, 500, 500)
    mg = ROOT.TMultiGraph()
    for cell in cells:
        x = array("d", cell.coordinates.vertices[:, 0] / 10)
        x.append(x[0])
        y = array("d", cell.coordinates.vertices[:, 1] / 10)
        y.append(y[0])
        gi = ROOT.TGraph(len(x), x, y)
        gi.SetMarkerStyle(1)
        mg.Add(gi, "AL")
    # draw multi graph and save image
    mg.Draw("A")
    for imgType in imgTypes:
        canvas.SaveAs("{}/{}.{}".format(outDir, tag, imgType))
    return mg

