import matplotlib
from matplotlib import colors
from .plotcluster import *
from .gene_enrich import *


class Colormaps():
    def __init__(self):
        self.colormaps = {
            "grey_red": colors.LinearSegmentedColormap.from_list("grey_red", ["lightgray", "red", "darkred"], N=256),
            "grey_green": colors.LinearSegmentedColormap.from_list("grey_green", ["lightgray", "limegreen", "forestgreen"], N=256),
            "grey_yellow": colors.LinearSegmentedColormap.from_list("grey_yellow", ["lightgray", "yellow", "gold"], N=256),
            "grey_violet": colors.LinearSegmentedColormap.from_list("grey_violet", ["lightgray", "mediumvioletred", "indigo"], N=256),
            "grey_blue": colors.LinearSegmentedColormap.from_list("grey_blue", ["lightgray", "cornflowerblue", "darkblue"], N=256),
            "magma": colors.LinearSegmentedColormap.from_list("magma", matplotlib.colormaps["magma"].colors+['white'], N=256),
            "magma_r": colors.LinearSegmentedColormap.from_list("magma_r", ['white'] + matplotlib.colormaps["magma_r"].colors, N=256)
        }

    def __getitem__(self, key):
        return self.colormaps.get(key)
    
    def __str__(self):
        return ", ".join(filter(lambda attr: not attr.startswith("__"), self.colormaps))

    def __repr__(self):
        return self.__str__()