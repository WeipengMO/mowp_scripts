import textwrap
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import re

color_names = {
    'aliceblue', 'antiquewhite', 'aqua', 'aquamarine', 'azure', 'beige', 'bisque', 'black', 'blanchedalmond', 
    'blue', 'blueviolet', 'brown', 'burlywood', 'cadetblue', 'chartreuse', 'chocolate', 'coral', 'cornflowerblue', 
    'cornsilk', 'crimson', 'cyan', 'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray', 'darkgreen', 'darkgrey', 
    'darkkhaki', 'darkmagenta', 'darkolivegreen', 'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen', 
    'darkslateblue', 'darkslategray', 'darkslategrey', 'darkturquoise', 'darkviolet', 'deeppink', 'deepskyblue', 
    'dimgray', 'dimgrey', 'dodgerblue', 'firebrick', 'floralwhite', 'forestgreen', 'fuchsia', 'gainsboro', 'ghostwhite', 
    'gold', 'goldenrod', 'gray', 'green', 'greenyellow', 'grey', 'honeydew', 'hotpink', 'indianred', 'indigo', 'ivory', 
    'khaki', 'lavender', 'lavenderblush', 'lawngreen', 'lemonchiffon', 'lightblue', 'lightcoral', 'lightcyan', 
    'lightgoldenrodyellow', 'lightgray', 'lightgreen', 'lightgrey', 'lightpink', 'lightsalmon', 'lightseagreen', 
    'lightskyblue', 'lightslategray', 'lightslategrey', 'lightsteelblue', 'lightyellow', 'lime', 'limegreen', 'linen', 
    'magenta', 'maroon', 'mediumaquamarine', 'mediumblue', 'mediumorchid', 'mediumpurple', 'mediumseagreen', 
    'mediumslateblue', 'mediumspringgreen', 'mediumturquoise', 'mediumvioletred', 'midnightblue', 'mintcream', 
    'mistyrose', 'moccasin', 'navajowhite', 'navy', 'oldlace', 'olive', 'olivedrab', 'orange', 'orangered', 'orchid', 
    'palegoldenrod', 'palegreen', 'paleturquoise', 'palevioletred', 'papayawhip', 'peachpuff', 'peru', 'pink', 'plum', 
    'powderblue', 'purple', 'red', 'rosybrown', 'royalblue', 'saddlebrown', 'salmon', 'sandybrown', 'seagreen', 
    'seashell', 'sienna', 'silver', 'skyblue', 'slateblue', 'slategray', 'slategrey', 'snow', 'springgreen', 'steelblue', 
    'tan', 'teal', 'thistle', 'tomato', 'turquoise', 'violet', 'wheat', 'white', 'whitesmoke', 'yellow', 'yellowgreen'
}


def is_color_string(s: str) -> bool:
    assert isinstance(s, str), "Color must be a string"

    if re.match(r'^#(?:[0-9a-fA-F]{3}){1,2}$', s):
        return True
    
    return s.lower() in color_names


class ColorMaps():
    def __init__(self, *colors):
        color_list = []
        for _color in colors:
            if is_color_string(_color):
                color_list.append(_color)
            elif _color in mpl.colormaps:
                if isinstance(mpl.colormaps[_color], mpl.colors.ListedColormap):
                    _color = mpl.colormaps[_color].colors
                elif isinstance(mpl.colormaps[_color], mpl.colors.LinearSegmentedColormap):
                    _color = mpl.colormaps[_color](np.linspace(0, 1, 256))

                color_list.extend(list(_color))
            else:
                raise ValueError(f"Invalid color: {_color}")
        
        self.cmap = mpl.colors.LinearSegmentedColormap.from_list('cmap', color_list, N=256)

    def to_array(self):
        """
        Convert the color palette to a numpy array.

        Returns:
        --------
        numpy.ndarray
            The color palette as a numpy array.
        """
        return self.cmap(np.linspace(0, 1, 256))
    
    def palplot(self):
        """
        Plot the color palette.
        """
        gradient = np.linspace(0, 1, 256)
        gradient = np.vstack((gradient, gradient))
        fig, ax = plt.subplots(figsize=(6, .5))
        ax.imshow(gradient, aspect='auto', cmap=self.cmap)
        ax.set_axis_off()
        plt.show()
    
    def reverse(self):
        """
        Reverse the color palette.
        """
        self.palette = self.palette.reversed()

    
    def __repr__(self):
        desc = textwrap.dedent(f'''\
            ColorMaps object.

            Available colormaps:
                color_names: {', '.join(color_names)}
                Sequential: 'viridis', 'plasma', 'inferno', 'magma', 'cividis', 'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'
                Diverging: 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'
            
            More cmap can be found at:
                https://matplotlib.org/stable/tutorials/colors/colormaps.html
            ''')
        return desc
        

class ColorPalette():
    """
    A class representing a color palette.

    Parameters:
    -----------
    n_colors : int
        The number of colors in the palette.
    colors : str or list, optional
        The name of the colormap or a list of colors. 
        If not provided, the 'tab20c' colormap will be used.

    Attributes:
    -----------
    n_colors : int
        The number of colors in the palette.

    palette : matplotlib.colors.LinearSegmentedColormap
        The colormap representing the color palette.
    color_names : str
        The name of the colormap or 'custom' if a list of colors is provided.

    Methods:
    --------
    palplot()
        Plot the color palette.
    to_array()
        Convert the color palette to a numpy array.

    Examples:
    ---------
    Create a color palette with 10 colors using the 'tab20' colormap:
    >>> palette = ColorPalette(10, 'tab20')

    Create a color palette with 5 custom colors:
    >>> palette = ColorPalette(5, '#FF0000', '#00FF00', '#0000FF')
    """

    def __init__(self, n_colors: int = 5, *colors):
        self.n_colors = n_colors
        if len(colors) == 0:
            _colors = mpl.colormaps['tab20c'].colors
            self.palette = mpl.colors.LinearSegmentedColormap.from_list('color_palette', _colors, N=n_colors)
            self.color_names = 'tab20c'
        elif len(colors) == 1:
            _colors = mpl.colormaps[colors[0]]
            if isinstance(_colors, mpl.colors.ListedColormap):
                _colors = _colors.colors
            elif isinstance(_colors, mpl.colors.LinearSegmentedColormap):
                _colors = _colors(np.linspace(0, 1, n_colors))
            self.palette = mpl.colors.LinearSegmentedColormap.from_list('color_palette', _colors, N=n_colors)
            self.color_names = colors[0]
        else:
            self.palette = mpl.colors.LinearSegmentedColormap.from_list('color_palette', colors, N=n_colors)
            self.color_names = 'custom'
    
    def palplot(self):
        """
        Plot the color palette.
        """
        sns.palplot(self.to_array())
    
    def reverse(self):
        """
        Reverse the color palette.
        """
        self.palette = self.palette.reversed()
    
    def to_array(self):
        """
        Convert the color palette to a numpy array.

        Returns:
        --------
        numpy.ndarray
            The color palette as a numpy array.
        """
        return self.palette(np.linspace(0, 1, self.n_colors))
    
    def __repr__(self):
        """
        Return a string representation of the ColorPalette object.

        Returns:
        --------
        str
            A string representation of the ColorPalette object.
        """
        desc = textwrap.dedent(f'''\
            ColorPalette object with "{self.color_names}" palette, {self.n_colors} colors

            Available palettes:
                Diverging: 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'
                Qualitative: 'Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2', 'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b','tab20c'
            
            More palettes can be found at:
                https://matplotlib.org/stable/tutorials/colors/colormaps.html
            ''')
        return desc