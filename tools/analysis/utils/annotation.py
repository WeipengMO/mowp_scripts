import re

def natural_keys(text):
    '''It splits a string into a list of strings and numbers, and then sorts the list by the numbers
    
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    
    Parameters
    ----------
    text
        the string to split
    
    Returns
    -------
        A list of integers and strings.
    
    '''

    return [ int(c) if c.isdigit() else c for c in re.split(r'(\d+)', text) ]
