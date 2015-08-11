"""
    Get size of population from width of node
"""
import math

def GetSize(width):
    if isinstance(width, str):
        width = float(width)
    return 100 ** (width/0.3)

def GetWidth(size):
    if isinstance(size, str):
        size = float(size)
    return math.log(size, 100) * 0.3


if __name__ == "__main__":
    print "size is:", GetSize(0.15)
    print "width is:", GetWidth(25)
    print "width is:", GetWidth(10)
    print "width is:", GetWidth(20)
    print "width is:", GetWidth(30)
    print "width is:", GetWidth(40)
    print "width is:", GetWidth(100)
