import points as pts
import convexhull as ch
from CGAL import CGAL_Convex_hull_2 as cgal
import time


# load point set
p1 = pts.load_points('samplepoints.txt')
p2 = p1.copy()
p3 = p1.copy()
p4 = pts.convert_to_CGAL(p1)
# benchmark algorthms
start = time.time()
ch1 = ch.graham_scan_incremental(p1)
end = time.time()
print('Graham Andrew:', end - start, 'seconds')
print(ch1, '\n')
start = time.time()
ch2 = ch.chans_hull(p2)
end = time.time()
print('Chans Algorithm:', end - start, 'seconds')
print(ch2, '\n')
start = time.time()
ch3 = ch.divcon_hull(p3)
end = time.time()
print('Divide and Conquer:', end - start, 'seconds')
print(ch3, '\n')
ch4 = list()
start = time.time()
cgal.ch_graham_andrew(p4, ch4)
end = time.time()
ch4 = pts.convert_to_python(ch4)
print('CGAL Graham Andrew:', end - start, 'seconds')
print(ch4)
