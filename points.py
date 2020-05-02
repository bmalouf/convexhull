import random
from CGAL.CGAL_Kernel import Point_2


def create_points(n, lower_bound=-100, upper_bound=100):
    if (upper_bound - lower_bound + 1)**2 < n:
        print('bounds too small for set')
        return None
    points = set()
    while len(points) < n:
        x = random.randint(lower_bound, upper_bound)
        y = random.randint(lower_bound, upper_bound)
        points.add((x, y))
    return list(points)


def create_points_circle(n, radius=100):
    if radius**2 * 4 < n:
        print('radius too small for set')
        return None
    points = set()
    while len(points) < n:
        x = random.randint(-radius, radius)
        y = random.randint(-radius, radius)
        if x**2 + y**2 < radius**2:
            points.add((x, y))
    return list(points)


def create_points_triangle(n, lower_bound=-100, upper_bound=100):
    if (upper_bound - lower_bound + 1)**2 < n:
        print('bounds too small for set')
        return None
    points = set()
    while len(points) < (n-3):
        x = random.randint(lower_bound, upper_bound)
        y = random.randint(lower_bound, upper_bound)
        points.add((x, y))
    points.add((3 * upper_bound, 3 * upper_bound))
    points.add((3 * upper_bound, -3 * upper_bound))
    points.add((-3 * upper_bound, 0))
    return list(points)


def save_points(points, filename):
    with open(filename, 'w') as f:
        for point in points:
            f.write(str(point[0]) + ' ' + str(point[1]) + '\n')


def load_points(filename):
    points = list()
    with open(filename, 'r') as f:
        for line in f:
            x, y = [int(i) for i in line.strip().split(' ')]
            points.append((x, y))
    return points


def convert_to_CGAL(points):
    cgal_points = list()
    for point in points:
        cgal_points.append(Point_2(point[0], point[1]))
    return cgal_points


def convert_to_python(cgal_points):
    points = list()
    for point in cgal_points:
        x, y = [int(i) for i in str(point).split(' ')]
        points.append((x, y))
    return points
