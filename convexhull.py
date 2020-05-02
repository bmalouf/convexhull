import math
from functools import cmp_to_key


# define constants
X = 0
Y = 1
LEFT = 0
RIGHT = 1
LOWER = 0
UPPER = 1


# determines whether 3 points make a left, right, or no turn
def check_turn(p, q, r, direction=None):
    if direction == LEFT:
        return ((r[Y] - p[Y]) * (q[X] - p[X])) - ((q[Y] - p[Y]) * (r[X] - p[X])) > 0
    if direction == RIGHT:
        return ((r[Y] - p[Y]) * (q[X] - p[X])) - ((q[Y] - p[Y]) * (r[X] - p[X])) < 0
    # left turn > 0, right turn < 0, straight == 0
    return ((r[Y] - p[Y]) * (q[X] - p[X])) - ((q[Y] - p[Y]) * (r[X] - p[X]))


# calculates squared distance between 2 points
def distance(p, q):
    return (p[X] - q[X])**2 + (p[Y] - q[Y])**2


# uses a binary search to find right tangent of a point to a convex hull (used by chans algorithm)
def right_tangent(points: list, p, type):
    low = 0
    high = len(points)
    if type == LOWER:
        while low < high:
            mid = (low + high) // 2
            if mid == 0:
                mid_pred = True
            else:
                mid_pred = not check_turn(p, points[mid], points[mid - 1], RIGHT)
            if mid == len(points) - 1:
                mid_suc = True
            else:
                mid_suc = check_turn(p, points[mid], points[mid + 1], LEFT)
            if not mid_suc:
                low = mid + 1
            elif not mid_pred:
                return right_tangent_l(points, p, LOWER)
            else:
                return points[mid]
    else:
        while low < high:
            mid = (low + high) // 2
            if mid == 0:
                mid_suc = True
            else:
                mid_suc = check_turn(p, points[mid], points[mid - 1], LEFT)
            if mid == len(points) - 1:
                mid_pred = True
            else:
                mid_pred = not check_turn(p, points[mid], points[mid + 1], RIGHT)
            if not mid_suc:
                high = mid
            elif not mid_pred:
                return right_tangent_l(points, p, UPPER)
            else:
                return points[mid]


def right_tangent_l(points, p, type):
    if type == LOWER:
        for i in range(1, len(points) - 1):
            if not check_turn(p, points[i], points[i - 1], RIGHT):
                if check_turn(p, points[i], points[i + 1], LEFT):
                    return points[i]
    else:
        for i in range(1, len(points) - 1):
            if check_turn(p, points[i], points[i - 1], LEFT):
                if not check_turn(p, points[i], points[i + 1], RIGHT):
                    return points[i]
    m = slope(points[0], points[-1])
    b = points[0][Y] - (m * points[0][X])
    p_int = p[Y] - (m * p[X])
    if p_int < b:
        return points[-1]
    else:
        return points[0]


# aki-toussaint heuristic
def aki_toussaint(points):
    left = min(points, key=lambda p: p[X])
    right = max(points, key=lambda p: p[X])
    bottom = min(points, key=lambda p: p[Y])
    top = max(points, key=lambda p: p[Y])
    tl_slope = slope(top, left)
    tr_slope = slope(top, right)
    bl_slope = slope(bottom, left)
    br_slope = slope(bottom, right)
    tl_intercept = top[Y] - (tl_slope * top[X])
    tr_intercept = top[Y] - (tr_slope * top[X])
    bl_intercept = bottom[Y] - (bl_slope * bottom[X])
    br_intercept = bottom[Y] - (br_slope * bottom[X])
    outside = list()
    for point in points:
        tl = point[Y] - (tl_slope * point[X])
        tr = point[Y] - (tr_slope * point[X])
        bl = point[Y] - (bl_slope * point[X])
        br = point[Y] - (br_slope * point[X])
        if tl >= tl_intercept:
            outside.append(point)
        elif tr >= tr_intercept:
            outside.append(point)
        elif bl <= bl_intercept:
            outside.append(point)
        elif br <= br_intercept:
            outside.append(point)
    return outside


def slope(p, q):
    return (p[Y] - q[Y]) / (p[X] - q[X])


def graham_scan_incremental(points, split=False):
    points = sorted(points, key=lambda x: x[0])  # sort only by x-coord
    if len(points) < 3:
        if split:
            return points, points
        else:
            return points
    upper_hull = list()
    lower_hull = list()
    upper_index = 0  # index of first item of upper hull
    lower_index = 0  # index of first item of lower hull
    start_index = 1  # index
    # since we can't assume general position, if multiple points share a minimum x-coord, we must find uppermost and lowermost
    while start_index < len(points) and points[start_index][X] == points[start_index - 1][X]:
        if points[upper_index][Y] < points[start_index][Y]:
            upper_index = start_index
        if points[lower_index][Y] > points[start_index][Y]:
            lower_index = start_index
        start_index += 1
    upper_hull.append(points[upper_index])
    lower_hull.append(points[lower_index])
    if start_index < len(points):
        upper_hull.append(points[start_index])
        lower_hull.append(points[start_index])
        start_index += 1
    # construct upper and lower hulls simultaneously
    for i in range(start_index, len(points)):
        while len(upper_hull) >= 2 and not check_turn(upper_hull[-2], upper_hull[-1], points[i], RIGHT):
            upper_hull.pop()
        upper_hull.append(points[i])
        while len(lower_hull) >= 2 and not check_turn(lower_hull[-2], lower_hull[-1], points[i], LEFT):
            lower_hull.pop()
        lower_hull.append(points[i])
    if split:  # return upper and lower hulls separately
        return lower_hull, upper_hull
    # concatenate upper and lower hulls while removing duplicate endpoints
    if upper_hull[-1] == lower_hull[-1]:
        upper_hull.pop()
        # also check colinearity on right side
        if check_turn(upper_hull[-1], lower_hull[-2], lower_hull[-1]) == 0:
            lower_hull.pop()
    upper_hull.reverse()
    convex_hull = lower_hull + upper_hull
    if convex_hull[0] == convex_hull[-1]:
        convex_hull.pop()
    return convex_hull


def graham_presorted(points):
    if len(points) < 3:
        return points
    upper_hull = list()
    lower_hull = list()
    upper_index = 0  # index of first item of upper hull
    lower_index = 0  # index of first item of lower hull
    start_index = 1  # index
    # since we can't assume general position, if multiple points share a minimum x-coord, we must find uppermost and lowermost
    while start_index < len(points) and points[start_index][X] == points[start_index - 1][X]:
        if points[upper_index][Y] < points[start_index][Y]:
            upper_index = start_index
        if points[lower_index][Y] > points[start_index][Y]:
            lower_index = start_index
        start_index += 1
    upper_hull.append(points[upper_index])
    lower_hull.append(points[lower_index])
    if start_index < len(points):
        upper_hull.append(points[start_index])
        lower_hull.append(points[start_index])
        start_index += 1
    # construct upper and lower hulls simultaneously
    for i in range(start_index, len(points)):
        while len(upper_hull) >= 2 and not check_turn(upper_hull[-2], upper_hull[-1], points[i], RIGHT):
            upper_hull.pop()
        upper_hull.append(points[i])
        while len(lower_hull) >= 2 and not check_turn(lower_hull[-2], lower_hull[-1], points[i], LEFT):
            lower_hull.pop()
        lower_hull.append(points[i])
    # concatenate upper and lower hulls while removing duplicate endpoints
    if upper_hull[-1] == lower_hull[-1]:
        upper_hull.pop()
        # also check colinearity on right side
        if check_turn(upper_hull[-1], lower_hull[-2], lower_hull[-1]) == 0:
            lower_hull.pop()
    upper_hull.reverse()
    convex_hull = lower_hull + upper_hull
    if convex_hull[0] == convex_hull[-1]:
        convex_hull.pop()
    return convex_hull


def graham_scan(points):
    if len(points) < 3:
        return sorted(points)
    # find bottom-leftmost point
    ch = list()
    p0 = min(points)
    ch.append(points.pop(points.index(p0)))
    points = sorted(points, key=cmp_to_key(lambda q, r: check_turn(p0, r, q)))
    ch.append(points[0])
    # scan points by polar angle, need to compare distance of colinear points
    for i in range(1, len(points)):
        # check for points colinear with p0
        if check_turn(ch[0], ch[-1], points[i]) == 0:
            d1 = distance(ch[0], ch[-1])
            d2 = distance(ch[0], points[i])
            if d2 > d1:  # new point is further away, pop old
                ch.pop()
            else:
                continue
        # check for colinear points on ch
        elif check_turn(ch[-2], ch[-1], points[i]) == 0:
            d1 = distance(ch[-2], ch[-1])
            d2 = distance(ch[-2], points[i])
            if d2 > d1:  # new point is further away, pop old
                ch.pop()
            else:
                continue
        while len(ch) >= 2 and not check_turn(ch[-2], ch[-1], points[i], LEFT):
            ch.pop()
        ch.append(points[i])
    return ch


def jarvis_march(points):
    ch = list()
    p0 = min(points)
    ch.append(p0)
    start = points.index(p0)
    last = start
    done = False
    while not done:
        current = (last + 1) % len(points)
        for i in range(len(points)):
            if i == last:
                continue
            turn = check_turn(ch[-1], points[current], points[i])
            if turn < 0:  # right turn
                current = i
            elif turn == 0:  # straight, us point that is further away
                d1 = distance(ch[-1], points[current])
                d2 = distance(ch[-1], points[i])
                if d2 > d1:
                    current = i
        if current == start:
            done = True
        else:
            ch.append(points[current])
            last = current
    return ch


def chans_restricted(points, k, partial_hull=None, prune_points=False):
    subsets = list()
    ch = list()
    if partial_hull is None:
        p0 = min(points)
        ch.append(p0)
    else:
        ch += partial_hull
    s = math.ceil(len(points) / k)
    for i in range(s):
        beg = i * k
        end = beg + k
        if end == len(points):
            end = len(points) + 1
        subsets.append(graham_scan_incremental(points[beg:end], True))
    i = len(ch) - 1  # adjust i if we pass in partial_hull
    while True:
        candidates = list()
        for subhull in subsets:
            candidates.append(right_tangent(subhull[LOWER], ch[-1], LOWER))  # check for tangent in upper hull
            candidates.append(right_tangent(subhull[UPPER], ch[-1], UPPER))  # check for tangent in lower hull
        next_point = candidates[0]
        for j in range(1, len(candidates)):
            if check_turn(ch[-1], next_point, candidates[j], RIGHT):
                next_point = candidates[j]
            elif check_turn(ch[-1], next_point, candidates[j]) == 0:
                d1 = distance(ch[-1], next_point)
                d2 = distance(ch[-1], candidates[j])
                if d2 > d1:
                    next_point = candidates[j]
        if next_point == ch[0]:
            if prune_points:
                return ch, True, subsets
            else:
                return ch, True
        ch.append(next_point)
        i += 1
        if i > k:
            if prune_points:
                return ch, False, subsets
            else:
                return ch, False


def chans_hull(points, save_hull=False, prune_points=False):
    k = 2
    done = False
    ch = list()
    prev_points = points
    while not done:
        if prune_points:
            k = min(k**2, len(prev_points))
        else:
            k = min(k**2, len(points))
        if save_hull and len(ch) > 0:
            if prune_points:
                ch, done, subsets = chans_restricted(prev_points, k, partial_hull=ch, prune_points=True)
            else:
                ch, done = chans_restricted(points, k, partial_hull=ch)
        else:
            if prune_points:
                ch, done, subsets = chans_restricted(points, k, prune_points=True)
            else:
                ch, done = chans_restricted(points, k)
        if prune_points and not done:
            prev_points = set()
            for subset in subsets:
                prev_points.update(subset[LOWER])
                prev_points.update(subset[UPPER])
            prev_points = list(prev_points)
    return ch


def divcon_hull(points, k=40):
    points = sorted(points, key=lambda x: x[0])
    lower_hull, upper_hull = divcon_divide(points, k)
    # merge upper and lower hulls
    if upper_hull[-1] == lower_hull[-1]:
        upper_hull.pop()
    # check colinearity on right side
    while check_turn(upper_hull[-1], lower_hull[-2], lower_hull[-1]) == 0:
        lower_hull.pop()
    while check_turn(upper_hull[-2], upper_hull[-1], lower_hull[-1]) == 0:
        upper_hull.pop()
    # check colinearity on left side
    while len(upper_hull) > 1 and upper_hull[0][X] == upper_hull[1][X]:
        upper_hull.pop(0)
    while len(lower_hull) > 1 and lower_hull[0][X] == lower_hull[1][X]:
        lower_hull.pop(0)
    upper_hull.reverse()
    convex_hull = lower_hull + upper_hull
    if convex_hull[0] == convex_hull[-1]:
        convex_hull.pop()
    return convex_hull


def divcon_divide(points, k):
    if len(points) <= k:
        return graham_scan_incremental(points, True)
    m = len(points) // 2
    left_hull = divcon_divide(points[0:m], k)
    right_hull = divcon_divide(points[m:], k)
    return divcon_merge(left_hull, right_hull)


def divcon_merge(left, right):
    # find upper tangent
    left_upper = len(left[UPPER]) - 1
    right_upper = 0
    next_left = left_upper - 1
    prev_right = right_upper + 1
    while True:
        while left_upper > 0 and not check_turn(right[UPPER][right_upper], left[UPPER][left_upper], left[UPPER][next_left], LEFT):
            left_upper -= 1
            next_left -= 1
        while right_upper < len(right[UPPER]) - 1 and not check_turn(right[UPPER][prev_right], right[UPPER][right_upper], left[UPPER][left_upper], LEFT):
            right_upper += 1
            prev_right += 1
        if left_upper == 0 or check_turn(right[UPPER][right_upper], left[UPPER][left_upper], left[UPPER][next_left], LEFT):
            break
    upper = left[UPPER][0:left_upper + 1] + right[UPPER][right_upper:]
    # find lower tangent
    left_lower = len(left[LOWER]) - 1
    right_lower = 0
    prev_left = left_lower - 1
    next_right = right_lower + 1
    while True:
        while left_lower > 0 and not check_turn(left[LOWER][prev_left], left[LOWER][left_lower], right[LOWER][right_lower], LEFT):
            left_lower -= 1
            prev_left -= 1
        while right_lower < len(right[LOWER]) - 1 and not check_turn(left[LOWER][left_lower], right[LOWER][right_lower], right[LOWER][next_right], LEFT):
            right_lower += 1
            next_right += 1
        if left_lower == 0 or check_turn(left[LOWER][prev_left], left[LOWER][left_lower], right[LOWER][right_lower], LEFT):
            break
    lower = left[LOWER][0:left_lower + 1] + right[LOWER][right_lower:]
    return [lower, upper]


def ccw(points, i):
    i += 1
    if i == len(points):
        i = 0
    return i


def cw(points, i):
    i -= 1
    if i <= -1:
        i = len(points) - 1
    return i
