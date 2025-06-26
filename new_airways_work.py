# current changes: n = 3.2, div = 1/1.5

from sympy import *
from scipy.spatial.transform import Rotation
from numpy.linalg import norm
import numpy
import pickle
from collections import defaultdict
from copy import deepcopy
import math

# model parameters
flow_rate_threshold = 0.00006
# flow_rate_threshold = 0.1
initial_length_to_diameter = 3
initial_rotation_angle = 90
min_distance_to_length = 3
max_distance_to_length = 6
volume_dividing_threshold = 0.05
increased_volume_dividing_threshold = 0.35
n = 2.8 # diameter exponent

savetofile = "test.pkl"

# supplementary rules
enable_4a = True
enable_6a = True
enable_7a = True
enable_8a = False

# for volume calculations
# min_grid_dimension = 6
# max_grid_dimension = 50
min_grid_dimension = 5
max_grid_dimension = 30

# evaluation
branch_count = 1
fail_flow_rate = 0
airway_volume = 0
flow_rates = [ ]
energy_loss = [ ]
terminal_volumes = [ ]
terminal_locations = [ ]

class Branch:
    def __init__(self, begin_point, end_point, flow_rate, diameter, length, direction_vector, branching_plane, grid_points, grid_scale, volume, depth):
        self.begin_point = begin_point
        self.end_point = end_point
        self.flow_rate = flow_rate
        self.diameter = diameter
        self.length = length
        self.direction_vector = direction_vector
        self.branching_plane = branching_plane
        self.grid_points = grid_points
        self.grid_scale = grid_scale
        self.volume = volume

        self.smaller_child = None
        self.larger_child = None

        self.depth = depth

segments_all = [ ]
segments_x = [ ]
segments_y = [ ] # ((y_begin, z_begin), (y_end, z_end), diameter)
diameters = [ ]

def get_segments_all(branch):
    global segments_all
    segments_all.append(((round(branch.begin_point.x, 3), round(branch.end_point.x, 3)), (round(branch.begin_point.y, 3), round(branch.end_point.y, 3)), (round(branch.begin_point.z, 3), round(branch.end_point.z, 3)), round(branch.diameter, 3)))
    if branch.smaller_child != None:
        get_segments_all(branch.smaller_child)
    if branch.larger_child != None:
        get_segments_all(branch.larger_child)

def get_segments_x(branch):
    global segments_x
    segments_x.append(((round(branch.begin_point.x, 3), round(branch.end_point.x, 3)), (round(branch.begin_point.z, 3), round(branch.end_point.z, 3)), round(branch.diameter, 3)))
    if branch.smaller_child != None:
        get_segments_x(branch.smaller_child)
    if branch.larger_child != None:
        get_segments_x(branch.larger_child)

def get_segments_y(branch):
    global segments_y
    segments_y.append(((round(branch.begin_point.y, 3), round(branch.end_point.y, 3)), (round(branch.begin_point.z, 3), round(branch.end_point.z, 3)), round(branch.diameter, 3)))
    if branch.smaller_child != None:
        get_segments_y(branch.smaller_child)
    if branch.larger_child != None:
        get_segments_y(branch.larger_child)

def generate_children(branch):
    global branch_count
    global fail_flow_rate
    global airway_volume
    global flow_rates
    global energy_loss
    global terminal_volumes
    global terminal_locations

    global n
    global count
    if branch.depth == 2:
        if count == 0: count += 1
        else: return
    # if branch.depth > 8:
    #     return

    if branch_count % 1000 == 0:
        global segments_x, segments_y, diameters, trunk
        segments_x, segments_y, segments_all, diameters = [ ], [ ], [ ], [ ]
        get_segments_x(trunk)
        get_segments_y(trunk)
        get_segments_all(trunk)
        segments = [segments_x, segments_y, segments_all]
        with open(savetofile, 'wb') as file:
            pickle.dump(segments, file)
        with open(savetofile, 'wb') as file:
            print([branch_count, fail_flow_rate])
            pickle.dump([branch_count, fail_flow_rate], file)

    # 9
    if branch.depth > 2 and not within_region(branch.grid_points, branch.end_point, branch.grid_scale):
        fail_flow_rate += branch.flow_rate
        branch = None # comment out to rmv gaps
        return
    if branch.flow_rate <= flow_rate_threshold:
        terminal_volumes.append(branch.volume)
        terminal_locations.append((branch.end_point[0], branch.end_point[1], branch.end_point[2]))
        return

    branch_count += 1
    if branch.depth > 18:
        print(branch_count, branch.depth+1)
    airway_volume += (branch.diameter/2)**2 * branch.length * pi
    flow_rates.append(branch.flow_rate)
    energy_loss.append(branch.length/(branch.diameter/2)**4*branch.flow_rate**2)

    # 8a
    rule_8a = True
    count_8a, max_r, rotation_angle = 0, 0, 0
    while rule_8a:
        # 4
        space_dividing_plane1 = Plane(branch.begin_point, branch.end_point, branch.end_point + branch.branching_plane.normal_vector) # create space-dividing plane using 3 points
        temp = space_dividing_plane1.normal_vector
        space_dividing_plane2 = Plane(branch.begin_point, normal_vector=(-temp[0], -temp[1], -temp[2]))

        # calculate volumes, volume-dividing ratio
        vol1 = get_volume(branch.grid_points, [space_dividing_plane1], False)
        vol2 = get_volume(branch.grid_points, [space_dividing_plane2], False)
        volumes = [vol1[1], vol2[1]]
        if volumes[0] <= 1 or volumes[1] <= 1:
            terminal_volumes.append(branch.volume)
            terminal_locations.append((branch.end_point[0], branch.end_point[1], branch.end_point[2]))
            return

        if (volumes[0] < volumes[1]):
            smaller_child_side = space_dividing_plane1.normal_vector
            smaller_points, larger_points = vol1[0], vol2[0]
            smaller_space_dividing_plane, larger_space_dividing_plane = space_dividing_plane1, space_dividing_plane2
            smaller_volume, larger_volume = vol1[1], vol2[1]
            r = volumes[0] / (volumes[0] + volumes[1])
        else:
            smaller_child_side = space_dividing_plane2.normal_vector
            smaller_points, larger_points = vol2[0], vol1[0]
            smaller_space_dividing_plane, larger_space_dividing_plane = space_dividing_plane2, space_dividing_plane1
            smaller_volume, larger_volume = vol2[1], vol1[1]
            r = volumes[1] / (volumes[0] + volumes[1])

        # eq 5 (rule 6)
        smaller_angle = (numpy.arccos((1+r**(4/n)-(1-r)**(4/n)) / (2*r**(2/n))))*180/numpy.pi
        larger_angle = (numpy.arccos((1+(1-r)**(4/n)-r**(4/n)) / (2*(1-r)**(2/n))))*180/numpy.pi
        # eq 3
        smaller_diameter = branch.diameter*r**(1/n)
        larger_diameter = branch.diameter*(1-r)**(1/n)
        # eq 2
        smaller_flow_rate = (smaller_diameter / trachea_diameter) ** n
        larger_flow_rate = (larger_diameter / trachea_diameter) ** n

        # 8a
        if enable_8a:
            if r > max_r:
                max_r = r
                max_r_angle = rotation_angle
            if r < 0.05 or (smaller_flow_rate < 1.5 * 0.00006 and r < 0.35):
                count_8a += 1
                if count_8a >= 40:
                    rotation_angle = max_r_angle
                    rule_8a = False
                else:
                    rotation_angle = float(9 * (-1)**count_8a * math.ceil(count_8a/2))
                temp = branch.branching_plane.normal_vector
                new_normal_vector = rotate_vector(branch.direction_vector, rotation_angle, (float(temp[0]), float(temp[1]), float(temp[2])))
                branch.branching_plane = Plane(branch.begin_point, normal_vector=new_normal_vector)
            else:
                rule_8a = False
        else:
            rule_8a = False 
    # end 8a loop


    # 4a
    if enable_4a and smaller_angle - larger_angle >= 10:
        # rotate the normal vector of the space-dividing plane
        angle_to_rotate = (smaller_angle + larger_angle) / 2 - larger_angle
        axis = branch.branching_plane.normal_vector
        new_normal_vector = rotate_vector(axis, angle_to_rotate, (float(smaller_child_side[0]), float(smaller_child_side[1]), float(smaller_child_side[2])))
        if dot(new_normal_vector, branch.direction_vector) < 0:
            new_normal_vector = rotate_vector(axis, -angle_to_rotate, (float(smaller_child_side[0]), float(smaller_child_side[1]), float(smaller_child_side[2])))
        # create new plane using normal vector
        smaller_half_dividing_plane = Plane(branch.end_point, normal_vector=new_normal_vector)
        larger_half_dividing_plane = Plane(branch.end_point, normal_vector=-new_normal_vector)
        # fix volumes and ratio
        vol1 = get_volume(branch.grid_points, [smaller_space_dividing_plane, smaller_half_dividing_plane], True)
        vol2 = get_volume(branch.grid_points, [larger_space_dividing_plane, larger_half_dividing_plane], False)
        smaller_points, larger_points = vol1[0], vol2[0]
        smaller_volume, larger_volume = vol1[1], vol2[1]
        if smaller_volume <= 1 or larger_volume <= 1:
            terminal_volumes.append(branch.volume)
            terminal_locations.append((branch.end_point[0], branch.end_point[1], branch.end_point[2]))
            return
        r = vol1[1] / (vol1[1]+vol2[1])
        smaller_angle = (numpy.arccos((1+r**(4/n)-(1-r)**(4/n)) / (2*r**(2/n))))*180/numpy.pi
        larger_angle = (numpy.arccos((1+(1-r)**(4/n)-r**(4/n)) / (2*(1-r)**(2/n))))*180/numpy.pi
    
    # fix
    if branch.depth == 0:
        smaller_angle += 10
        larger_angle += 10
    if branch.depth == 2:
        smaller_angle += 15
    if branch.depth == 3:
        larger_angle += 15
        smaller_angle += 15
    if branch.depth == 4:
        if branch.end_point.z < 15:
            smaller_angle -= 15
    if branch.depth >= 5:
        larger_angle += 10

    # 6
    # eq 3
    n = 3.2
    smaller_diameter = branch.diameter*r**(1/n)
    larger_diameter = branch.diameter*(1-r)**(1/n)
    # new
    origsumdiam = smaller_diameter + larger_diameter
    sumdiam = smaller_diameter**(1/1.5) + larger_diameter**(1/1.5)
    smaller_diameter = smaller_diameter**(1/1.5) / sumdiam * origsumdiam
    larger_diameter = larger_diameter**(1/1.5) / sumdiam * origsumdiam
    # eq 2
    smaller_flow_rate = (smaller_diameter / trachea_diameter) ** n
    larger_flow_rate = (larger_diameter / trachea_diameter) ** n
    n = 2.8

    temp1, temp2 = redefine_grid(smaller_points, branch.grid_scale), redefine_grid(larger_points, branch.grid_scale)
    smaller_grid_points, larger_grid_points = temp1[0], temp2[0]
    smaller_scale, larger_scale = temp1[1], temp2[1]

    # 6a (fix branching angle)
    if enable_6a:
        smaller_bisecting_angle = bisecting_angle(branch, smaller_grid_points, smaller_space_dividing_plane)
        larger_bisecting_angle = bisecting_angle(branch, larger_grid_points, larger_space_dividing_plane)
        if smaller_bisecting_angle > smaller_angle:
            smaller_angle = (smaller_angle + smaller_bisecting_angle) / 2
        if larger_bisecting_angle > larger_angle:
            larger_angle = (larger_angle + larger_bisecting_angle) / 2
    # rotate on plane for direction vector
    axis = branch.branching_plane.normal_vector
    smaller_direction_vector = rotate_vector(axis, smaller_angle, branch.direction_vector)
    larger_direction_vector = rotate_vector(axis, larger_angle, branch.direction_vector)
    # fix direction of rotation by testing which dot product is negative
    if dot(smaller_direction_vector, smaller_child_side) < 0:
        smaller_direction_vector = rotate_vector(axis, -smaller_angle, branch.direction_vector)
    else:
        larger_direction_vector = rotate_vector(axis, -larger_angle, branch.direction_vector)

    # 7 (lengths)
    smaller_length = smaller_diameter * initial_length_to_diameter
    larger_length = larger_diameter * initial_length_to_diameter

    # 7a
    if enable_7a:
        smaller_length = fix_length(branch, smaller_grid_points, smaller_length, smaller_diameter, smaller_direction_vector, smaller_scale)
        larger_length = fix_length(branch, larger_grid_points, larger_length, larger_diameter, larger_direction_vector, larger_scale)

    # 8 (rotation angle)
    if initial_rotation_angle == 90:
        smaller_branching_plane = branch.branching_plane.perpendicular_plane(branch.end_point, branch.end_point + smaller_direction_vector)
        larger_branching_plane = branch.branching_plane.perpendicular_plane(branch.end_point, branch.end_point + larger_direction_vector)
    else:
        temp = branch.branching_plane.normal_vector
        smaller_normal_vector = rotate_vector(smaller_direction_vector, initial_rotation_angle, (float(temp[0]), float(temp[1]), float(temp[2])))
        larger_normal_vector = rotate_vector(larger_direction_vector, initial_rotation_angle, (float(temp[0]), float(temp[1]), float(temp[2])))
        smaller_branching_plane = Plane(branch.end_point, normal_vector=smaller_normal_vector)
        larger_branching_plane = Plane(branch.end_point, normal_vector=larger_normal_vector)

    # endpoints
    smaller_end_point = branch.end_point + smaller_length * smaller_direction_vector
    larger_end_point = branch.end_point + larger_length * larger_direction_vector

    # save memory
    branch.grid_points = None
    branch.flow_rate = None
    branch.length = None
    branch.direction_vector = None
    branch.branching_plane = None
    branch.grid_scale = None
    # generate children
    smaller_child = Branch(branch.end_point, smaller_end_point, smaller_flow_rate, smaller_diameter, smaller_length, smaller_direction_vector, smaller_branching_plane, smaller_grid_points, smaller_scale, smaller_volume/smaller_scale**3, branch.depth+1)
    larger_child = Branch(branch.end_point, larger_end_point, larger_flow_rate, larger_diameter, larger_length, larger_direction_vector, larger_branching_plane, larger_grid_points, larger_scale, larger_volume/larger_scale**3, branch.depth+1)
    branch.smaller_child = smaller_child
    branch.larger_child = larger_child
    generate_children(smaller_child)
    generate_children(larger_child)

    if smaller_child == None and larger_child == None:
        terminal_volumes.append(branch.volume)
        terminal_locations.append((branch.end_point[0], branch.end_point[1], branch.end_point[2]))
    branch.volume = None

def dot(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

# takes the longest
def get_volume(points, planes, convex): # inclusive of the points on the plane
    new_points = deepcopy(points)
    volume = 0
    for x, v1 in new_points.items():
        for y, v2 in v1.items():
            for z, p in v2.items():
                if p:
                    flag = False
                    for pl in planes:
                        plane_point_vector = (x-pl.p1.x, y-pl.p1.y, z-pl.p1.z)
                        if dot(pl.normal_vector, plane_point_vector) < 0:
                            if convex:
                                new_points[x][y][z] = 0
                                break
                        else: 
                            flag = True
                    if not convex and not flag:
                        new_points[x][y][z] = 0
    for x, v1 in new_points.items():
        for y, v2 in v1.items():
            for z, p in v2.items():
                if p:
                    volume += 1
    return (new_points, volume)

def redefine_grid(points, scale):
    x_min, x_max, y_min, y_max, z_min, z_max = bounds_xyz(points)
    grid = auto_dict()
    for x in range(int(x_min*scale), int(x_max*scale+1)):
        for y in range(int(y_min*scale), int(y_max*scale+1)):
            for z in range(int(z_min*scale), int(z_max*scale+1)):
                grid[x/scale][y/scale][z/scale] = points[x/scale][y/scale][z/scale]
    while ((x_max-x_min)*scale < min_grid_dimension or (y_max-y_min)*scale < min_grid_dimension or (z_max-z_min)*scale < min_grid_dimension):
        if (x_max-x_min)*scale >= max_grid_dimension/2 or (y_max-y_min)*scale >= max_grid_dimension/2 or (z_max-z_min)*scale >= max_grid_dimension/2:
            break
        grid = auto_dict()
        for x in range(int(x_min*scale), int(x_max*scale+1)):
            for y in range(int(y_min*scale), int(y_max*scale+1)):
                for z in range(int(z_min*scale), int(z_max*scale+1)):
                    grid[x/scale][y/scale][z/scale] = points[x/scale][y/scale][z/scale]
        # adding intermediate points
        scale *= 2
        for i in range(int(x_min*scale)-1, int(x_max*scale+1)+1):
            for j in range(int(y_min*scale)-1, int(y_max*scale+1)+1):
                for k in range(int(z_min*scale)-1, int(z_max*scale+1)+1):
                    x = i/scale
                    y = j/scale
                    z = k/scale
                    if i%2 == 0 and j%2 == 0 and k%2 == 0 and points[x][y][z] == 1:
                        x1, x2, y1, y2, z1, z2 = x-2/scale, x+2/scale, y-2/scale, y+2/scale, z-2/scale, z+2/scale
                        test_coords = [[x1,y1,z1], [x1,y1,z], [x1,y1,z2], [x1,y,z1], [x1,y,z], [x1,y,z2], [x1,y2,z1], [x1,y2,z], [x1,y2,z2], [x,y1,z1], [x,y1,z], [x,y1,z2], [x,y,z1], [x,y,z2], [x,y2,z1], [x,y2,z], [x,y2,z2], [x2,y1,z1], [x2,y1,z], [x2,y1,z2], [x2,y,z1], [x2,y,z], [x2,y,z2], [x2,y2,z1], [x2,y2,z], [x2,y2,z2]]
                        for coord in test_coords:
                            if x_min <= coord[0] <= x_max and y_min <= coord[1] <= y_max and z_min <= coord[2] <= z_max:
                                if grid[coord[0]][coord[1]][coord[2]] == 1:
                                    grid[(x+coord[0])/2][(y+coord[1])/2][(z+coord[2])/2] = 1
                    # fill in rest with 0s
                    if not (x in grid and y in grid[x] and z in grid[x][y]) or grid[x][y][z] != 1:
                        grid[x][y][z] = 0
        points = grid
    return [grid, scale]

def bounds_xyz(points):
    x_min, x_max, y_min, y_max, z_min, z_max = 1e6, -1e6, 1e6, -1e6, 1e6, -1e6
    for x, v1 in points.items():
        for y, v2 in v1.items():
            for z, v3 in v2.items():
                if v3 == 1:
                    x_min, x_max = min(x_min, x), max(x_max, x)
                    y_min, y_max = min(y_min, y), max(y_max, y)
                    z_min, z_max = min(z_min, z), max(z_max, z)
    return [x_min, x_max, y_min, y_max, z_min, z_max]

def within_region(points, point, scale):
    try:
        return points[float(round(point.x*scale)/scale)][float(round(point.y*scale)/scale)][float(round(point.z*scale)/scale)] == 1
    except:
        return False

def to_radians(angle_degrees):
    return angle_degrees * numpy.pi / 180

def rotate_vector(axis, angle, vector): # clockwise
    axis = axis / norm([float(axis[0]), float(axis[1]), float(axis[2])])
    return Rotation.from_rotvec(to_radians(-angle) * axis).apply(vector)

def bisecting_angle(branch, points, plane):
    x_sum, y_sum, z_sum = 0, 0, 0
    num_points = 0
    for x, v1 in points.items():
        for y, v2 in v1.items():
            for z, p in v2.items():
                if p == 1:
                    num_points += 1
                    x_sum += x
                    y_sum += y
                    z_sum += z
    center_of_mass = Point3D(x_sum/num_points, y_sum/num_points, z_sum/num_points)
    bisecting_line = Line(center_of_mass, branch.end_point)
    bisecting_angle = abs(float(plane.angle_between(bisecting_line))*180/numpy.pi)
    return bisecting_angle

def fix_length(branch, points, length, diameter, direction_vector, scale):
    distance_to_bound = distance_to_edge(branch.end_point, points, direction_vector, scale)
    length_to_diameter = initial_length_to_diameter
    if distance_to_bound == 0:
        distance_to_bound = 1 / scale
    while distance_to_bound / length < 3.0:
        length_to_diameter -= 0.25
        if length_to_diameter == 0:
            break
        length = diameter * length_to_diameter
    while distance_to_bound / length > 6.0:
        length_to_diameter += 0.25
        length = diameter * length_to_diameter
    return length

def distance_to_edge(point, points, direction_vector, scale):
    # binary search to find max distance from branch.end_point to edge of region
    low = 0
    high = max_grid_dimension * scale
    mid = (low + high + 1) // 2
    while low != high:
        end_pt = point + ((mid/scale)*direction_vector[0], (mid/scale)*direction_vector[1], (mid/scale)*direction_vector[2])
        if within_region(points, end_pt, scale):
            low = mid
        else:
            high = mid - 1
        mid = (low + high + 1) // 2
    return low / scale

def auto_dict(): return defaultdict(auto_dict)
points = auto_dict()
total_volume = 0
for x in range(-72, 73):
    if not (-8 <= x <= 8):
        for y in range(-44, 45):
            # for z in range(121):
            for z in range(141):
                if z/4 >= 1*15**(-3)*(((x/4)**2+(1.7*(y/4))**2)**2):
                    points[x/4][y/4][z/4] = 1
                    total_volume += 1/4
                else:
                    points[x/4][y/4][z/4] = 0

# input the first few branches stemming from trachea
trachea_diameter = 1.8
trunk = Branch(Point3D(0, 0, 0), Point3D(0, 0, 7.5), 1, trachea_diameter, trachea_diameter*initial_length_to_diameter, (0, 0, 1), Plane(Point3D(0, 0, 0), normal_vector=(0, 1, 0)), points, 4, total_volume, 0)

count = 0
generate_children(trunk)

#######################
segments_x, segments_y, segments_all, diameters = [ ], [ ], [ ], [ ]
get_segments_x(trunk)
get_segments_y(trunk)
get_segments_all(trunk)
segments = [segments_x, segments_y, segments_all]
with open(savetofile, 'wb') as file:
    pickle.dump(segments, file)


# data = [branch_count, fail_flow_rate, float(airway_volume), total_volume, flow_rates, energy_loss, terminal_volumes, terminal_locations]
# print(data)

# def compute_metrics(data):
#     cubes = [[[0 for z in range(10)] for y in range(7)] for x in range(10)]

#     for t in terminal_locations:
#         cubes[math.floor((t[0]+15)/3)][math.floor((t[1]+10.5)/3)][math.floor(t[2]/3)] += 1
#     nt = [ ]
#     for x in range(len(cubes)):
#         for y in range(len(cubes[0])):
#             for z in range(len(cubes[0][0])):
#                 if cubes[x][y][z] > 0:
#                     nt.append(cubes[x][y][z]/27)
#     print(len(nt))

#     metrics = [data[0], round(data[1], 3), round(numpy.mean(nt)/numpy.std(nt), 2), round(numpy.mean(terminal_volumes)/numpy.std(terminal_volumes), 2), round(data[2]/data[3], 3), round(sum(energy_loss)/len(terminal_volumes), 3)]
#     p_index = sum(metrics) - data[0]
#     metrics.append(round(p_index, 3))

#     return metrics

# metrics = compute_metrics(data)
# print(metrics)

# with open('new_results/coeff_1_05_M.pkl', 'wb') as file:
#     pickle.dump(metrics, file)


########################################################

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

def plot_branches(segments):
    plt.figure(0)
    x = np.linspace(-18, 18) 
    y = 1*15**(-3)*x**4
    plt.plot(x, y)
    for s in segments[0]:
        x, y = s[0], s[1]
        plt.plot(x, y, linewidth=s[2]*8)
    plt.axis('equal')
    plt.gca().invert_yaxis()
    plt.show()

    plt.figure(1)
    x = np.linspace(-10, 10) 
    y = 1*15**(-3)*(1.7*x)**4
    plt.plot(x, y)
    for s in segments[1]:
        x, y = s[0], s[1]
        plt.plot(x, y, linewidth=s[2]*8)
    plt.axis('equal')
    plt.gca().invert_yaxis()
    plt.show()

plot_branches(segments)
