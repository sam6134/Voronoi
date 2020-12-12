# import modules
from bokeh.plotting import figure, output_notebook, show
from bokeh.models import Ray, markers
from math import pi as pi
from numpy import arctan


def add_point(p, x, y):
    p.circle([x], [y], size=7, fill_color="red")


def add_line_segment(p, pt1, pt2):
    (x1, y1), (x2, y2) = pt1, pt2
    p.line([x1, x2], [y1, y2], line_width=2)


def add_ray(p, pt1, unit_vector):
    (x0, y0) = pt1
    (x1, y1) = unit_vector
    if x1 == 0:
        rad_angle = pi / 2
    else:
        rad_angle = arctan(y1 / x1)

    if( x1 <= 0 and y1 > 0):
        rad_angle += pi
    elif( x1 <= 0 and y1 < 0):
        rad_angle += pi

    p.ray(x=x0, y=y0, length=0, angle=rad_angle)


# plot line of form ax+by=c
def add_line(p, a, b, c):

    if(a == 0):
        (x0, y0) = (0,-c/b)
    else:
        (x0, y0) = (-c / a, 0)
    if b == 0:
        rad_angle = pi / 2
    else:
        rad_angle = arctan(-a / b)
    p.ray(x0, y0, 0, rad_angle)
    p.ray(x0, y0, 0, pi + rad_angle)


# if using in jupyter notebook, uncomment the below line to see interactive plot
# output_notebook()

p = figure(plot_width=600, plot_height=600)


with open("shapes.txt") as f:
    shapes = [line.strip() for line in f]
f.close()

for s in shapes:
    # Line
    if s[0] == "L":
        temp = s.split()
        a, b, c = temp[1:4]
        add_line(p, float(a), float(b), float(c))
        print("line", a, b, c)
    # Line segment
    if s[0] == "S":
        temp = s.split()
        x1, y1, x2, y2 = temp[1:5]
        add_line_segment(p, (float(x1), float(y1)), (float(x2), float(y2)))
        print("segment:", (float(x1), float(y1)), " -> ", (float(x2), float(y2)))
    # Ray
    if s[0] == "R":
        temp = s.split()
        x1, y1, x2, y2 = temp[1:5]
        add_ray(p, (float(x1), float(y1)), (float(x2), float(y2)))
        print((float(x1), float(y1)), (float(x2), float(y2)))

    if s[0] == "P":
        temp = s.split()
        x, y = temp[1:3]
        add_point(p, float(x), float(y))

show(p)
