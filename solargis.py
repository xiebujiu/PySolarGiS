__author__ = 'pxiex'

import math
from osgeo import ogr
import numpy as np
from scipy import spatial



__DAY__ = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])


def get_pt_near(wkb_point, map_file, nb=1):
    assert wkb_point.GetGeometryType() == ogr.wkbPoint
    assert map_file[-3:] == 'map'
    point_array = get_pt2d_map(map_file)
    ckd_tree = spatial.cKDTree(point_array)
    return ckd_tree.query(np.asarray(wkb_point.GetPoint_2D(), dtype=np.double), nb)


def unique_point(pt_array):
    array_tmp = pt_array[:, 0] + pt_array[:, 1] * 1j
    array_unique, index = np.unique(array_tmp, return_index=True)
    return pt_array[index]


def get_pt2d_map(map_file):
    assert map_file[-3:] == 'map'
    dri = ogr.GetDriverByName('WAsP')
    map_read = dri.Open(map_file, 0)
    layer = map_read.GetLayer()
    pt_2d_list = []
    for feature in layer:
        geom = feature.GetGeometryRef()
        for i in range(0, geom.GetPointCount()):
            pt = geom.GetPoint_2D(i)
            pt_2d = ogr.Geometry(ogr.wkbPoint)
            pt_2d.AddPoint(pt[0], pt[1])
            pt_2d_list.append(pt_2d.GetPoint())
    pt_2d_array = np.asarray(pt_2d_list, dtype=np.double).reshape(len(np.ravel(pt_2d_list)) / 2, 2)
    return unique_point(pt_2d_array)


def get_pt3d_map(map_file):
    assert map_file[-3:] == 'map'
    dri = ogr.GetDriverByName('WAsP')
    map_read = dri.Open(map_file, 0)
    layer = map_read.GetLayer()
    pt_3d_list = []
    for feature in layer:
        geom = feature.GetGeometryRef()
        elev = feature.GetField("ELEVATION")
        for i in range(0, geom.GetPointCount()):
            pt = geom.GetPoint_2D(i)
            pt_3d = ogr.Geometry(ogr.wkbPoint)
            pt_3d.AddPoint(pt[0], pt[1], elev)
            pt_3d_list.append(pt_3d.GetPoint())
    pt_3d_array = np.asarray(pt_3d_list).reshape(len(np.ravel(pt_3d_list)) / 3, 3)
    pt_3d_list_unique = unique_point(pt_3d_array)
    return np.asarray(pt_3d_list_unique, dtype=np.double)


def get_envelope_map(map_file):
    assert map_file[-3:] == 'map'
    dri = ogr.GetDriverByName('WAsP')
    map_read = dri.Open(map_file, 0)
    layer = map_read.GetLayer()
    geom_collection = ogr.Geometry(ogr.wkbGeometryCollection)
    for feature in layer:
        geom = feature.GetGeometryRef()
        geom_collection.AddGeometry(geom)
    return np.asarray(ogr.Geometry.GetEnvelope(geom_collection), dtype=np.double)


def get_centroid_map(map_file):
    (x_min, x_max, y_min, y_max) = get_envelope_map(map_file)
    print [(x_min + x_max) / 2, (y_min + y_max) / 2]
    exit()
    return np.asarray([(x_min + x_max) / 2, (y_min + y_max) / 2], dtype=np.double)


def get_grid_nb(map_file, res=50):
    assert map_file[-3:] == 'map'
    (x_min, x_max, y_min, y_max) = get_envelope_map(map_file)
    nb_grid_x = int(2 * math.ceil(((x_max - x_min - res) / 2) / res) + 1)
    nb_grid_y = int(2 * math.ceil(((y_max - y_min - res) / 2) / res) + 1)
    return nb_grid_x, nb_grid_y


def get_grid_extreme(map_file, res=50):
    (x_c, y_c) = get_centroid_map(map_file)
    (nb_grid_x, nb_grid_y) = get_grid_nb(map_file, res)
    (grid_x_min, grid_x_max, grid_y_min, grid_y_max) = \
        (x_c - (nb_grid_x - 1) / 2 * res, x_c + (nb_grid_x - 1) / 2 * res,
         y_c - (nb_grid_y - 1) / 2 * res, y_c + (nb_grid_y - 1) / 2 * res)
    return grid_x_min, grid_x_max, grid_y_min, grid_y_max


def gene_ele_grid(map_file, res=50.0, nb_near=5):
    assert map_file[-3:] == 'map'
    (nb_grid_x, nb_grid_y) = get_grid_nb(map_file, res)
    (grid_x_min, grid_x_max, grid_y_min, grid_y_max) = get_grid_extreme(map_file, res)
    pt_3d_array = get_pt3d_map(map_file)
    (x_space, y_space) = (
        np.linspace(grid_x_min, grid_x_max, nb_grid_x), np.linspace(grid_y_min, grid_y_max, nb_grid_y))
    xv, yv = np.meshgrid(x_space, y_space)
    ckd_tree = spatial.cKDTree(pt_3d_array[:, [0, 1]])

    # f_xyz = open(map_file[-3:]+'xyz', 'w')

    def get_ptn_kdt_fix(wkb_point, nb=nb_near):
        assert wkb_point.GetGeometryType() == ogr.wkbPoint
        return ckd_tree.query(np.asarray(wkb_point.GetPoint_2D(), dtype=np.double), nb)

    def get_elev(x, y):
        grid_point = ogr.Geometry(ogr.wkbPoint)
        grid_point.AddPoint_2D(x, y)
        dis_array, index_array = get_ptn_kdt_fix(grid_point, nb_near)
        elev_array = pt_3d_array[index_array, 2]
        if np.any(dis_array == 0):
            return np.mean(pt_3d_array[index_array[dis_array == 0], 2])
        else:
            return np.sum(elev_array * (1 / (dis_array * dis_array)) / np.sum(1 / (dis_array * dis_array)))

    get_elev_grid = np.vectorize(get_elev)
    elev_gpt_array = get_elev_grid(xv, yv)
    array_write = np.vstack((xv.ravel(), yv.ravel(), elev_gpt_array.ravel())).T
    np.savetxt(map_file[:-3] + 'xyz', array_write)
    return np.meshgrid(x_space, y_space), elev_gpt_array


def gene_slope_grid(map_file, res=50.0):
    (x_grid, y_grid), elev_grid = gene_ele_grid(map_file, res)
    size_v = np.size(elev_grid, 0)
    size_h = np.size(elev_grid, 1)

    def get_slope_pt(v, h):
        assert 0 <= v < size_v
        assert 0 <= h < size_h
        #              h-1  h  h+1
        #             === === ===
        #        v+1:| a | b | c |
        #             === === ===
        #        v:  | d | e | f |
        #             === === ===
        #        v-1:| g | h | i |
        #             === === ===
        #        dz/dx = ((c + 2f + i) - (a + 2d + g)) / (8 * dx)
        #        dz/dy = ((g + 2h + i) - (a + 2b + c)) / (8 * dy)
        (v_up, v_down, h_left, h_right) = \
            (v + 1 if v != size_v - 1 else v, v - 1 if v != 0 else v,
             h - 1 if h != 0 else h, h + 1 if h != size_h - 1 else h)
        (a, b, c, d, e, f, g, h, i) = \
            (elev_grid[v_up][h_left], elev_grid[v_up][h], elev_grid[v_up][h_right], \
             elev_grid[v][h_left], elev_grid[v][h], elev_grid[v][h_right],
             elev_grid[v_down][h_left], elev_grid[v_down][h], elev_grid[v_down][h_right])
        dz_dx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * res)
        dz_dy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8 * res)
        rise_run = np.sqrt(dz_dx * dz_dx + dz_dy * dz_dy)
        return np.rad2deg(np.arctan(rise_run))

    get_slope_map = np.vectorize(get_slope_pt)
    h_grid, v_grid = np.meshgrid(np.arange(size_h), np.arange(size_v))
    slope_array = get_slope_map(v_grid, h_grid)
    array_write = np.vstack((x_grid.ravel(), y_grid.ravel(), slope_array.ravel())).T
    np.savetxt(map_file[:-4] + '_slope.txt', array_write)
    return x_grid, y_grid, slope_array


def gene_aspect_grid(map_file, res=50.0):
    (x_grid, y_grid), elev_grid = gene_ele_grid(map_file, res)
    print x_grid.shape
    size_v = np.size(elev_grid, 0)
    size_h = np.size(elev_grid, 1)

    def get_aspect_pt(v, h):
        assert 0 <= v < size_v
        assert 0 <= h < size_h
        #              h-1  h  h+1
        #             === === ===
        #        v+1:| a | b | c |
        #             === === ===
        #        v:  | d | e | f |
        #             === === ===
        #        v-1:| g | h | i |
        #             === === ===
        #        dz/dx = ((c + 2f + i) - (a + 2d + g)) / (8 * dx)
        #        dz/dy = ((g + 2h + i) - (a + 2b + c)) / (8 * dy)
        (v_up, v_down, h_left, h_right) = \
            (v + 1 if v != size_v - 1 else v, v - 1 if v != 0 else v,
             h - 1 if h != 0 else h, h + 1 if h != size_h - 1 else h)
        (a, b, c, d, e, f, g, h, i) = \
            (elev_grid[v_up][h_left], elev_grid[v_up][h], elev_grid[v_up][h_right], \
             elev_grid[v][h_left], elev_grid[v][h], elev_grid[v][h_right],
             elev_grid[v_down][h_left], elev_grid[v_down][h], elev_grid[v_down][h_right])
        dz_dx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * res)
        dz_dy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8 * res)
        print a, b, c, d, e, f, g, h, i
        print dz_dy, dz_dx
        aspect_deg = np.rad2deg(np.arctan2(dz_dy, -dz_dx))

        if aspect_deg <= 90.0:
            aspect_deg = 90.0 - aspect_deg
        else:
            aspect_deg = 360.0 - aspect_deg + 90.0
        return aspect_deg

        '''
        if 0 <= aspect_deg < 22.5 or 337.5 <= aspect_deg <= 360:
            return 0  # north
        elif 22.5 <= aspect_deg < 67.5:
            return 1  # north east
        elif 67.5 <= aspect_deg < 112.5:
            return 2  # east
        elif 112.5 <= aspect_deg < 157.5:
            return 3  # southeast
        elif 157.5 <= aspect_deg < 202.5:
            return 4  # south
        elif 202.5 <= aspect_deg < 247.5:
            return 5  # south west
        elif 247.5 <= aspect_deg < 292.5:
            return 6  # west
        elif 292.5 <= aspect_deg < 337.5:
            return 7  # north
        '''

    get_aspect_map = np.vectorize(get_aspect_pt)
    h_grid, v_grid = np.meshgrid(np.arange(size_h), np.arange(size_v))
    slope_array = get_aspect_map(v_grid, h_grid)
    array_write = np.vstack((x_grid.ravel(), y_grid.ravel(), slope_array.ravel())).T
    np.savetxt(map_file[:-4] + '_aspect.txt', array_write)
    return x_grid, y_grid, slope_array


def xyz2array(xyz_file):
    assert xyz_file[-3:] == 'xyz'
    array_pt = np.loadtxt(xyz_file, np.double)
    x, y, z = np.hsplit(array_pt, 3)
    return x.ravel(), y.ravel(), z.ravel()


def txt2array(txt_file):
    assert txt_file[-3:] == 'txt'
    array_pt = np.loadtxt(txt_file, np.double)
    x, y, z = np.hsplit(array_pt, 3)
    return x.ravel(), y.ravel(), z.ravel()


def get_view_shed_pt(x, y, xyz_file, res=50):
    assert xyz_file[-3:] == 'xyz'
    (array_x, array_y, array_z) = xyz2array(xyz_file)
    ckd_tree = spatial.cKDTree(zip(array_x, array_y))
    index_array = ckd_tree.query(np.asarray([x, y], dtype=np.double))[1]
    (pt_ref_x, pt_ref_y, pt_ref_z) = \
        (array_x[index_array], array_y[index_array], array_z[index_array])
    azi_array = np.array(np.linspace(0, 359, 360))
    array_c = array_x + array_y * 1j

    array_x = np.delete(array_x, index_array)
    array_y = np.delete(array_y, index_array)
    array_z = np.delete(array_z, index_array)
    array_c = np.delete(array_c, index_array)

    zone1 = np.argwhere(np.logical_and(array_x >= pt_ref_x, array_y >= pt_ref_y)).reshape(-1)
    zone2 = np.argwhere(np.logical_and(array_x <= pt_ref_x, array_y >= pt_ref_y)).reshape(-1)
    zone3 = np.argwhere(np.logical_and(array_x <= pt_ref_x, array_y <= pt_ref_y)).reshape(-1)
    zone4 = np.argwhere(np.logical_and(array_x >= pt_ref_x, array_y <= pt_ref_y)).reshape(-1)

    def create_line(x_tmp, y_tmp):
        l = ogr.Geometry(ogr.wkbLineString)
        if x_tmp >= pt_ref_x and y_tmp >= pt_ref_y:
            l.AddPoint_2D(x_tmp - res / 2, y_tmp + res / 2)
            l.AddPoint_2D(x_tmp + res / 2, y_tmp - res / 2)
        elif x_tmp <= pt_ref_x and y_tmp >= pt_ref_y:
            l.AddPoint_2D(x_tmp - res / 2, y_tmp - res / 2)
            l.AddPoint_2D(x_tmp + res / 2, y_tmp + res / 2)
        elif x_tmp <= pt_ref_x and y_tmp <= pt_ref_y:
            l.AddPoint_2D(x_tmp - res / 2, y_tmp + res / 2)
            l.AddPoint_2D(x_tmp + res / 2, y_tmp - res / 2)
        elif x_tmp >= pt_ref_x and y_tmp <= pt_ref_y:
            l.AddPoint_2D(x_tmp - res / 2, y_tmp - res / 2)
            l.AddPoint_2D(x_tmp + res / 2, y_tmp + res / 2)
        return l

    create_line_z = np.vectorize(create_line)
    line_array_z = create_line_z(array_x, array_y)
    r = np.max(np.abs(array_c - (pt_ref_x + pt_ref_y * 1j)))

    def get_azi_vs(azi_tmp):
        line_ref = ogr.Geometry(ogr.wkbLineString)
        line_ref.AddPoint_2D(pt_ref_x, pt_ref_y)
        line_ref.AddPoint_2D(pt_ref_x + r * math.cos(np.deg2rad(azi_tmp)),
                             pt_ref_y + r * math.sin(np.deg2rad(azi_tmp)))

        def get_atan_laz(zone):
            ip_array = np.array([])
            for line in line_array_z[zone]:
                if ogr.Geometry.Intersects(line, line_ref):
                    ip_array = np.append(ip_array, np.argwhere(line_array_z == line)).astype(int)
            ptc_array = array_x[ip_array] + array_y[ip_array] * 1j
            ptz_array = array_z[ip_array] - pt_ref_z
            ptz_array[np.argwhere(ptz_array < 0)] = 0.0
            return np.rad2deg(np.max(np.arctan(ptz_array / np.abs(ptc_array - (pt_ref_x + pt_ref_y * 1j)))))

        if 0.0 <= azi_tmp < 90.0:
            return get_atan_laz(zone1)
        elif 90.0 <= azi_tmp < 180.0:
            return get_atan_laz(zone2)
        elif 180.0 <= azi_tmp < 270.0:
            return get_atan_laz(zone3)
        elif 270.0 <= azi_tmp < 360.0:
            return get_atan_laz(zone4)
        line_ref.Destroy()

    get_vs = np.vectorize(get_azi_vs)
    view_shed_array = get_vs(azi_array)
    return azi_array, view_shed_array


def get_declination_angle(d_index):
    return np.deg2rad(23.45 * np.sin(2 * np.pi * (284 + d_index) / 365))


def get_hour_angle(h_index):
    return np.deg2rad(15 * (h_index - 12))


def get_sunset_angle(d_index, lat):
    def get_sunset_ele(d):
        if (-np.tan(get_declination_angle(d)) * np.tan(np.deg2rad(lat))) > 1.0:
            return np.arccos(1.0)
        elif (-np.tan(get_declination_angle(d)) * np.tan(np.deg2rad(lat))) < -1.0:
            return np.arccos(-1.0)
        else:
            return np.arccos(-np.tan(get_declination_angle(d)) * np.tan(np.deg2rad(lat)))

    f = np.vectorize(get_sunset_ele)
    return f(d_index)


def get_sunset_time(d_index, lat):
    def get_sunset_time_ele(d):
        sunset = np.rad2deg(get_sunset_angle(d, lat))
        return sunset / 15 + 12

    f = np.vectorize(get_sunset_time_ele)
    return f(d_index)


def get_zenith_angle(d_index, h_index, lat):
    def get_zenith_ele(d, h):
        term1 = np.sin(np.deg2rad(lat)) * np.sin(get_declination_angle(d))
        term2 = np.cos(np.deg2rad(lat)) * np.cos(get_declination_angle(d)) * np.cos(get_hour_angle(h))

        if np.arccos(term1 + term2) <= np.pi / 2:
            return np.arccos(term1 + term2)

        else:
            return np.pi / 2

            # return np.arccos(term1 + term2)

    f = np.vectorize(get_zenith_ele)
    return f(d_index, h_index)


def get_azi_angle(d_index, h_index, lat):
    def get_azimuth_angle_ele(d, h):
        term1 = np.sin(np.deg2rad(lat)) * np.cos(get_zenith_angle(d, h, lat))
        term2 = np.sin(get_declination_angle(d))
        term3 = np.cos(np.deg2rad(lat)) * np.sin(get_zenith_angle(d, h, lat))
        if (term1 - term2) / term3 > 1:
            azi = np.arccos(1.0)
        elif (term1 - term2) / term3 < -1:
            azi = np.arccos(-1.0)
        else:
            azi = np.arccos((term1 - term2) / term3)
        if (get_hour_angle(h)) > 0:
            azi = np.mod((azi + np.pi), (2 * np.pi))
        else:
            azi = np.mod((3 * np.pi - azi), (2 * np.pi))
        if azi <= np.pi:
            azi = azi
        else:
            azi -= 2 * np.pi
        if azi >= 0:
            return azi - np.pi
        else:
            return azi + np.pi
            # return azi

    f = np.vectorize(get_azimuth_angle_ele)

    return f(d_index, h_index)


def get_azi_math(azi):
    azi = -azi
    azi -= np.pi / 2
    return azi


def get_azi_sunset(d_index, lat):
    def get_azi_sunset_ele(d):
        term1 = 0.0
        term2 = np.sin(get_declination_angle(d))
        term3 = np.cos(np.deg2rad(lat))
        if (term1 - term2) / term3 > 1:
            azi = np.arccos(1.0)
        elif (term1 - term2) / term3 < -1:
            azi = np.arccos(-1.0)
        else:
            azi = np.arccos((term1 - term2) / term3)
        azi = np.mod((3 * np.pi - azi), (2 * np.pi))
        if azi <= np.pi:
            azi = azi
        else:
            azi -= 2 * np.pi
        if azi >= 0:
            return azi - np.pi
        else:
            return azi + np.pi
            # return azi

    f = np.vectorize(get_azi_sunset_ele)

    return f(d_index)


def get_sun_map(lat):
    day_interval = np.arange(173, 356, 14)
    hour_interval = np.arange(6, 18, 0.5)
    hour_v, day_v = np.meshgrid(hour_interval, day_interval)

    zenith_v = get_zenith_angle(day_v, hour_v, lat)
    zenith_add = np.array([[np.pi / 2] * len(day_interval)]).reshape(len(day_interval), 1)
    zenith_v = np.hstack((zenith_add, zenith_v, zenith_add))
    azimuth_v = get_azi_angle(day_v, hour_v, lat)
    azimuth_add = get_azi_sunset(day_interval, lat).reshape(len(day_interval), 1)
    azimuth_v = np.hstack((azimuth_add, azimuth_v, -azimuth_add))

    hour_add = get_sunset_time(day_interval, lat).reshape(len(day_interval), 1)
    hour_v = np.hstack((hour_add - 2 * (hour_add - 12), hour_v, hour_add))
    for hour_line in hour_v:
        hour_line[hour_line < hour_line[0]] = hour_line[0]
        hour_line[hour_line > hour_line[-1]] = hour_line[-1]

    day_add = day_interval.reshape(len(day_interval), 1)
    day_v = np.hstack((day_add, day_v, day_add))

    return np.rad2deg(zenith_v), get_azi_math(azimuth_v), day_v, hour_v


def get_sun_sec(lat):
    (zenith_v, azimuth_v, day_v, hour_v) = get_sun_map(lat)
    zenith_c, azimuth_c = np.meshgrid(np.zeros(len(zenith_v[0]) - 1),
                                      np.zeros(len(zenith_v) - 1))
    hour_p, day_p = np.meshgrid(np.zeros(len(zenith_v[0]) - 1),
                                np.zeros(len(zenith_v) - 1))

    for i in range(len(zenith_c)):
        for j in range(len(zenith_c[0])):
            zenith_c[i][j] = \
                (zenith_v[i][j] + zenith_v[i + 1][j] + zenith_v[i][j + 1] + zenith_v[i + 1][j + 1]) / 4
            azimuth_c[i][j] = \
                (azimuth_v[i][j] + azimuth_v[i + 1][j] + azimuth_v[i][j + 1] + azimuth_v[i + 1][j + 1]) / 4
            hour_p[i][j] = \
                ((hour_v[i][j + 1] - hour_v[i][j]) + (hour_v[i + 1][j + 1] - hour_v[i + 1][j])) / 2
            day_p[i][j] = \
                ((day_v[i + 1][j] - day_v[i][j]) + (day_v[i + 1][j + 1] - day_v[i][j + 1])) / 2

    time_p = hour_p * day_p
    return zenith_c, azimuth_c, time_p


def get_incidence_angle(zeni, azi, array_slope, array_aspect):
    incidence = np.arccos(np.cos(np.deg2rad(zeni) * np.cos(np.deg2rad(array_slope))
                                 + np.sin(np.deg2rad(zeni)) * np.sin(np.deg2rad(array_slope))
                                 * np.cos(azi - np.deg2rad(array_aspect))))

    return incidence


def get_sky_map(sec_zenith=8, sec_azimuth=16):
    zeni_sky = np.linspace(0.0, 90.0, sec_zenith + 1)
    azi_sky = np.deg2rad(np.linspace(0.0, 360.0, sec_azimuth + 1))
    azi_v, zeni_v = np.meshgrid(azi_sky, zeni_sky)
    return zeni_v, azi_v


def get_sky_sec(sec_zenith=8, sec_azimuth=16):
    zeni_c_s = 90.0 / sec_zenith / 2
    azi_c_s = 360.0 / sec_azimuth / 2
    zeni_c = np.arange(zeni_c_s, 90.0, 90.0 / sec_zenith)
    azi_c = np.arange(azi_c_s, 360.0, 360.0 / sec_azimuth)
    azi_c_v, zeni_c_v = np.meshgrid(azi_c, zeni_c)
    return zeni_c_v, azi_c_v


def get_dir_map(xyz_file, slope_file, aspect_file, lat, beta=0.70, res=50):
    assert xyz_file[-3:] == 'xyz'
    x_array, y_array, z_array = xyz2array(xyz_file)
    x_array, y_array, slope_array = txt2array(slope_file)
    x_array, y_array, aspect_array = txt2array(aspect_file)
    zeni_c, azi_c, time_p = get_sun_sec(lat)
    ckd_tree = spatial.cKDTree(zip(x_array, y_array))
    zeni_v = get_sun_map(lat)[0]
    azi_v = get_sun_map(lat)[1]

    def get_sun_gap_ele(x_tmp, y_tmp):

        def get_view_shed_ele(x_tmp2, y_tmp2):
            index_array = ckd_tree.query(np.asarray([x_tmp2, y_tmp2], dtype=np.double))[1]
            (pt_ref_x, pt_ref_y, pt_ref_z) = \
                (x_array[index_array], y_array[index_array], z_array[index_array])
            azi_array = np.array(np.linspace(0, 359, 360))
            array_c = x_array + y_array * 1j

            array_x = np.delete(x_array, index_array)
            array_y = np.delete(y_array, index_array)
            array_z = np.delete(z_array, index_array)
            array_c = np.delete(array_c, index_array)

            zone1 = np.argwhere(np.logical_and(array_x >= pt_ref_x, array_y >= pt_ref_y)).reshape(-1)
            zone2 = np.argwhere(np.logical_and(array_x <= pt_ref_x, array_y >= pt_ref_y)).reshape(-1)
            zone3 = np.argwhere(np.logical_and(array_x <= pt_ref_x, array_y <= pt_ref_y)).reshape(-1)
            zone4 = np.argwhere(np.logical_and(array_x >= pt_ref_x, array_y <= pt_ref_y)).reshape(-1)

            def create_line(x_temp, y_temp):
                l = ogr.Geometry(ogr.wkbLineString)
                if x_temp >= pt_ref_x and y_temp >= pt_ref_y:
                    l.AddPoint_2D(x_temp - res / 2, y_temp + res / 2)
                    l.AddPoint_2D(x_temp + res / 2, y_temp - res / 2)
                elif x_temp <= pt_ref_x and y_temp >= pt_ref_y:
                    l.AddPoint_2D(x_temp - res / 2, y_temp - res / 2)
                    l.AddPoint_2D(x_temp + res / 2, y_temp + res / 2)
                elif x_temp <= pt_ref_x and y_temp <= pt_ref_y:
                    l.AddPoint_2D(x_temp - res / 2, y_temp + res / 2)
                    l.AddPoint_2D(x_temp + res / 2, y_temp - res / 2)
                elif x_temp >= pt_ref_x and y_temp <= pt_ref_y:
                    l.AddPoint_2D(x_temp - res / 2, y_temp - res / 2)
                    l.AddPoint_2D(x_temp + res / 2, y_temp + res / 2)
                return l

            create_line_z = np.vectorize(create_line)
            line_array_z = create_line_z(array_x, array_y)
            r = np.max(np.abs(array_c - (pt_ref_x + pt_ref_y * 1j)))

            def get_azi_vs(azi_tmp):
                line_ref = ogr.Geometry(ogr.wkbLineString)
                line_ref.AddPoint_2D(pt_ref_x, pt_ref_y)
                line_ref.AddPoint_2D(pt_ref_x + 2 * r * math.cos(np.deg2rad(azi_tmp)),
                                     pt_ref_y + 2 * r * math.sin(np.deg2rad(azi_tmp)))

                def get_atan_laz(zone):
                    ip_array = np.array([])
                    for line in line_array_z[zone]:
                        if ogr.Geometry.Intersects(line, line_ref):
                            ip_array = np.append(ip_array, np.argwhere(line_array_z == line)).astype(int)
                    if len(ip_array) != 0:
                        ptc_array = array_x[ip_array] + array_y[ip_array] * 1j
                        ptz_array = array_z[ip_array] - pt_ref_z
                        ptz_array[np.argwhere(ptz_array < 0)] = 0.0
                        return np.rad2deg(np.max(np.arctan(ptz_array / np.abs(ptc_array - (pt_ref_x + pt_ref_y * 1j)))))
                    else:
                        # print ip_array, len(ip_array)
                        return 0.0

                if 0.0 <= azi_tmp < 90.0:
                    return get_atan_laz(zone1)
                elif 90.0 <= azi_tmp < 180.0:
                    return get_atan_laz(zone2)
                elif 180.0 <= azi_tmp < 270.0:
                    return get_atan_laz(zone3)
                elif 270.0 <= azi_tmp < 360.0:
                    return get_atan_laz(zone4)
                line_ref.Destroy()

            get_vs = np.vectorize(get_azi_vs)
            view_shed_array = get_vs(azi_array)
            # print azi_array, view_shed_array
            return azi_array, view_shed_array

        azi_vs, vs_array = get_view_shed_ele(x_tmp, y_tmp)
        zenith_vs = 90 - vs_array

        poly_ref = ogr.Geometry(ogr.wkbPolygon)
        linear_ring = ogr.Geometry(ogr.wkbLinearRing)
        for k_tmp in range(len(zenith_vs)):
            linear_ring.AddPoint_2D(zenith_vs[k_tmp] * np.cos(np.deg2rad(azi_vs[k_tmp])),
                                    zenith_vs[k_tmp] * np.sin(np.deg2rad(azi_vs[k_tmp])))
        linear_ring.AddPoint_2D(zenith_vs[0] * np.cos(np.deg2rad(azi_vs[0])),
                                zenith_vs[0] * np.sin(np.deg2rad(azi_vs[0])))
        poly_ref.AddGeometry(linear_ring)

        def get_sun_gap_position(i_tmp, j_tmp):
            poly = ogr.Geometry(ogr.wkbPolygon)
            ls = ogr.Geometry(ogr.wkbLinearRing)
            ls.AddPoint_2D(zeni_v[i_tmp][j_tmp], azi_v[i_tmp][j_tmp])
            ls.AddPoint_2D(zeni_v[i_tmp][j_tmp + 1], azi_v[i_tmp][j_tmp + 1])
            ls.AddPoint_2D(zeni_v[i_tmp + 1][j_tmp + 1], azi_v[i_tmp + 1][j_tmp + 1])
            ls.AddPoint_2D(zeni_v[i_tmp + 1][j_tmp], azi_v[i_tmp + 1][j_tmp])
            ls.AddPoint_2D(zeni_v[i_tmp][j_tmp], azi_v[i_tmp][j_tmp])
            poly.AddGeometry(ls)

            intersection = poly_ref.Intersection(poly)
            if poly_ref.Intersects(poly):
                area = intersection.GetArea()
            else:
                area = 0.0
            if poly.GetArea() == 0.0:
                sun_gap_tmp = 0.0
            else:
                sun_gap_tmp = area / poly.GetArea()
            if sun_gap_tmp < 0:
                return 0.0
            else:
                return sun_gap_tmp

        array_i = np.arange(len(zeni_v[0]) - 1)
        array_j = np.arange(len(zeni_v) - 1)
        j_v, i_v = np.meshgrid(array_i, array_j)
        f_sun_gap = np.vectorize(get_sun_gap_position)
        sun_gap_return = f_sun_gap(i_v, j_v)
        # print sun_gap_return
        return sun_gap_return

    def get_dir_tot_ele(x, y, z, slope, aspect):
        sun_gap = get_sun_gap_ele(x, y)
        ele = z

        def get_incidence_ele(zeni, azi, slope_in, aspect_in):
            inc = np.arccos(np.cos(np.deg2rad(zeni)) * np.cos(np.deg2rad(slope_in))
                            + np.sin(np.deg2rad(zeni)) * np.sin(np.deg2rad(slope_in))
                            * np.cos(azi - np.deg2rad(aspect_in)))
            return inc

        f_incidence = np.vectorize(get_incidence_ele, excluded=['slope_in', 'aspect_in'])
        incidence = f_incidence(zeni_c, azi_c, slope, aspect)

        def get_trans_ele(zeni, elevation, transmittance):
            m = np.exp(-0.000118 * elevation - 1.638 * np.power(10, -9) * elevation * elevation) / \
                (np.cos(np.deg2rad(zeni)) + 0.50572 * np.power((96.07995 - zeni), -1.6364))
            return np.power(transmittance, m)

        f_trans = np.vectorize(get_trans_ele, excluded=['elevation', 'transmittance'])
        trans = f_trans(zeni_c, ele, beta)

        def get_dir_ele(time_p_ele, sun_gap_ele, incidence_ele, trans_ele):
            # incidence = get_incidence_angle(zeni_c_ele, azi_c_ele, slope, aspect)
            dir_ele = 1.367 * trans_ele * time_p_ele * sun_gap_ele * np.cos(incidence_ele)
            if dir_ele < 0.0:
                return 0.0
            return dir_ele

        f_dir = np.vectorize(get_dir_ele)
        dir_array = f_dir(time_p, sun_gap, incidence, trans)
        dir_tot = np.sum(dir_array)
        return dir_tot

    f_dir_tot = np.vectorize(get_dir_tot_ele)
    dir_tot_array = f_dir_tot(x_array, y_array, z_array, slope_array, aspect_array)

    np.savetxt(xyz_file[:-3] + '.txt', np.vstack((x_array.ravel(),
                                                  y_array.ravel(),
                                                  dir_tot_array.ravel())).T)
    return x_array, y_array, dir_tot_array


def get_dif_map(xyz_file, slope_file, aspect_file, lat,
                beta=0.70, p_dif=0.70, res=50, sec_zenith=8, sec_azimuth=16):
    assert xyz_file[-3:] == 'xyz'
    x_array, y_array, z_array = xyz2array(xyz_file)
    x_array, y_array, slope_array = txt2array(slope_file)
    x_array, y_array, aspect_array = txt2array(aspect_file)
    zeni_c_sky, azi_c_sky = get_sky_sec(sec_zenith, sec_azimuth)
    zeni_c, azi_c, time_p = get_sun_sec(lat)
    ckd_tree = spatial.cKDTree(zip(x_array, y_array))
    zeni_sky, azi_sky = get_sky_map(sec_zenith, sec_azimuth)

    def get_sky_gap_ele(x_tmp, y_tmp):

        def get_view_shed_ele(x_tmp2, y_tmp2):
            index_array = ckd_tree.query(np.asarray([x_tmp2, y_tmp2], dtype=np.double))[1]
            (pt_ref_x, pt_ref_y, pt_ref_z) = \
                (x_array[index_array], y_array[index_array], z_array[index_array])
            azi_array = np.array(np.linspace(0, 359, 360))
            array_c = x_array + y_array * 1j

            array_x = np.delete(x_array, index_array)
            array_y = np.delete(y_array, index_array)
            array_z = np.delete(z_array, index_array)
            array_c = np.delete(array_c, index_array)

            zone1 = np.argwhere(np.logical_and(array_x >= pt_ref_x, array_y >= pt_ref_y)).reshape(-1)
            zone2 = np.argwhere(np.logical_and(array_x <= pt_ref_x, array_y >= pt_ref_y)).reshape(-1)
            zone3 = np.argwhere(np.logical_and(array_x <= pt_ref_x, array_y <= pt_ref_y)).reshape(-1)
            zone4 = np.argwhere(np.logical_and(array_x >= pt_ref_x, array_y <= pt_ref_y)).reshape(-1)

            def create_line(x_temp, y_temp):
                l = ogr.Geometry(ogr.wkbLineString)
                if x_temp >= pt_ref_x and y_temp >= pt_ref_y:
                    l.AddPoint_2D(x_temp - res / 2, y_temp + res / 2)
                    l.AddPoint_2D(x_temp + res / 2, y_temp - res / 2)
                elif x_temp <= pt_ref_x and y_temp >= pt_ref_y:
                    l.AddPoint_2D(x_temp - res / 2, y_temp - res / 2)
                    l.AddPoint_2D(x_temp + res / 2, y_temp + res / 2)
                elif x_temp <= pt_ref_x and y_temp <= pt_ref_y:
                    l.AddPoint_2D(x_temp - res / 2, y_temp + res / 2)
                    l.AddPoint_2D(x_temp + res / 2, y_temp - res / 2)
                elif x_temp >= pt_ref_x and y_temp <= pt_ref_y:
                    l.AddPoint_2D(x_temp - res / 2, y_temp - res / 2)
                    l.AddPoint_2D(x_temp + res / 2, y_temp + res / 2)
                return l

            create_line_z = np.vectorize(create_line)
            line_array_z = create_line_z(array_x, array_y)
            r = np.max(np.abs(array_c - (pt_ref_x + pt_ref_y * 1j)))

            def get_azi_vs(azi_tmp):
                line_ref = ogr.Geometry(ogr.wkbLineString)
                line_ref.AddPoint_2D(pt_ref_x, pt_ref_y)
                line_ref.AddPoint_2D(pt_ref_x + 2 * r * math.cos(np.deg2rad(azi_tmp)),
                                     pt_ref_y + 2 * r * math.sin(np.deg2rad(azi_tmp)))

                def get_atan_laz(zone):
                    ip_array = np.array([])
                    for line in line_array_z[zone]:
                        if ogr.Geometry.Intersects(line, line_ref):
                            ip_array = np.append(ip_array, np.argwhere(line_array_z == line)).astype(int)
                    if len(ip_array) != 0:
                        ptc_array = array_x[ip_array] + array_y[ip_array] * 1j
                        ptz_array = array_z[ip_array] - pt_ref_z
                        ptz_array[np.argwhere(ptz_array < 0)] = 0.0
                        return np.rad2deg(np.max(np.arctan(ptz_array / np.abs(ptc_array - (pt_ref_x + pt_ref_y * 1j)))))
                    else:
                        # print ip_array, len(ip_array)
                        return 0.0

                if 0.0 <= azi_tmp < 90.0:
                    return get_atan_laz(zone1)
                elif 90.0 <= azi_tmp < 180.0:
                    return get_atan_laz(zone2)
                elif 180.0 <= azi_tmp < 270.0:
                    return get_atan_laz(zone3)
                elif 270.0 <= azi_tmp < 360.0:
                    return get_atan_laz(zone4)
                line_ref.Destroy()

            get_vs = np.vectorize(get_azi_vs)
            view_shed_array = get_vs(azi_array)
            # print azi_array, view_shed_array
            return azi_array, view_shed_array

        azi_vs, vs_array = get_view_shed_ele(x_tmp, y_tmp)
        zenith_vs = 90 - vs_array

        poly_ref = ogr.Geometry(ogr.wkbPolygon)
        linear_ring = ogr.Geometry(ogr.wkbLinearRing)
        for k_tmp in range(len(zenith_vs)):
            linear_ring.AddPoint_2D(zenith_vs[k_tmp] * np.cos(np.deg2rad(azi_vs[k_tmp])),
                                    zenith_vs[k_tmp] * np.sin(np.deg2rad(azi_vs[k_tmp])))
        linear_ring.AddPoint_2D(zenith_vs[0] * np.cos(np.deg2rad(azi_vs[0])),
                                zenith_vs[0] * np.sin(np.deg2rad(azi_vs[0])))
        poly_ref.AddGeometry(linear_ring)

        def get_sky_gap_position(i_tmp, j_tmp):
            poly = ogr.Geometry(ogr.wkbPolygon)
            ls = ogr.Geometry(ogr.wkbLinearRing)
            ls.AddPoint_2D(zeni_sky[i_tmp][j_tmp], azi_sky[i_tmp][j_tmp])
            ls.AddPoint_2D(zeni_sky[i_tmp][j_tmp + 1], azi_sky[i_tmp][j_tmp + 1])
            ls.AddPoint_2D(zeni_sky[i_tmp + 1][j_tmp + 1], azi_sky[i_tmp + 1][j_tmp + 1])
            ls.AddPoint_2D(zeni_sky[i_tmp + 1][j_tmp], azi_sky[i_tmp + 1][j_tmp])
            ls.AddPoint_2D(zeni_sky[i_tmp][j_tmp], azi_sky[i_tmp][j_tmp])
            poly.AddGeometry(ls)

            intersection = poly_ref.Intersection(poly)
            if poly_ref.Intersects(poly):
                area = intersection.GetArea()
            else:
                area = 0.0
            if poly.GetArea() == 0.0:
                sky_gap_tmp = 0.0
            else:
                sky_gap_tmp = area / poly.GetArea()
            if sky_gap_tmp < 0:
                return 0.0
            else:
                return sky_gap_tmp

        array_i = np.arange(len(zeni_c_sky[0]))
        array_j = np.arange(len(zeni_c_sky))
        j_v, i_v = np.meshgrid(array_i, array_j)
        f_sun_gap = np.vectorize(get_sky_gap_position)
        sky_gap_return = f_sun_gap(i_v, j_v)
        # print sun_gap_return
        return sky_gap_return

    def get_dif_tot_ele(x, y, z, slope, aspect):
        sky_gap = get_sky_gap_ele(x, y)
        ele = z

        def get_trans_ele(zeni, elevation, transmittance):
            m = np.exp(-0.000118 * elevation - 1.638 * np.power(10, -9) * elevation * elevation) / \
                (np.cos(np.deg2rad(zeni)) + 0.50572 * np.power((96.07995 - zeni), -1.6364))
            return np.power(transmittance, m)

        f_trans = np.vectorize(get_trans_ele, excluded=['elevation', 'transmittance'])
        trans_s = f_trans(zeni_c, ele, beta)

        def get_glb_ele(trans_s_ele, time_p_ele):
            glb_ele = 1.376 * (1 - trans_s_ele) * time_p_ele
            return glb_ele

        f_glb = np.vectorize(get_glb_ele)
        glb = np.sum(f_glb(trans_s, time_p))

        def get_incidence_ele(zeni, azi, slope_in, aspect_in):
            inc = np.arccos(np.cos(np.deg2rad(zeni)) * np.cos(np.deg2rad(slope_in))
                            + np.sin(np.deg2rad(zeni)) * np.sin(np.deg2rad(slope_in))
                            * np.cos(azi - np.deg2rad(aspect_in)))
            return inc

        f_incidence = np.vectorize(get_incidence_ele, excluded=['slope_in', 'aspect_in'])
        incidence = f_incidence(zeni_c_sky, azi_c_sky, slope, aspect)

        def get_weight_ele(i, j):
            return (np.cos(np.deg2rad(zeni_sky[i][j])) - np.cos(np.deg2rad(zeni_sky[i + 1][j]))) / sec_azimuth

        f_weight = np.vectorize(get_weight_ele)
        i_array = np.arange(len(zeni_c_sky[0]))
        j_array = np.arange(len(zeni_c_sky))
        j_v, i_v = np.meshgrid(i_array, j_array)
        weight = f_weight(i_v, j_v)

        def get_dif_ele(sky_gap_ele, incidence_ele, weight_ele):
            # incidence = get_incidence_angle(zeni_c_ele, azi_c_ele, slope, aspect)
            dif_ele = glb * p_dif * weight_ele * sky_gap_ele * np.cos(incidence_ele)
            if dif_ele < 0.0:
                return 0.0
            else:
                return dif_ele

        f_dif = np.vectorize(get_dif_ele)
        dif_array = f_dif(sky_gap, incidence, weight)
        dif_tot = np.sum(dif_array)
        return dif_tot

    f_dif_tot = np.vectorize(get_dif_tot_ele)
    dif_tot_array = f_dif_tot(x_array, y_array, z_array, slope_array, aspect_array)

    np.savetxt(xyz_file[:-3] + '.txt', np.vstack((x_array.ravel(),
                                                  y_array.ravel(),
                                                  dif_tot_array.ravel())).T)
    return x_array, y_array, dif_tot_array


def get_proximity_map(xyz_file, element_shp):
    assert xyz_file[-3:] == 'xyz'
    assert element_shp[-3:] == 'shp'
    x_array, y_array, z_array = xyz2array(xyz_file)
    driver = ogr.GetDriverByName("ESRI Shapefile")
    element_data = driver.Open(element_shp, 0)
    element_layer = element_data.GetLayer(0)

    def get_distance_ele(x, y, fea):
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint_2D(x, y)
        distance = point.Distance(fea)
        return distance

    f_distance = np.vectorize(get_distance_ele, excluded=['fea'])
    dis_mat_min = np.array([10000000] * len(x_array))

    for feature in element_layer:
        geom = feature.GetGeometryRef()
        dis = f_distance(x_array, y_array, geom)
        dis_min = np.minimum(dis_mat_min, dis)

    np.savetxt(xyz_file[:-3] + '_' + element_shp[:-3] + '_distance.txt', np.vstack((x_array.ravel(),
                                                                                    y_array.ravel(),
                                                                                    dis_min.ravel())).T)
    return x_array, y_array, dis_min


def get_constraint_map(xyz_file, element_shp):
    assert xyz_file[-3:] == 'xyz'
    assert element_shp[-3:] == 'shp'
    x_array, y_array, z_array = xyz2array(xyz_file)
    driver = ogr.GetDriverByName("ESRI Shapefile")
    element_data = driver.Open(element_shp, 0)
    element_layer = element_data.GetLayer(0)

    def get_constraint_ele(x, y, fea):
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint_2D(x, y)
        if point.Intersects(fea):
            return 0
        else:
            return 1

    f_constraint = np.vectorize(get_constraint_ele, excluded=['fea'])
    con_union = np.array([1] * len(x_array))

    for feature in element_layer:
        geom = feature.GetGeometryRef()
        constraint = f_constraint(x_array, y_array, geom)
        con_union = np.minimum(con_union, constraint)

    np.savetxt(xyz_file[:-3] + '_' + element_shp[:-3] + '_constraint.txt', np.vstack((x_array.ravel(),
                                                                                      y_array.ravel(),
                                                                                      con_union.ravel())).T)
    return x_array, y_array, con_union


def get_site_map(dir_file, dif_file, aspect_file, slope_file, dis_file):
    x_dir, y_dir, z_dir = txt2array(dir_file)
    x_dif, y_dif, z_dif = txt2array(dif_file)
    x_as, y_as, z_as = txt2array(aspect_file)
    x_sl, y_sl, z_sl = txt2array(slope_file)
    x_dis, y_dis, z_dis = txt2array(dis_file)
    z_tot = z_dir + z_dif
    z_as = np.abs(z_as - 180)
    w = np.array([0.66898093, 0.05147005, 0.19953907, 0.08000995])

    tot_range = np.max(z_tot) - np.min(z_tot)
    as_range = np.max(z_as) - np.min(z_as)
    sl_range = np.max(z_sl) - np.min(z_sl)
    dis_range = np.max(z_dis) - np.min(z_dis)

    tot_score = (z_tot - np.min(z_tot)) / (tot_range / 10)
    print tot_score[:100]
    as_score = 100 - (z_as - np.min(z_as)) / (as_range / 10)
    sl_score = 100 - (z_sl - np.min(z_sl)) / (sl_range / 10)
    dis_score = (z_dis - np.min(z_dis)) / (dis_range / 10)

    final_score = np.vstack((tot_score.ravel(), dis_score.ravel(), as_score.ravel(), sl_score.ravel())).T
    final_score = np.dot(final_score, w)
    final_range = np.max(final_score) - np.min(final_score)
    final_score = (final_score - np.min(final_score)) / (final_range / 10)
    np.savetxt('site_selection.txt', np.vstack((x_dir.ravel(), y_dir.ravel(), final_score.ravel())).T)
    return x_dir, y_dir, final_score


gene_aspect_grid('./Topology.map', res=10.0)
