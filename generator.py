import math
import random
import os
import json

from models import Ownship, Intruder
from daa_logic_docker import *

'''
    Формирование траекторий .daa со столкновением

    Структура выходных данных:
    "NAME, lat, lon, alt, vx, vy, vz, time"
    "unitless, [deg], [deg], [ft], [knot], [knot], [fpm], [s]"
'''
# Параметры WGS-84
EARTH_A = 6378137.0
EARTH_E = 0.0818191908426

#  Константы безопасности
NM_TO_M = 1852.0  # морской мили в метры
type_conflict = "NMAC"  # или "LOWC"

import math


def vx_vy_to_heading(vx, vy):

    gs = math.sqrt(vx ** 2 + vy ** 2)

    heading_rad = math.atan2(vx, vy)

    heading_deg = math.degrees(heading_rad)
    if heading_deg < 0:
        heading_deg += 360

    return {'track': heading_deg, 'gs' :gs}
if type_conflict == "LOWC":
    SAFE_HORIZ_M = 609.9
    SAFE_VERT_FT = 500.0
elif type_conflict == "NMAC":
    SAFE_HORIZ_M = 150.0
    SAFE_VERT_FT = 100.0
else:
    raise ValueError("type_conflict must be 'NMAC' or 'LOWC'")

NM_TO_DEG = 1.0 / 60.0
KNOT_TO_MS = 0.514444


def heading_to_vx_vy(heading_deg, speed_knot):
    rad = math.radians(heading_deg)
    vx = speed_knot * math.sin(rad)  # Восток (KNOT)
    vy = speed_knot * math.cos(rad)  # Север (KNOT)
    return vx, vy


EARTH_A = 6378137.0
EARTH_E = 0.0818191908426
EARTH_F = 1 / 298.257223563  # сжатие
EARTH_E2 = EARTH_E ** 2
EARTH_B_A = math.sqrt(1 - EARTH_E2)  # b/a = √(1-e²)


def geodetic_to_reduced(lat_deg):

    φ_rad = math.radians(lat_deg)
    tan_φ = math.tan(φ_rad)
    tan_β = EARTH_B_A * tan_φ
    β_rad = math.atan(tan_β)
    return β_rad


def reduced_to_geodetic(beta_deg):

    β_rad = math.radians(beta_deg)
    tan_β = math.tan(β_rad)
    tan_φ = tan_β / EARTH_B_A
    φ_rad = math.atan(tan_φ)
    return math.degrees(φ_rad)


def andoyer_distance(lat1_deg, lon1_deg, lat2_deg, lon2_deg):

    u1 = geodetic_to_reduced(lat1_deg)  # u = β
    u2 = geodetic_to_reduced(lat2_deg)

    delta_lambda = math.radians(lon2_deg - lon1_deg)

    sin_u1 = math.sin(u1)
    sin_u2 = math.sin(u2)
    cos_u1 = math.cos(u1)
    cos_u2 = math.cos(u2)

    sin_dl = math.sin(delta_lambda)
    cos_dl = math.cos(delta_lambda)

    p = cos_u2 * sin_dl  # sinσ * sinA
    q = cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_dl  # sinσ * cosA
    n = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_dl  # cosσ

    sigma = math.atan2(math.sqrt(p ** 2 + q ** 2), n)

    if sigma < 1e-12:
        return 0.0

    sin_sigma = math.sin(sigma)
    cos_sigma = math.cos(sigma)

    if abs(1 + cos_sigma) < 1e-12:
        M = 0
    else:
        M = (sigma - 3 * sin_sigma) / (1 + cos_sigma)

    if abs(1 - cos_sigma) < 1e-12:
        N = 0
    else:
        N = (sigma + 3 * sin_sigma) / (1 - cos_sigma)

    U = (sin_u1 + sin_u2) ** 2
    V = (sin_u1 - sin_u2) ** 2

    alpha = EARTH_E2
    delta_sigma = -0.25 * alpha * (M * U + N * V)

    S = EARTH_A * (sigma + delta_sigma)

    return S

def andoyer_azimuth(lat1_deg, lon1_deg, lat2_deg, lon2_deg):
    B1 = math.radians(lat1_deg)
    B2 = math.radians(lat2_deg)
    delta_L = math.radians(lon2_deg - lon1_deg)

    y = math.cos(B2) * math.sin(delta_L)
    x = math.cos(B1) * math.sin(B2) - math.sin(B1) * math.cos(B2) * math.cos(delta_L)
    azimuth_rad = math.atan2(y, x)

    azimuth_deg = math.degrees(azimuth_rad)
    if azimuth_deg < 0:
        azimuth_deg += 360

    return azimuth_deg


def generate_single_daa(filename, ownship_init, ac1_init, duration_sec=10):
    lat_o = math.radians(ownship_init['lat'])
    lon_o = math.radians(ownship_init['lon'])
    alt_o = float(ownship_init['alt'])
    vx_o = ownship_init['vx']
    vy_o = ownship_init['vy']
    vz_o_func = ownship_init['vz_func']

    lat_a = math.radians(ac1_init['lat'])
    lon_a = math.radians(ac1_init['lon'])
    alt_a = float(ac1_init['alt'])
    vx_a = ac1_init['vx']
    vy_a = ac1_init['vy']
    vz_a_func = ac1_init['vz_func']

    with open(filename, 'w') as f:
        f.write("NAME, lat, lon, alt, vx, vy, vz, time\n")
        f.write("unitless, [deg], [deg], [ft], [knot], [knot], [fpm], [s]\n")

        # итеративный расчёт координат
        for t in range(0, duration_sec + 1):
            vz_o = vz_o_func(t)
            vz_a = vz_a_func(t)
            if t > 0:
                # Ownship
                dx_nm = (vx_o * KNOT_TO_MS)
                dy_nm = (vy_o * KNOT_TO_MS)

                lat_o += dy_nm / (EARTH_A)
                lon_o += dx_nm / (EARTH_A * math.sin(lat_o))
                alt_o += vz_o / 60.0

                # AC1
                dx_nm_a = vx_a * KNOT_TO_MS
                dy_nm_a = vy_a * KNOT_TO_MS
                lat_a += dy_nm_a / EARTH_A
                lon_a += dx_nm_a / (EARTH_A * math.sin(lat_a))
                alt_a += vz_a / 60.0

            lat_o_print = math.degrees(lat_o)
            lon_o_print = math.degrees(lon_o)
            lat_a_print = math.degrees(lat_a)
            lon_a_print = math.degrees(lon_a)
            f.write(
                f"Ownship, {lat_o_print:.8f}, {lon_o_print:.8f}, {alt_o:.0f}, {vx_o:.1f}, {vy_o:.1f}, {vz_o}, {t}\n")
            f.write(f"AC1, {lat_a_print:.8f}, {lon_a_print:.8f}, {alt_a:.0f}, {vx_a:.1f}, {vy_a:.1f}, {vz_a}, {t}\n")


            track_info_ownship = vx_vy_to_heading(vx_o,vy_o)
            track_info_intruder = vx_vy_to_heading(vx_a,vy_a)
            horiz_dist = andoyer_distance(lat_o_print, lon_o_print, lat_a_print, lon_a_print)
            vert_dist = abs(alt_a - alt_o)
            print(f"t={t}: S={horiz_dist:.1f}м, Δh={vert_dist:.1f}фт")

            if horiz_dist < SAFE_HORIZ_M and vert_dist < SAFE_VERT_FT:
                print(f">>> РЕАЛЬНЫЙ КОНФЛИКТ! <<<")
            else:
                print(f"Конфликта нет: горизонтальное {horiz_dist:.1f}м > {SAFE_HORIZ_M}м")
            json_template = {
                "time": t,
                "ownship": {
                    "lat": lat_o_print,
                    "lon": lon_o_print,
                    "alt": alt_o,
                    "track": track_info_ownship['track'],
                    "gs": track_info_ownship['gs'],
                    "vs": vz_o
                },
                "traffic": [{
                    "id": 'AC1',  # Идентификатор объекта
                    "lat": lat_a_print,  # Широта в градусах
                    "lon": lon_a_print,  # Долгота в градусах
                    "alt": alt_a,  # Высота
                    "track": track_info_intruder['track'],
                    "gs": track_info_intruder['gs'],
                    "vs": vz_a
                }]
            }

            # print(json_template)
            run_simulation(json_template)


def rewind_trajectory(lat_end, lon_end, alt_end, vx, vy, vz_fpm, steps):
    """
    Численно отматывает траекторию НАЗАД на `steps` секунд,
    используя ТУ ЖЕ логику, что и generate_single_daa.
    """

    alt = alt_end
    # Координаты объекта в радианах
    lat = math.radians(lat_end)
    lon = math.radians(lon_end)

    # Преобразование скоростей из узлов в м/с
    dx = vx * KNOT_TO_MS
    dy = vy * KNOT_TO_MS

    # Итеративно находим точку старта с шагом в 1 секунду
    for _ in range(steps):
        # Откатываем широту
        lat_prev = lat - dy / EARTH_A
        lon_prev = lon - dx / (EARTH_A * math.sin(lat))
        alt_prev = alt - vz_fpm / 60.0
        lat, lon, alt = lat_prev, lon_prev, alt_prev

    lat = math.degrees(lat)
    lon = math.degrees(lon)
    return lat, lon, alt


def generate_conflict_scenario(
        filename,
        duration_sec=10,
        conflict_time=None,
        region=None,
        ownship_config=None,
        ac1_config=None,
):
    if conflict_time is None:
        conflict_time = random.randint(2, duration_sec - 2)
    if not (2 <= conflict_time <= duration_sec - 2):
        raise ValueError("conflict_time must be between 2 and duration_sec-2")

    if region is None:
        region = {
            'lat_center': 33.792,
            'lat_spread': 0.01,
            'lon_center': 117.015,
            'lon_spread': 0.02,
            'alt_center': 4500,
            'alt_spread': 1000,
        }

    if ownship_config is None:
        ownship_config = {
            'speed_knot': (90, 180),
            'heading_range': (0, 359),
            'altitude_range': (4000, 6000),
            'vz_range': (-500, 500),
        }

    if ac1_config is None:
        ac1_config = {
            'speed_knot': (140, 160),
            'heading_range': (0, 359),
            'vz_range': (-500, 500),
        }

    # 1. Позиция Ownship в момент конфликта
    lat_c = random.uniform(region['lat_center'] - region['lat_spread'],
                           region['lat_center'] + region['lat_spread'])
    lon_c = random.uniform(region['lon_center'] - region['lon_spread'],
                           region['lon_center'] + region['lon_spread'])
    alt_c = ownship_config.altitude  # Генерация случайной высоты в футах

    # 2. Позиция AC1: ГАРАНТИРОВАННОЕ нарушение
    # Определение случайного положения нарушитля внутрри безопасного круга
    # Положение нарушителя определяется с помощьью случайного радиуса и угла
    # Горизонтальное смещение: внутри круга радиусом SAFE_HORIZ_M

    r = random.uniform(0, SAFE_HORIZ_M)
    angle = random.uniform(0, 2 * math.pi)
    dx_m = r * math.cos(angle)  # Метры
    dy_m = r * math.sin(angle)

    # Определение радиусов кривизны
    ro1 = (EARTH_A * (1 - EARTH_E ** 2)) / (math.sqrt((1 - EARTH_E ** 2 * (math.sin(lat_c)) ** 2) ** 3)) + alt_c
    ro2 = (EARTH_A) / (math.sqrt(1 - EARTH_E ** 2 * (math.sin(lat_c)) ** 2)) + alt_c

    # Перевод из метра в градусы
    dlat = dy_m / ro1  # Радианы
    dlat = math.degrees(dlat)
    ac_lat_c = lat_c + dlat  # Градусы
    dlon = dx_m / (ro2 * math.cos(math.radians(ac_lat_c)))  # Радианы
    dlon = math.degrees(dlon)
    ac_lon_c = lon_c + dlon  # Градусы

    # Вертикальное смещение: ГАРАНТИРОВАННО внутри порога
    alt_offset = random.uniform(-SAFE_VERT_FT + 1, SAFE_VERT_FT - 1)
    ac_alt_c = alt_c + alt_offset

    # 3. Движение
    own_heading = ownship_config.heading_deg
    own_speed = ownship_config.speed_knot
    own_vx, own_vy = heading_to_vx_vy(own_heading, own_speed)  # Проекции скоростей в узлах
    own_vz = ownship_config.vz  # Вертикальная скорость в футах

    ac_speed = ac1_config.speed_knot
    ac_heading = ac1_config.heading_deg
    ac_vx, ac_vy = heading_to_vx_vy(ac_heading, ac_speed)
    ac_vz = ac1_config.vz

    #  4. Численная отмотка назад
    lat_o0, lon_o0, alt_o0 = rewind_trajectory(lat_c, lon_c, alt_c, own_vx, own_vy, own_vz, conflict_time)
    lat_a0, lon_a0, alt_a0 = rewind_trajectory(ac_lat_c, ac_lon_c, ac_alt_c, ac_vx, ac_vy, ac_vz, conflict_time)

    # 5. Функции vz (постоянные)
    def make_vz_func(vz_val):
        return lambda t: int(round(vz_val))

    ownship_init = {
        'lat': lat_o0, 'lon': lon_o0, 'alt': alt_o0,
        'vx': own_vx, 'vy': own_vy, 'vz_func': make_vz_func(own_vz)
    }
    ac1_init = {
        'lat': lat_a0, 'lon': lon_a0, 'alt': alt_a0,
        'vx': ac_vx, 'vy': ac_vy, 'vz_func': make_vz_func(ac_vz)
    }

    generate_single_daa(filename, ownship_init, ac1_init, duration_sec)


def generate_multiple_scenarios(
        output_dir="conflict_scenarios",
        num_scenarios=20,
        duration_sec=10,
        conflict_time=None,
        seed=42,
        region=None,
        ownship_config=None,
        ac1_config=None
):
    if seed is not None:
        random.seed(seed)
    os.makedirs(output_dir, exist_ok=True)

    for i in range(1, num_scenarios + 1):
        fname = os.path.join(output_dir, f"conflict_{i:04d}.daa")
        generate_conflict_scenario(
            filename=fname,
            duration_sec=duration_sec,
            conflict_time=conflict_time,
            region=region,
            ownship_config=ownship_config,
            ac1_config=ac1_config
        )
        print(f"Создан: {fname}")


#
#
#
if __name__ == "__main__":
    ownship_params = Ownship()

    ac1_params = Intruder()

    region = {
        'lat_center': 45,
        'lat_spread': 0.01,
        'lon_center': 45,
        'lon_spread': 0.02,
        'alt_center': 5000,
        'alt_spread': 1000,
    }
    # Параметры конфликта
    generate_multiple_scenarios(
        output_dir="random_conflicts",
        num_scenarios=50000,  # Количество генерируемых сценариев
        duration_sec=120,  # Продолжительность полёта
        conflict_time=100,  # Примерное время конфликта
        seed=2025,
        region=region,
        ownship_config=ownship_params,
        ac1_config=ac1_params
    )