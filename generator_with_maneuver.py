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

# ============================================================
# 1. ПАРАМЕТРЫ WGS-84
# ============================================================
EARTH_A = 6378137.0  # большая полуось (м)
EARTH_F = 1 / 298.257223563  # сжатие
EARTH_B = EARTH_A * (1 - EARTH_F)  # малая полуось (м)
EARTH_E2 = 2 * EARTH_F - EARTH_F ** 2  # квадрат эксцентриситета
EARTH_E = math.sqrt(EARTH_E2)  # эксцентриситет
EARTH_B_A = math.sqrt(1 - EARTH_E2)  # b/a

# Константы безопасности
type_conflict = "NMAC"  # или "LOWC"

if type_conflict == "LOWC":
    SAFE_HORIZ_M = 609.9
    SAFE_VERT_FT = 500.0
elif type_conflict == "NMAC":
    SAFE_HORIZ_M = 150.0
    SAFE_VERT_FT = 100.0
else:
    raise ValueError("type_conflict must be 'NMAC' or 'LOWC'")

KNOT_TO_MS = 0.514444


# ============================================================
# 2. ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
# ============================================================

def heading_to_vx_vy(heading_deg, speed_knot):
    """Преобразует курс и скорость в компоненты vx, vy (узлы)"""
    rad = math.radians(heading_deg)
    vx = speed_knot * math.sin(rad)
    vy = speed_knot * math.cos(rad)
    return vx, vy


def vx_vy_to_heading(vx, vy):
    """Преобразует компоненты vx, vy в курс и путевую скорость"""
    gs = math.sqrt(vx ** 2 + vy ** 2)
    heading_rad = math.atan2(vx, vy)
    heading_deg = math.degrees(heading_rad)
    if heading_deg < 0:
        heading_deg += 360
    return {'track': heading_deg, 'gs': gs}


# ============================================================
# 3. МЕТОД АНДУАЙЕ (обратная геодезическая задача)
# ============================================================

def geodetic_to_reduced(lat_deg):
    """Геодезическая → приведенная широта"""
    φ_rad = math.radians(lat_deg)
    tan_β = EARTH_B_A * math.tan(φ_rad)
    return math.atan(tan_β)


def andoyer_distance(lat1_deg, lon1_deg, lat2_deg, lon2_deg):
    """Расстояние между точками по поверхности эллипсоида (метры)"""
    u1 = geodetic_to_reduced(lat1_deg)
    u2 = geodetic_to_reduced(lat2_deg)
    delta_lambda = math.radians(lon2_deg - lon1_deg)

    sin_u1, sin_u2 = math.sin(u1), math.sin(u2)
    cos_u1, cos_u2 = math.cos(u1), math.cos(u2)
    sin_dl, cos_dl = math.sin(delta_lambda), math.cos(delta_lambda)

    p = cos_u2 * sin_dl
    q = cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_dl
    n = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_dl

    sigma = math.atan2(math.sqrt(p ** 2 + q ** 2), n)

    if sigma < 1e-12:
        return 0.0

    sin_sigma, cos_sigma = math.sin(sigma), math.cos(sigma)

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

    delta_sigma = -0.25 * EARTH_E2 * (M * U + N * V)
    return EARTH_A * (sigma + delta_sigma)


def andoyer_azimuth(lat1_deg, lon1_deg, lat2_deg, lon2_deg):
    """Азимут из точки 1 в точку 2 (градусы от севера)"""
    B1 = math.radians(lat1_deg)
    B2 = math.radians(lat2_deg)
    delta_L = math.radians(lon2_deg - lon1_deg)

    y = math.cos(B2) * math.sin(delta_L)
    x = math.cos(B1) * math.sin(B2) - math.sin(B1) * math.cos(B2) * math.cos(delta_L)
    azimuth_rad = math.atan2(y, x)

    azimuth_deg = math.degrees(azimuth_rad)
    return azimuth_deg if azimuth_deg >= 0 else azimuth_deg + 360


# ============================================================
# 4. ПРЯМАЯ ГЕОДЕЗИЧЕСКАЯ ЗАДАЧА (алгоритм Винсента)
# ============================================================

def direct_geodetic_problem(lat1_deg, lon1_deg, azimuth_deg, distance_m):
    """
    Прямая геодезическая задача.
    По точке, азимуту и расстоянию → конечная точка.
    """
    lat1 = math.radians(lat1_deg)
    lon1 = math.radians(lon1_deg)
    alpha1 = math.radians(azimuth_deg)

    # Приведенная широта
    U1 = math.atan((1 - EARTH_F) * math.tan(lat1))
    sin_U1 = math.sin(U1)
    cos_U1 = math.cos(U1)

    sigma1 = math.atan2(math.tan(U1), math.cos(alpha1))
    sin_alpha = cos_U1 * math.sin(alpha1)
    cos2_alpha = 1 - sin_alpha ** 2
    u2 = cos2_alpha * (EARTH_A ** 2 - EARTH_B ** 2) / (EARTH_B ** 2)

    # Коэффициенты
    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))

    sigma = distance_m / (EARTH_B * A)

    for _ in range(100):
        cos_2sigma_m = math.cos(2 * sigma1 + sigma)
        sin_sigma = math.sin(sigma)
        cos_sigma = math.cos(sigma)

        delta_sigma = B * sin_sigma * (
                cos_2sigma_m + B / 4 * (
                cos_sigma * (-1 + 2 * cos_2sigma_m ** 2) -
                B / 6 * cos_2sigma_m * (-3 + 4 * sin_sigma ** 2) * (-3 + 4 * cos_2sigma_m ** 2)
        )
        )

        sigma_new = distance_m / (EARTH_B * A) + delta_sigma

        if abs(sigma_new - sigma) < 1e-12:
            sigma = sigma_new
            break
        sigma = sigma_new

    sin_sigma = math.sin(sigma)
    cos_sigma = math.cos(sigma)

    lat2 = math.atan2(
        sin_U1 * cos_sigma + cos_U1 * sin_sigma * math.cos(alpha1),
        (1 - EARTH_F) * math.sqrt(sin_alpha ** 2 + (sin_U1 * sin_sigma - cos_U1 * cos_sigma * math.cos(alpha1)) ** 2)
    )

    lambda_rad = math.atan2(
        sin_sigma * math.sin(alpha1),
        cos_U1 * cos_sigma - sin_U1 * sin_sigma * math.cos(alpha1)
    )

    C = EARTH_F / 16 * cos2_alpha * (4 + EARTH_F * (4 - 3 * cos2_alpha))
    L = lambda_rad - (1 - C) * EARTH_F * sin_alpha * (
            sigma + C * sin_sigma * (
            cos_2sigma_m + C * cos_sigma * (-1 + 2 * cos_2sigma_m ** 2)
    )
    )

    lon2 = lon1 + L
    return math.degrees(lat2), math.degrees(lon2)


# ============================================================
# 5. ОТМОТКА ТРАЕКТОРИИ
# ============================================================

def rewind_trajectory(lat_end_deg, lon_end_deg, alt_end_ft,
                      heading_deg, speed_knot, vz_fpm, steps):

    lat = lat_end_deg
    lon = lon_end_deg
    alt = alt_end_ft

    current_heading = heading_deg
    distance_per_step_m = speed_knot * KNOT_TO_MS
    alt_per_step_ft = vz_fpm / 60.0

    for _ in range(steps):
        back_azimuth = (current_heading + 180) % 360
        print(current_heading,back_azimuth)
        prev_lat, prev_lon = direct_geodetic_problem(lat, lon, back_azimuth, distance_per_step_m)
        prev_alt = alt - alt_per_step_ft


        new_heading = andoyer_azimuth(prev_lat, prev_lon, lat, lon)

        lat, lon, alt = prev_lat, prev_lon, prev_alt
        current_heading = new_heading

    return lat, lon, alt


def apply_motion(lat_deg, lon_deg, alt_ft, heading_deg, speed_knot, vz_fpm, dt=1.0):
    """
    Прямое движение ВС за время dt.
    Возвращает новые координаты.
    """
    distance_m = speed_knot * KNOT_TO_MS * dt
    new_lat, new_lon = direct_geodetic_problem(lat_deg, lon_deg, heading_deg, distance_m)
    new_alt = alt_ft + vz_fpm / 60.0 * dt
    return new_lat, new_lon, new_alt


def get_safe_heading(current_heading, horizontal_bands):


    safe_sectors = []
    for band in horizontal_bands:
        if band.get('Bands_Type') == 'RECOVERY':
            safe_sectors.append({
                'low': band['low'],
                'high': band['high'],
                'center': (band['low'] + band['high']) / 2
            })

    if safe_sectors:
        best_sector = None
        min_distance = 360

        for sector in safe_sectors:
            dist = min(abs(current_heading - sector['center']),
                       360 - abs(current_heading - sector['center']))
            if dist < min_distance:
                min_distance = dist
                best_sector = sector

        if best_sector:
            min_safe_angle = best_sector['low']
            return min_safe_angle, f"RECOVERY сектор {best_sector['low']:.0f}°-{best_sector['high']:.0f}°"

    none_sectors = []
    for band in horizontal_bands:
        if band.get('Bands_Type') == 'NONE':
            none_sectors.append({
                'low': band['low'],
                'high': band['high'],
                'center': (band['low'] + band['high']) / 2
            })

    if none_sectors:
        best_sector = None
        min_distance = 360
        for sector in none_sectors:
            dist = min(abs(current_heading - sector['center']),
                       360 - abs(current_heading - sector['center']))
            if dist < min_distance:
                min_distance = dist
                best_sector = sector

        if best_sector:
            return best_sector['low'], f"NONE сектор {best_sector['low']:.0f}°-{best_sector['high']:.0f}°"

    return (current_heading + 90) % 360, "Нет безопасных секторов, разворот на 90°"

def generate_single_daa(filename, ownship_init, ac1_init, duration_sec=10):
    # Начальные координаты
    set_scenario(filename)

    lat_o = ownship_init['lat']
    lon_o = ownship_init['lon']
    alt_o = ownship_init['alt']
    vx_o = ownship_init['vx']
    vy_o = ownship_init['vy']
    vz_o_func = ownship_init['vz_func']

    lat_a = ac1_init['lat']
    lon_a = ac1_init['lon']
    alt_a = ac1_init['alt']
    vx_a = ac1_init['vx']
    vy_a = ac1_init['vy']
    vz_a_func = ac1_init['vz_func']

    # Начальные курсы
    ownship_heading = vx_vy_to_heading(vx_o, vy_o)['track']
    ownship_speed = vx_vy_to_heading(vx_o, vy_o)['gs']
    intruder_heading = vx_vy_to_heading(vx_a, vy_a)['track']
    intruder_speed = vx_vy_to_heading(vx_a, vy_a)['gs']

    with open(filename, 'w') as f:
        f.write("NAME, lat, lon, alt, vx, vy, vz, time\n")
        f.write("unitless, [deg], [deg], [ft], [knot], [knot], [fpm], [s]\n")

        for t in range(0, duration_sec + 1):
            vz_o = vz_o_func(t)
            vz_a = vz_a_func(t)

            # Запись текущих координат
            f.write(f"Ownship, {lat_o:.8f}, {lon_o:.8f}, {alt_o:.0f}, "
                    f"{vx_o:.1f}, {vy_o:.1f}, {vz_o}, {t}\n")
            f.write(f"AC1, {lat_a:.8f}, {lon_a:.8f}, {alt_a:.0f}, "
                    f"{vx_a:.1f}, {vy_a:.1f}, {vz_a}, {t}\n")

            # Проверка конфликта
            horiz_dist = andoyer_distance(lat_o, lon_o, lat_a, lon_a)
            vert_dist = abs(alt_a - alt_o)
            print(f"t={t}: S={horiz_dist:.1f}м, Δh={vert_dist:.1f}фт")

            if horiz_dist < SAFE_HORIZ_M and vert_dist < SAFE_VERT_FT:
                print(f">>> РЕАЛЬНЫЙ КОНФЛИКТ! <<<")
            else:
                print(f"Конфликта нет: горизонтальное {horiz_dist:.1f}м > {SAFE_HORIZ_M}м")
            noise_std = 1.0
            NOISE_STD = {
                'lat_m': 5.0 * noise_std,
                'lon_m': 5.0 * noise_std,
                'alt_ft': 10.0 * noise_std,
                'track_deg': 2.0 * noise_std,
                'gs_knot': 2.0 * noise_std,
                'vs_fpm': 5.0 * noise_std
            }


            # Перевод СКО из метров в градусы
            lat_rad = math.radians(lat_o)
            std_lat_deg = NOISE_STD['lat_m'] / 111320.0
            std_lon_deg = NOISE_STD['lon_m'] / (111320.0 * math.cos(lat_rad))

            # # Отправка в DAIDALUS
            # json_template = {
            #     "time": t,
            #     "ownship": {
            #         "lat": lat_o, "lon": lon_o, "alt": alt_o,
            #         "track": ownship_heading, "gs": ownship_speed, "vs": vz_o
            #     },
            #     "traffic": [{
            #         "id": 'AC1',
            #         "lat": lat_a, "lon": lon_a, "alt": alt_a,
            #         "track": intruder_heading, "gs": intruder_speed, "vs": vz_a
            #     }]
            # }
            # Зашумленные данные для отправки
            json_template = {
                "time": t,
                "ownship": {
                    "lat": lat_o + random.gauss(0, std_lat_deg),
                    "lon": lon_o + random.gauss(0, std_lon_deg),
                    "alt": alt_o + random.gauss(0, NOISE_STD['alt_ft']),
                    "track": (ownship_heading + random.gauss(0, NOISE_STD['track_deg'])) % 360,
                    "gs": max(0, ownship_speed + random.gauss(0, NOISE_STD['gs_knot'])),
                    "vs": vz_o + random.gauss(0, NOISE_STD['vs_fpm'])
                },
                "traffic": [{
                    "id": 'AC1',
                    "lat": lat_a + random.gauss(0, std_lat_deg),
                    "lon": lon_a + random.gauss(0, std_lon_deg),
                    "alt": alt_a + random.gauss(0, NOISE_STD['alt_ft']),
                    "track": (intruder_heading + random.gauss(0, NOISE_STD['track_deg'])) % 360,
                    "gs": max(0, intruder_speed + random.gauss(0, NOISE_STD['gs_knot'])),
                    "vs": vz_a + random.gauss(0, NOISE_STD['vs_fpm'])
                }]
            }
            horizontal_bands = run_simulation(json_template)
            # min_safety_angle = get_safe_heading(ownship_heading,horizontal_bands)
            # for band in horizontal_bands:
            #     if band['Bands_Type'] == 'RECOVERY':
            #         min_safety_angle = min(band['high'],band['low'])
            #         print(min_safety_angle)
            # print(json_template)
            # print(horizontal_bands)
            if horizontal_bands:
                for band in horizontal_bands:
                    if band['Bands_Type'] in ['RECOVERY']:
                        low, high = band['low'], band['high']
                        if low <= ownship_heading <= high:
                            ownship_heading = ownship_heading
                        dist_to_low = abs(ownship_heading - low)
                        dist_to_low = min(dist_to_low, 360 - dist_to_low)

                        dist_to_high = abs(ownship_heading - high)
                        dist_to_high = min(dist_to_high, 360 - dist_to_high)

                        if dist_to_low <= dist_to_high:
                            ownship_heading= low + 5
                        else:
                            ownship_heading =  high - 5
            print(ownship_heading)
            if t < duration_sec:
                # 1. Вычисляем новые координаты
                distance_o = ownship_speed * KNOT_TO_MS * 1.0
                distance_a = intruder_speed * KNOT_TO_MS * 1.0

                new_lat_o, new_lon_o = direct_geodetic_problem(
                    lat_o, lon_o, ownship_heading, distance_o
                )
                new_lat_a, new_lon_a = direct_geodetic_problem(
                    lat_a, lon_a, intruder_heading, distance_a
                )
                new_alt_o = alt_o + vz_o / 60.0
                new_alt_a = alt_a + vz_a / 60.0

                # 2. Обновляем курсы (из старой точки в новую)
                #    Это критически важно для точного движения по геодезической!
                new_ownship_heading = andoyer_azimuth(lat_o, lon_o, new_lat_o, new_lon_o)
                new_intruder_heading = andoyer_azimuth(lat_a, lon_a, new_lat_a, new_lon_a)

                # 3. Присваиваем новые значения
                lat_o, lon_o, alt_o = new_lat_o, new_lon_o, new_alt_o
                lat_a, lon_a, alt_a = new_lat_a, new_lon_a, new_alt_a
                ownship_heading = new_ownship_heading
                intruder_heading = new_intruder_heading

                # Обновляем vx, vy для записи в следующий шаг (опционально)
                ownship_speed = vx_vy_to_heading(vx_o, vy_o)['gs']
                intruder_speed = vx_vy_to_heading(vx_a, vy_a)['gs']
    finalize_scenario()


# ============================================================
# 7. ГЕНЕРАЦИЯ СЦЕНАРИЯ КОНФЛИКТА
# ============================================================

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

    # 1. Позиция Ownship в момент конфликта
    lat_c = random.uniform(region['lat_center'] - region['lat_spread'],
                           region['lat_center'] + region['lat_spread'])
    lon_c = random.uniform(region['lon_center'] - region['lon_spread'],
                           region['lon_center'] + region['lon_spread'])
    alt_c = ownship_config.altitude

    # 2. Позиция AC1 внутри безопасного круга
    r = random.uniform(0, SAFE_HORIZ_M)
    angle = random.uniform(0, 2 * math.pi)
    dx_m = r * math.cos(angle)
    dy_m = r * math.sin(angle)

    # Радиусы кривизны (с учетом высоты)
    lat_c_rad = math.radians(lat_c)
    sin_lat = math.sin(lat_c_rad)
    ro1 = (EARTH_A * (1 - EARTH_E2)) / (math.sqrt((1 - EARTH_E2 * sin_lat ** 2) ** 3)) + alt_c
    ro2 = EARTH_A / math.sqrt(1 - EARTH_E2 * sin_lat ** 2) + alt_c

    # Метры → градусы
    dlat = math.degrees(dy_m / ro1)
    ac_lat_c = lat_c + dlat
    dlon = dx_m / (ro2 * math.cos(math.radians(ac_lat_c)))
    dlon = math.degrees(dlon)
    ac_lon_c = lon_c + dlon

    # Вертикальное смещение (гарантированно внутри порога)
    # alt_offset = random.uniform(-SAFE_VERT_FT + 1, SAFE_VERT_FT - 1)
    alt_offset = 0
    ac_alt_c = alt_c + alt_offset

    # 3. Параметры движения
    own_heading = ownship_config.heading_deg
    own_speed = ownship_config.speed_knot
    own_vz = ownship_config.vz
    # Для записи в .daa нужны vx, vy
    own_vx, own_vy = heading_to_vx_vy(own_heading, own_speed)

    ac_heading = ac1_config.heading_deg
    ac_speed = ac1_config.speed_knot
    ac_vz = ac1_config.vz
    ac_vx, ac_vy = heading_to_vx_vy(ac_heading, ac_speed)

    # 4. ОТМОТКА НАЗАД (используем heading и speed)
    lat_o0, lon_o0, alt_o0 = rewind_trajectory(
        lat_c, lon_c, alt_c,
        own_heading, own_speed, own_vz,
        conflict_time
    )
    lat_a0, lon_a0, alt_a0 = rewind_trajectory(
        ac_lat_c, ac_lon_c, ac_alt_c,
        ac_heading, ac_speed, ac_vz,
        conflict_time
    )

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


# ============================================================
# 8. ГЕНЕРАЦИЯ МНОЖЕСТВА СЦЕНАРИЕВ
# ============================================================

def generate_multiple_scenarios(
        output_dir="conflict_scenarios_with_maneuver",
        num_scenarios=20,
        duration_sec=10,
        conflict_time=None,
        seed=42,
        region=None,
        ownship_config=None,
        ac1_config=None
):
    # ownship_params = Ownship()  # ← новый экземпляр
    # ac1_params = Intruder()
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


# ============================================================
# 9. ТОЧКА ВХОДА
# ============================================================

if __name__ == "__main__":
    ownship_params = Ownship(speed_knot=120, altitude=5000, vz=0)
    ac1_params = Intruder(speed_knot=150, vz=0)

    region = {
        'lat_center': 45,
        'lat_spread': 0.01,
        'lon_center': 45,
        'lon_spread': 0.02,
        'alt_center': 5000,
        'alt_spread': 1000,
    }

    generate_multiple_scenarios(
        output_dir="random_conflicts_with_maneuver",
        num_scenarios=10,
        duration_sec=30,
        conflict_time=20,
        seed=2025,
        region=region,
        ownship_config=ownship_params,
        ac1_config=ac1_params
    )